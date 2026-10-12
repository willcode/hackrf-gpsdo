/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Runs the GPS app on the host for gps/test_gps_app.py, with the sample
 * buffer, the data region and the register service simulated.
 *
 * Input lines: "m <hex>" a host message; "w <hex words>" 16-bit words the
 * M0 writes to the sample buffer, followed by one tick; "t" a tick;
 * "f <path> <offset> <bytes>" that many bytes of a saved stream from a byte
 * offset, 2,048 bytes a tick (what a tick reads), then 16 idle ticks, with the
 * register writes not printed; a path of "-" goes on where the last one ended.
 * "e <path> <offset>" makes the harness an engine for snapshots: a setting-1
 * command the app writes to register 6 is answered with that many samples of
 * the 8-bit I/Q recording at path (offset samples in) and the trailer, the
 * engine's sample index advancing 100 samples a tick; "n <ticks>" runs
 * that many ticks.
 * Output lines: "r <hex>" a reply; "6 <hex>" each register-6 write; "b a b c"
 * each PLL B ratio the app asks for (always printed). A ratio asked for is
 * written, as the firmware would, after the call that asked.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "app.h"

uint8_t app_host_data[0xA000];
uint8_t app_host_samples[0x8000];
uint32_t app_host_m0_count;
static int quiet;
static FILE* stream;

/* The engine for snapshots: the recording, the sample index, and the command being read. */
static struct {
	FILE* f;
	long offset;
	uint32_t now;
	int op, nbytes, got;
	uint32_t payload;
	uint32_t pending_count; /* a snapshot asked for */
} eng;

static void put_word(uint16_t w)
{
	app_host_samples[app_host_m0_count++ % 0x8000] = w & 0xFF;
	app_host_samples[app_host_m0_count++ % 0x8000] = w >> 8;
}

/* The snapshot asked for: count samples from the engine's index now, then the trailer (the index
 * at the trailer's second word is the last sample's plus 2; gps_engine.decode_snapshot). */
static void snapshot(void)
{
	const uint32_t count = eng.pending_count, first = eng.now;
	eng.pending_count = 0;
	int8_t iq[2];
	fseek(eng.f, 2 * (eng.offset + (long) first), SEEK_SET);
	for (uint32_t k = 0; k < count; k++) {
		if (fread(iq, 1, 2, eng.f) != 2) {
			iq[0] = iq[1] = 0;
		}
		put_word((uint8_t) iq[0] | ((uint16_t) (uint8_t) iq[1] << 8));
	}
	const uint32_t second = first + count - 1 + 2;
	const uint32_t lo = second & 0xFFFF;
	put_word(0xD7 | ((count / 1024 - 1) << 8));
	put_word(lo);
	put_word(((second >> 16) + (lo == 0xFFFF)) & 0xFFFF);
	eng.now += count;
}

/* Register-6 bytes (dsp/gps_cmd.py): a header byte (bit 6) with the op, then 6 bits a byte. */
static void engine_byte(uint32_t b)
{
	static const int lengths[4] = {20, 11, 0, 2};
	if (b & 0x40) {
		eng.op = (b >> 4) & 3;
		eng.nbytes = lengths[eng.op];
		eng.got = 0;
		eng.payload = 0;
		return;
	}
	if (eng.got < eng.nbytes) {
		if (eng.got < 5) {
			eng.payload |= (b & 0x3F) << (6 * eng.got);
		}
		if ((++eng.got == eng.nbytes) && (eng.op == 3) && ((eng.payload & 15) == 1)) {
			eng.pending_count = ((eng.payload >> 4) + 1) * 1024;
		}
	}
}

uint32_t app_fpga_write(uint32_t reg, uint32_t value)
{
	if (reg != 6) {
		return 2;
	}
	if (eng.f) {
		engine_byte(value & 0xFF);
	}
	if (!quiet) {
		printf("6 %02x\n", value & 0xFF);
	}
	return 0;
}

static struct {
	int pending;
	uint32_t a, b, c, writes;
} pll_b;

uint32_t app_pll_b_set(uint32_t a, uint32_t b, uint32_t c)
{
	if (pll_b.pending) {
		return 1;
	}
	printf("b %u %u %u\n", a, b, c);
	pll_b.a = a;
	pll_b.b = b;
	pll_b.c = c;
	pll_b.pending = 1;
	return 0;
}

void app_pll_b_status(uint32_t out[4])
{
	out[0] = pll_b.writes;
	out[1] = pll_b.writes * 1000u; /* cycle counts stand in */
	out[2] = pll_b.writes * 1000u + 500u;
	out[3] = (uint32_t) pll_b.pending;
}

/* After each call, as the firmware's main loop does; with the engine on, each tick advances it and
 * delivers a snapshot asked for. */
static uint32_t call(uint32_t op, uint32_t arg)
{
	if (eng.f && (op == APP_OP_TICK)) {
		eng.now += 100; /* 25 us a tick */
		if (eng.pending_count) {
			snapshot();
		}
	}
	const uint32_t r = app_main(op, arg);
	if (pll_b.pending) {
		pll_b.pending = 0;
		pll_b.writes++;
	}
	return r;
}

int main(void)
{
	static char line[1 << 16];
	call(APP_OP_INIT, 0);
	call(APP_OP_TICK, 0); /* syncs to the (empty) buffer */
	while (fgets(line, sizeof(line), stdin)) {
		if (line[0] == 'm') {
			unsigned n = 0, v;
			for (const char* p = line + 2; sscanf(p, "%2x", &v) == 1; p += 2) {
				if (n == 64) { /* the firmware's APP_MESSAGE_MAX */
					fprintf(stderr, "a message over 64 bytes\n");
					return 1;
				}
				app_host_data[0xA000 - 256 + n++] = v;
			}
			const uint32_t r = call(APP_OP_MESSAGE, n);
			printf("r ");
			for (uint32_t i = 0; i < r; i++) {
				printf("%02x", app_host_data[0xA000 - 128 + i]);
			}
			printf("\n");
		} else if (line[0] == 'w') {
			unsigned v;
			int k;
			for (const char* p = line + 2; sscanf(p, "%x%n", &v, &k) == 1; p += k) {
				app_host_samples[app_host_m0_count % 0x8000] = v;
				app_host_samples[(app_host_m0_count + 1) % 0x8000] = v >> 8;
				app_host_m0_count += 2;
			}
			call(APP_OP_TICK, 0);
		} else if (line[0] == 't') {
			call(APP_OP_TICK, 0);
		} else if (line[0] == 'e') {
			char path[1024];
			if ((sscanf(line + 2, "%1023s %ld", path, &eng.offset) != 2) ||
			    !(eng.f = fopen(path, "rb"))) {
				return 1;
			}
		} else if (line[0] == 'n') {
			long n = atol(line + 2);
			while (n-- > 0) {
				call(APP_OP_TICK, 0);
			}
		} else if (line[0] == 'f') {
			char path[1024];
			long offset, bytes;
			if (sscanf(line + 2, "%1023s %ld %ld", path, &offset, &bytes) != 3) {
				return 1;
			}
			if (strcmp(path, "-") != 0) {
				if (stream) {
					fclose(stream);
				}
				stream = fopen(path, "rb");
				if (!stream || fseek(stream, offset, SEEK_SET) != 0) {
					return 1;
				}
			}
			quiet = 1;
			uint8_t chunk[2048];
			size_t n;
			while ((bytes > 0) &&
			       (n = fread(chunk, 1, bytes < 2048 ? (size_t) bytes : 2048, stream)) > 0) {
				for (size_t i = 0; i < n; i++) {
					app_host_samples[app_host_m0_count++ % 0x8000] = chunk[i];
				}
				bytes -= (long) n;
				call(APP_OP_TICK, 0);
			}
			/* Ticks with nothing new, so a fix in progress finishes. */
			for (int i = 0; i < 16; i++) {
				call(APP_OP_TICK, 0);
			}
			quiet = 0;
		}
	}
	return 0;
}
