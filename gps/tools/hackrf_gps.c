/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * hackrf_gps runs a HackRF Pro as a GPS frequency reference. It loads the GPS image
 * (fpga-build/top/gps_rx.py) and the GPS app (apps/gps/gps_app.c) into RAM and sets them up. It
 * then prints what the app reports once a second: the channels, the fix and the oscillator's
 * frequency offset. With -B or -D it steers PLL B, which clocks CLKOUT and the RF references. The
 * app does the search, the tracking, the fixes and the estimate on the board.
 *
 *   hackrf_gps [-o out.bin] [-t seconds] [-l lna] [-g vga] [-a amp] [-p antenna_power]
 *              [-P x,y,z | -P run.fix | -P grid[,height]] [-E height] [-B t:a,b,c;t:a,b,c...]
 *              [-D t] [-S s]
 *              [-G gps_rx.pack] [-A gps_app.bin] [-N] [-K] [-W state | -W none]
 *              [-H [-k snapshots] [-n channels] [-c first_channel]]
 *   hackrf_gps -f recording.iq [-k snapshots] [-s offset_ms]
 *
 * -h prints every option. -t sets the seconds of tracking (default 60); -t 0 runs until Ctrl-C.
 * -p 1 powers an active antenna; turn it on only when that antenna is the only thing on the port.
 * The RF amplifier is on by default (-a 1): without it the receiver's own noise swamps an active
 * antenna's output.
 *
 * The tool runs on the one HackRF Pro whose firmware takes apps. It leaves any other HackRF alone,
 * and it refuses to start when two such Pros are attached.
 *
 * A run, in order:
 * 1. Switches the board to GPS: stops any app, loads the image and the app over USB, and prints how
 *    long that took. -G and -A name the files. The defaults are the repository's own, found from
 *    this program's place in it. -N keeps what is loaded.
 * 2. Tunes to L1 at 4 Msps, sets the image to send records (setting 0 = 1) and starts receiving.
 * 3. Sends the app the tuner's offset and quarter shift ('F'), the height a three-satellite fix
 *    holds ('H') and the held position or its survey ('U').
 * 4. Starts the search on the box ('C'). The app takes its own snapshots, searches them and arms
 *    the satellites it finds. Its manager then keeps or stops each channel and searches again while
 *    a channel is free. The search's progress ('K') is printed once a second until a pass arms a
 *    satellite, for at most -S seconds (default 120).
 * 5. Tracks for -t seconds. Each second it prints the app's status ('Q'): every tracking channel
 *    with its PRN, a search under way, and the M4's load. The load compares the main loop's passes
 *    a second with the count before any channel was armed. A line follows for each channel armed,
 *    re-armed or stopped. The fix ('P') and the frequency estimate ('R') follow, or the reason
 *    there is none.
 * 6. Ends at -t, Ctrl-C, SIGTERM or SIGHUP alike. It saves the warm start, turns the steering loop
 *    off if it ran, puts PLL B back to nominal if the run moved it, stops the channels, turns
 *    antenna power off and stops receiving. It then stops the app and loads stock image 0 from
 *    flash, and the stock tools work at once. -K keeps GPS loaded.
 *    A process killed outright leaves the radio receiving with antenna power on, until the next
 *    stop or USB reconfiguration.
 *
 * The frequency estimate compares each satellite's carrier with the one expected at a held
 * position. With a position held, one satellite gives an estimate. The held position is one of:
 * - -P x,y,z: an ECEF position in meters.
 * - -P run.fix: the mean of the fixes in a previous run's .fix file.
 * - a survey: the mean of the first 60 kept fixes. A kept fix needs four satellites.
 * - -P grid[,height]: a Maidenhead square's center (4 to 10 characters, such as FN20xr), held
 *   until the survey's first kept fix. The search predicts from the almanac from the start.
 * - the warm start's position, held the same way as a grid square's center.
 * Until the first kept fix, a fix from three satellites holds the height and serves as the held
 * position. It replaces a grid square's center; a warm start's position stays. -E gives that
 * height in meters above the WGS 84 ellipsoid, which lies within about 100 m of sea level. The
 * default is 0, or the grid square's or the warm start's height.
 *
 * With no position at all, a search narrows itself once the almanac is read (12.5 minutes of data,
 * or the warm start's). It searches only the PRNs and Doppler bins visible from where the tracked
 * satellites' Dopplers agree.
 *
 * With -o out.bin the tool writes three files:
 * - out.bin: the stream as received, records and snapshots.
 * - out.bin.fix: a line per fix: sample index, GPS time of week, ECEF x, y, z (m), satellites and
 *   GDOP. The printed output carries no position.
 * - out.bin.freq: a line per estimate: the sample index its second of samples ended at, the
 *   estimate and the used satellites' rms about it (ppb), the channels used and usable (bit
 *   masks), the held position's state (0 none, 1 surveying, 2 surveyed, 3 given), each channel's
 *   offset (ppb, nan where not usable), each channel's symbol-sum phase lock, and the PLL B
 *   steering asked after this estimate (ppb). That steering holds over the next line's window:
 *   CLKOUT's error in a line is its estimate plus the line before's steering.
 *
 * -B t:a,b,c;... sets PLL B's feedback ratio to a + b/c at second t of tracking (nominal 32), for
 * testing steps. The firmware keeps the ratio within 5 ppm of nominal. Each write is printed from
 * the app's log ('J').
 * -D t turns the app's steering loop on at second t of tracking. PLL B then follows the TCXO's
 * predicted error each second: CLKOUT's error is the TCXO's less its prediction.
 *
 * The warm start (-W, default gps/data/hackrf_gps.state; -W none for none) keeps what the app
 * learned: its almanac, the held position once surveyed or given, and the carrier offset of a
 * tracked satellite. A run within four weeks hands it back with GPS time from this computer's
 * clock and the almanac moved to the current week. With a held position and the carrier offset,
 * the first search predicts the visible satellites at three Dopplers each, at 128 blocks.
 *
 * -H runs the search on the host instead of the box, as a reference. It takes K snapshots of 32 ms
 * (-k, default 6) 100 ms apart and searches them as gps_model.acquire_all does: +-6 kHz in 250 Hz
 * bins, 1 ms blocks summed noncoherently, with its threshold and cross-correlation rule. It
 * searches each detection again in K fresh snapshots, within 250 Hz in 25 Hz steps. It arms up to
 * seven channels at a code epoch about a second ahead. -n fills up to n channels with the
 * next-best PRNs, as a load test; -c names the first channel. This search takes seconds of CPU per
 * snapshot.
 *
 * The second form runs that search alone on a recording of 8-bit I/Q at 4 Msps, as
 * hackrf_transfer writes it. It takes K slices 100 ms apart from -s milliseconds in (default 100)
 * and prints the arms a live run would send. It opens no device.
 */

#include <complex.h>
#include <ctype.h>
#include <fcntl.h>
#include <fftw3.h>
#include <hackrf.h>
#include <libusb.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define FS             4000000.0
#define N              4000 /* samples per code period, per millisecond */
#define L1_HZ          1575.42e6
#define CHIP_RATE      1.023e6
#define GPS_WEEK       604800.0
#define SNAP_VALUE     127 /* (127 + 1) * 1024 samples, 32.8 ms */
#define SNAP_LEN       ((SNAP_VALUE + 1) * 1024)
#define BLOCKS         32
#define MAX_SNAPS      16
#define SNAP_GAP_MS    100
#define CHANNELS       7
#define REQ_APP_MSG    80
#define REQ_APP_REPLY  81
#define REQ_BEGIN      70 /* the RAM stage: begin, write, load an image */
#define REQ_WRITE      71
#define REQ_LOAD       72
#define REQ_ABORT      73
#define REQ_APP_START  77
#define REQ_APP_STOP   78
#define REQ_APP_STATUS 79
#define REQ_BOARD_ID   14
#define REQ_FPGA_READ  50
#define REQ_RADIO_READ 60   /* value: register, index: bank; 8 bytes */
#define REG_DROPPED    0x07 /* the command port's dropped count, wrapping at 256 */

static double now_s(void);

static volatile sig_atomic_t interrupted;

static void on_signal(int sig)
{
	(void) sig;
	interrupted = 1;
}

static const int prns[] = {1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,
			   12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,
			   23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  120,
			   121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131,
			   132, 133, 134, 135, 136, 137, 138};
#define NPRN ((int) (sizeof(prns) / sizeof(prns[0])))

/* gps_model.G2_DELAY */
static int g2_delay(int prn)
{
	static const int gps[33] = {0,   5,   6,   7,   8,   17,  18,  139, 140,
				    141, 251, 252, 254, 255, 256, 257, 258, 469,
				    470, 471, 472, 473, 474, 509, 512, 513, 514,
				    515, 516, 859, 860, 861, 862};
	static const int sbas[19] = {
		145,
		175,
		52,
		21,
		237,
		235,
		886,
		657,
		634,
		762,
		355,
		1012,
		176,
		603,
		130,
		359,
		595,
		68,
		386};
	return prn <= 32 ? gps[prn] : sbas[prn - 120];
}

/* gps_model.lfsr_out: the outputs of a 10-stage generator from all ones. */
static void lfsr(const int* taps, int ntaps, uint8_t* out)
{
	int reg[10];
	for (int i = 0; i < 10; i++) {
		reg[i] = 1;
	}
	for (int k = 0; k < 1023; k++) {
		out[k] = reg[9];
		int f = 0;
		for (int t = 0; t < ntaps; t++) {
			f ^= reg[taps[t] - 1];
		}
		memmove(&reg[1], &reg[0], 9 * sizeof(int));
		reg[0] = f;
	}
}

static const int G1_TAPS[] = {3, 10};
static const int G2_FEEDBACK[] = {2, 3, 6, 8, 9, 10};

/* The code (0/1) and the arm's generator values, as gps_track_model.arm_values
 * gives them. */
static void code_and_arm(
	int prn,
	uint8_t* code,
	uint32_t* g1v,
	uint32_t* g2v,
	uint32_t* win)
{
	uint8_t g1o[1023], g2o[1023], g2d[1023];
	lfsr(G1_TAPS, 2, g1o);
	lfsr(G2_FEEDBACK, 6, g2o);
	const int d = g2_delay(prn);
	for (int k = 0; k < 1023; k++) {
		g2d[k] = g2o[(k - d + 1023 * 2) % 1023];
		code[k] = g1o[k] ^ g2d[k];
	}
	/* g2_start: the first ten delayed G2 chips, the first one bit 9 (stage 10) */
	int g1[10], g2[10];
	for (int i = 0; i < 10; i++) {
		g1[i] = 1;
		g2[i] = g2d[9 - i];
	}
	for (int s = 0; s < 2; s++) {
		const int f1 = g1[2] ^ g1[9];
		const int f2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9];
		memmove(&g1[1], &g1[0], 9 * sizeof(int));
		memmove(&g2[1], &g2[0], 9 * sizeof(int));
		g1[0] = f1;
		g2[0] = f2;
	}
	*g1v = *g2v = 0;
	for (int i = 0; i < 10; i++) {
		*g1v |= (uint32_t) g1[i] << i;
		*g2v |= (uint32_t) g2[i] << i;
	}
	*win = code[1022] | (code[0] << 1) | (code[1] << 2);
}

/* ---- the stream */

static FILE* out_file;
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static uint16_t* capture; /* words since the snapshot was asked for */
static size_t capture_n, capture_max;
static volatile int capturing;
static uint64_t stream_bytes;
static int stream_failed; /* a write to the -o file failed: the recording is short */

static int rx_callback(hackrf_transfer* t)
{
	const size_t w = fwrite(t->buffer, 1, t->valid_length, out_file);
	pthread_mutex_lock(&lock);
	stream_failed |= w != (size_t) t->valid_length;
	stream_bytes += w;
	if (capturing) {
		for (int i = 0; i + 1 < t->valid_length && capture_n < capture_max; i +=
		     2) {
			capture[capture_n++] = t->buffer[i] | (t->buffer[i + 1] << 8);
		}
	}
	pthread_mutex_unlock(&lock);
	return 0;
}

/* ---- messages to the app */

static libusb_device_handle* usb;

static int app_message(const uint8_t* msg, int len, uint8_t* reply, int max)
{
	const uint8_t out = LIBUSB_ENDPOINT_OUT | LIBUSB_REQUEST_TYPE_VENDOR |
		LIBUSB_RECIPIENT_DEVICE;
	const uint8_t in =
		LIBUSB_ENDPOINT_IN | LIBUSB_REQUEST_TYPE_VENDOR | LIBUSB_RECIPIENT_DEVICE;
	if (libusb_control_transfer(
		    usb,
		    out,
		    REQ_APP_MSG,
		    0,
		    0,
		    (uint8_t*) msg,
		    len,
		    1000) != len) {
		return -1;
	}
	for (int i = 0; i < 200; i++) {
		usleep(5000);
		const int r = libusb_control_transfer(
			usb,
			in,
			REQ_APP_REPLY,
			0,
			0,
			reply,
			max,
			1000);
		if (r != 0) {
			return r;
		}
	}
	return -1;
}

/* ---- the switch to GPS: the image and the app loaded from files over USB */

static int ctrl(int in, uint8_t req, uint32_t arg, uint8_t* data, uint16_t len, unsigned ms)
{
	const uint8_t type = (in ? LIBUSB_ENDPOINT_IN : LIBUSB_ENDPOINT_OUT) |
		LIBUSB_REQUEST_TYPE_VENDOR | LIBUSB_RECIPIENT_DEVICE;
	return libusb_control_transfer(usb, type, req, arg >> 16, arg & 0xFFFF, data, len, ms);
}

static uint32_t crc32(const uint8_t* data, size_t len)
{
	uint32_t r = 0xffffffff;
	for (size_t i = 0; i < len; i++) {
		r ^= data[i];
		for (int b = 0; b < 8; b++) {
			r = (r >> 1) ^ ((r & 1) ? 0xedb88320 : 0);
		}
	}
	return r ^ 0xffffffff;
}

/* A file into the stage (BEGIN, then WRITE in 4 KB pieces); its CRC-32, or 0 on failure. */
static uint32_t stage_file(const char* path)
{
	FILE* f = fopen(path, "rb");
	static uint8_t buf[256 * 1024];
	const size_t len = f ? fread(buf, 1, sizeof(buf), f) : 0;
	const int more = f && (fgetc(f) != EOF);
	if (f) {
		fclose(f);
	}
	if (!len || more) {
		fprintf(stderr, "%s: %s\n", path, more ? "larger than 256 KB" : "cannot read");
		return 0;
	}
	if (ctrl(0, REQ_BEGIN, (uint32_t) len, NULL, 0, 1000) != 0) {
		fprintf(stderr, "%s: the stage refused it (receiving, or an app still running)\n", path);
		return 0;
	}
	for (size_t off = 0; off < len; off += 4096) {
		const uint16_t n = (uint16_t) ((len - off < 4096) ? len - off : 4096);
		if (ctrl(0, REQ_WRITE, (uint32_t) off, &buf[off], n, 1000) != n) {
			fprintf(stderr, "%s: stage write at %zu failed\n", path, off);
			ctrl(0, REQ_ABORT, 0, NULL, 0, 1000);
			return 0;
		}
	}
	return crc32(buf, len);
}

/* Any app stopped: the firmware ends it at its next poll. */
static void app_stop_wait(void)
{
	ctrl(0, REQ_APP_STOP, 0, NULL, 0, 1000); /* refused when none runs */
	for (int i = 0; i < 100; i++) {
		uint8_t st[44];
		if ((ctrl(1, REQ_APP_STATUS, 0, st, sizeof(st), 1000) == (int) sizeof(st)) && (st[0] != 1)) {
			break;
		}
		usleep(10000);
	}
}

/*
 * Whatever runs, an app or another image, gives way to GPS: the app stopped, the GPS image
 * loaded, the GPS app started. Prints each step's time.
 */
static int switch_to_gps(const char* pack, const char* app)
{
	const double t0 = now_s();
	app_stop_wait();
	const double t1 = now_s();
	uint32_t crc = stage_file(pack);
	uint8_t result = 0;
	if (!crc || (ctrl(1, REQ_LOAD, crc, &result, 1, 3000) != 1) || (result != 1)) {
		fprintf(stderr, "%s: the image did not load (result %u)\n", pack, result);
		ctrl(0, REQ_ABORT, 0, NULL, 0, 1000);
		return 1;
	}
	const double t2 = now_s();
	crc = stage_file(app);
	uint8_t cause = 0xFF;
	if (!crc || (ctrl(1, REQ_APP_START, crc, &cause, 1, 3000) != 1) || (cause != 0)) {
		fprintf(stderr, "%s: the app did not start (cause %u)\n", app, cause);
		ctrl(0, REQ_ABORT, 0, NULL, 0, 1000);
		return 1;
	}
	const double t3 = now_s();
	printf("switched to GPS in %.0f ms: app stopped %.0f ms, image %.0f ms, app %.0f ms\n",
	       (t3 - t0) * 1e3,
	       (t1 - t0) * 1e3,
	       (t2 - t1) * 1e3,
	       (t3 - t2) * 1e3);
	return 0;
}

/* An FPGA register, or -1: firmware before the wider read stalls above register 6. */
static int fpga_read(uint8_t reg)
{
	const uint8_t in =
		LIBUSB_ENDPOINT_IN | LIBUSB_REQUEST_TYPE_VENDOR | LIBUSB_RECIPIENT_DEVICE;
	uint8_t v;
	if (libusb_control_transfer(usb, in, REQ_FPGA_READ, 0, reg, &v, 1, 1000) != 1) {
		return -1;
	}
	return v;
}

/* A radio register's applied value (bank 0); 0 with *ok cleared when unreadable. */
static uint64_t radio_read(int reg, int* ok)
{
	const uint8_t in =
		LIBUSB_ENDPOINT_IN | LIBUSB_REQUEST_TYPE_VENDOR | LIBUSB_RECIPIENT_DEVICE;
	uint8_t b[8];
	if (libusb_control_transfer(usb, in, REQ_RADIO_READ, reg, 0, b, 8, 1000) != 8) {
		*ok = 0;
		return 0;
	}
	uint64_t v = 0;
	for (int i = 7; i >= 0; i--) {
		v = (v << 8) | b[i];
	}
	return v;
}

/*
 * The tuner's offset from 1575.42 MHz in hertz, from the RF frequency the firmware applied
 * (register 1 of the applied bank, in 2^-24 Hz); 0 with the reason printed if it cannot be read.
 */
static double tuner_offset(void)
{
	int ok = 1;
	const uint64_t v = radio_read(1, &ok); /* the digital center, Q40.24 */
	if (!ok) {
		printf("tuner offset: not readable, taken as 0\n");
		return 0;
	}
	return (double) (int64_t) (v - (uint64_t) 1575420000 * (1u << 24)) / (1u << 24);
}

/*
 * The quarter shift in hertz: the digital center less the analog (radio.c,
 * compute_offset and digital_from_analog_rf): the AFE rate (the sample rate,
 * Q28.36, times 2^n for the receive decimation n) times the rotation's top
 * byte / 256, positive for a rotation under half a turn. The analog LO is on
 * PLL B, the shift on PLL A.
 */
static int32_t quarter_shift(void)
{
	int ok = 1;
	const uint64_t rate = radio_read(6, &ok), n = radio_read(8, &ok);
	uint64_t rot = radio_read(5, &ok) & 0xFFFFFFFFu;
	if (!ok || (rot == 0) || (n > 12)) {
		return 0;
	}
	const int negative = rot > 0x80000000u;
	if (rot & 0x80000000u) {
		rot = 0x100000000u - rot;
	}
	const double afe = (double) rate / 68719476736.0 * (double) (1u << n); /* 2^36 */
	const double q = afe * (double) (rot >> 24) / 256.0;
	return (int32_t) llround(negative ? -q : q);
}

static int app_simple(const uint8_t* msg, int len, const char* what)
{
	uint8_t r[128];
	const int n = app_message(msg, len, r, sizeof(r));
	if ((n != 1) || (r[0] != 0)) {
		fprintf(stderr,
			"%s: app answered %d (status %d)\n",
			what,
			n,
			n > 0 ? r[0] : -1);
		return -1;
	}
	return 0;
}

/* Stops every channel: an engine channel armed by an earlier run tracks on without its loops. */
static void stop_all(void)
{
	for (int c = 0; c < CHANNELS; c++) {
		const uint8_t stop[2] = {'X', c};
		app_simple(stop, 2, "stop");
	}
}

static void put32(uint8_t* p, uint32_t v)
{
	p[0] = v;
	p[1] = v >> 8;
	p[2] = v >> 16;
	p[3] = v >> 24;
}

static uint32_t get32(const uint8_t* p)
{ return p[0] | (p[1] << 8) | (p[2] << 16) | ((uint32_t) p[3] << 24); }

static double get_double(const uint8_t* p)
{
	double v;
	memcpy(&v, p, sizeof(v));
	return v;
}

/* The app's fix ('P'), printed and written when it is new. */
static struct {
	uint32_t count, n;
	double t, x[3], x0[3];
	int have, have_first;
} last_fix;

static void report_fix(FILE* f)
{
	uint8_t r[128];
	const uint8_t q = 'P';
	if (app_message(&q, 1, r, sizeof(r)) != 128) {
		return;
	}
	static const char* const results[6] = {
		"kept",
		"fewer than four usable satellites",
		"not converged",
		"not near the Earth",
		"GDOP too high",
		"three satellites with the height held: the held position until a fix"};
	const uint32_t count = get32(&r[0]);
	if ((count == 0) || (count == last_fix.count)) {
		printf("  fix: none (%u usable satellites; %s; %u rejected)\n",
		       r[11],
		       r[124] < 6 ? results[r[124]] : "?",
		       get32(&r[120]));
		return;
	}
	const uint32_t n = get32(&r[4]);
	const double t = get_double(&r[12]);
	double x[3];
	for (int i = 0; i < 3; i++) {
		x[i] = get_double(&r[20 + 8 * i]);
	}
	float dop[5];
	memcpy(dop, &r[44], sizeof(dop));
	if (!last_fix.have_first) {
		memcpy(last_fix.x0, x, sizeof(x));
		last_fix.have_first = 1;
	}
	const double dx = x[0] - last_fix.x0[0], dy = x[1] - last_fix.x0[1],
		     dz = x[2] - last_fix.x0[2];
	printf("  fix %u: %u satellites, GDOP %.1f, GPS time of week %.9f s",
	       count,
	       r[8],
	       dop[0],
	       t);
	if (last_fix.have) {
		/* GPS seconds against seconds of samples: the sample clock's offset. */
		const double samples = (double) (uint32_t) (n - last_fix.n) / FS;
		double dt = t - last_fix.t;
		if (dt < -GPS_WEEK / 2) {
			dt += GPS_WEEK; /* across the week's end */
		}
		printf(", sample clock %+.3f ppm", ((samples / dt) - 1) * 1e6);
	}
	printf(", %.1f m from the first fix\n", sqrt(dx * dx + dy * dy + dz * dz));
	fprintf(f,
		"%u %.9f %.3f %.3f %.3f %u %.2f\n",
		n,
		t,
		x[0],
		x[1],
		x[2],
		r[8],
		dop[0]);
	fflush(f);
	last_fix.count = count;
	last_fix.n = n;
	last_fix.t = t;
	last_fix.have = 1;
}

/* PLL B ratios asked for at seconds of tracking (-B). */
#define MAX_STEPS 32
static struct {
	int sec;
	uint32_t a, b, c;
} steps[MAX_STEPS];
static int nsteps, steered; /* steered: the last ratio asked for is not nominal */
static int loop_at, loop_on; /* -D: the second to turn the steering loop on; on */
static uint32_t last_writes;

static int parse_steps(const char* arg)
{
	const char* p = arg;
	while (*p && (nsteps < MAX_STEPS)) {
		int n = 0;
		if ((sscanf(p,
			   "%d:%u,%u,%u%n",
			   &steps[nsteps].sec,
			   &steps[nsteps].a,
			   &steps[nsteps].b,
			   &steps[nsteps].c,
			   &n) != 4) || (steps[nsteps].sec < 1)) {
			return -1;
		}
		nsteps++;
		p += n;
		if (*p == ';') {
			p++;
		} else if (*p) {
			return -1;
		}
	}
	return *p ? -1 : 0;
}

static int ask_ratio(uint32_t a, uint32_t b, uint32_t c)
{
	uint8_t msg[13] = {'B'};
	put32(&msg[1], a);
	put32(&msg[5], b);
	put32(&msg[9], c);
	printf("  PLL B: asking for %u + %u/%u\n", a, b, c);
	const int rc = app_simple(msg, 13, "PLL B ratio");
	if (rc == 0) {
		steered = (a != 32) || (b != 0); /* a refused ratio leaves the last one in effect */
	}
	return rc;
}

/* The app's PLL B writes ('J') not yet printed. */
static void report_writes(void)
{
	uint8_t r[128];
	const uint8_t q = 'J';
	if (app_message(&q, 1, r, sizeof(r)) != 116) {
		return;
	}
	const uint32_t count = get32(&r[0]);
	if ((count != last_writes) && !get32(&r[4 + 20])) {
		return; /* the latest write's cycle counts not in yet: printed with them */
	}
	for (uint32_t i = 0; (i < 4) && (count - i > last_writes); i++) {
		const uint8_t* e = &r[4 + 28 * i];
		const uint32_t start = get32(&e[20]), end = get32(&e[24]);
		printf("  PLL B write %u: %u + %u/%u, a record at sample %u, the write %.1f us "
		       "(%.1f us after the tick that asked)\n",
		       count - i,
		       get32(&e[0]),
		       get32(&e[4]),
		       get32(&e[8]),
		       get32(&e[16]),
		       (end - start) / 204.0,
		       (start - get32(&e[12])) / 204.0);
	}
	last_writes = count;
}

/* The steering in effect, in ppb, from the app's frequency reply; 0 when unreadable. */
static float steering_now(void)
{
	uint8_t r[128];
	const uint8_t q = 'R';
	float d = 0;
	if (app_message(&q, 1, r, sizeof(r)) == 64) {
		memcpy(&d, &r[60], 4);
	}
	return d;
}

/* Back to nominal at the end of a run left steered (the loop off first), waiting up to 2 s for
 * the write. */
static void restore_ratio(void)
{
	if (loop_on) {
		const uint8_t off[2] = {'D', 0};
		app_simple(off, 2, "steering loop off");
		loop_on = 0;
		usleep(200000); /* a request the loop made may still be written */
		report_writes();
		steered = steering_now() != 0;
	}
	if (!steered) {
		return;
	}
	const uint32_t before = last_writes;
	if (ask_ratio(32, 0, 1) != 0) {
		return;
	}
	for (int i = 0; i < 20; i++) {
		usleep(100000);
		report_writes();
		if (last_writes != before) {
			return;
		}
	}
	fprintf(stderr, "PLL B: the write back to nominal was not seen\n");
}

/* The search's mode in its reply ('K'): a cold pass over the almanac's Dopplers alone is
 * "predicted". */
static const char* search_mode(const uint8_t* k)
{
	if (k[0] == 1) {
		return k[20] ? "predicted" : "cold";
	}
	return k[0] == 2 ? "refining" : k[0] == 3 ? "predicting" : k[0] == 4 ? "locating" : "done";
}

/*
 * What changed since the last second's status ('Q'): a channel armed, re-armed or stopped, with
 * its Doppler before and a second after. A re-arm empties the channel for a status or two, or
 * starts its record count again; the same PRN back within 3 s is one.
 */
static void report_channels(const uint8_t* st)
{
	static uint8_t prn[CHANNELS], away[CHANNELS];
	static uint32_t rec[CHANNELS];
	static int32_t mhz[CHANNELS];
	for (int c = 0; c < CHANNELS; c++) {
		const uint8_t* e = &st[16 * c];
		const uint32_t n = get32(&e[4]);
		const int32_t f = (int32_t) get32(&e[8]);
		if (!e[0]) {
			if (prn[c] && (++away[c] == 3)) {
				printf("  channel %d: PRN %u stopped at %+.1f Hz\n", c, prn[c], mhz[c] * 1e-3);
				prn[c] = 0;
			}
			continue;
		}
		if (e[0] == prn[c] && (away[c] || (n < rec[c]))) {
			printf("  channel %d: PRN %u re-armed, %+.1f Hz before, %+.1f Hz a second on\n",
			       c,
			       e[0],
			       mhz[c] * 1e-3,
			       f * 1e-3);
		} else if (e[0] != prn[c]) {
			if (prn[c]) {
				printf("  channel %d: PRN %u stopped at %+.1f Hz\n", c, prn[c], mhz[c] * 1e-3);
			}
			printf("  channel %d: PRN %u armed, %+.1f Hz\n", c, e[0], f * 1e-3);
		}
		prn[c] = e[0];
		away[c] = 0;
		rec[c] = n;
		mhz[c] = f;
	}
}

/* The app's frequency estimate ('R'), printed and written when it is new. */
static uint32_t last_estimates;
static float steer_in_effect; /* asked after the estimate before: it held over this one's window */

static void report_freq(FILE* f)
{
	uint8_t r[128];
	const uint8_t q = 'R';
	if (app_message(&q, 1, r, sizeof(r)) != 64) {
		return;
	}
	static const char* const results[5] = {
		"estimated",
		"no position yet",
		"no usable satellite",
		"no window yet",
		"the satellites disagree"};
	static const char* const held[4] = {"none", "surveying", "surveyed", "given"};
	const uint32_t count = get32(&r[0]);
	if ((count == 0) || (count == last_estimates)) {
		printf("  frequency: none (%s)\n", r[18] < 5 ? results[r[18]] : "?");
		return;
	}
	last_estimates = count;
	float v[2 + CHANNELS];
	memcpy(v, &r[8], 8);
	memcpy(&v[2], &r[24], 4 * CHANNELS);
	float steer;
	memcpy(&steer, &r[60], 4);
	printf("  frequency: TCXO %+.3f ppb from %d satellites (rms %.3f ppb), position %s (%u fixes)",
	       v[0],
	       __builtin_popcount(r[17]),
	       v[1],
	       r[19] < 4 ? held[r[19]] : "?",
	       get32(&r[20]));
	if ((steer != 0) || (steer_in_effect != 0)) {
		printf("; PLL B steered %+.3f ppb, CLKOUT %+.3f ppb", steer_in_effect, v[0] + steer_in_effect);
	}
	steer_in_effect = steer;
	printf("\n");
	fprintf(f, "%u %.4f %.4f %u %u %u", get32(&r[4]), v[0], v[1], r[17], r[16], r[19]);
	for (int c = 0; c < CHANNELS; c++) {
		if (r[16] & (1u << c)) {
			fprintf(f, " %.4f", v[2 + c]);
		} else {
			fprintf(f, " nan");
		}
	}
	for (int c = 0; c < CHANNELS; c++) {
		fprintf(f, " %.2f", (int8_t) r[52 + c] / 100.0);
	}
	fprintf(f, " %.4f\n", steer);
	fflush(f);
}

/* ---- acquisition */

struct found {
	int prn;
	int epoch; /* first sample of chip 0 in snapshot 0, mod N */
	double doppler;
	double metric;
};

/* A snapshot: its samples and the engine's index of the first. */
struct snap {
	const float complex* x;
	uint32_t first;
};

static fftwf_complex *fin, *fout;
static fftwf_plan fwd, inv;
static float complex* locals;  /* per PRN: the local code's conjugate spectrum */
static float complex* spectra; /* per snapshot and block: the mixed block's spectrum */

/* The local code's conjugate spectrum (gps_model.code_samples at fs, phase 0).
 */
static void make_local(int prn, float complex* dst)
{
	uint8_t code[1023];
	uint32_t a, b, c;
	code_and_arm(prn, code, &a, &b, &c);
	for (int k = 0; k < N; k++) {
		const int chip = (int) floor(k * CHIP_RATE / FS) % 1023;
		((float complex*) fin)[k] = 1.0f - 2.0f * code[chip];
	}
	fftwf_execute(fwd);
	for (int k = 0; k < N; k++) {
		dst[k] = conjf(((float complex*) fout)[k]);
	}
}

/* Mix block b of snapshot j down by fd and transform it into spectra. */
static void mix_block(
	const struct snap* s,
	int j,
	int b,
	double fd,
	fftwf_complex* in,
	fftwf_complex* out)
{
	const double complex rot = cexp(-2 * M_PI * I * fd / FS);
	/* the phase at the block's first sample, exact; then a rotator within it */
	double complex z = cexp(-2 * M_PI * I * fd * (b * N) / FS);
	for (int k = 0; k < N; k++) {
		((float complex*) in)[k] = s[j].x[b * N + k] * (float complex) z;
		z *= rot;
	}
	fftwf_execute_dft(fwd, in, out);
	memcpy(&spectra[(j * BLOCKS + b) * N], out, N * sizeof(float complex));
}

/*
 * |correlation|^2 summed over the blocks of every snapshot, for one PRN at the mixed Doppler fd,
 * indexed by the epoch in snapshot 0. Snapshot j's epoch lies (F_j - F_0) mod P samples earlier,
 * P the code period at fd: its row is shifted by that much, rounded, before it is added.
 */
static void search_row(
	const float complex* loc,
	const struct snap* s,
	int count,
	double fd,
	fftwf_complex* in,
	fftwf_complex* out,
	float* row)
{
	const double period = 1023 / (CHIP_RATE * (1 + fd / L1_HZ)) * FS;
	for (int k = 0; k < N; k++) {
		row[k] = 0;
	}
	for (int j = 0; j < count; j++) {
		const double since = (double) (uint32_t) (s[j].first - s[0].first);
		const int shift = (int) llround(fmod(since, period)) % N;
		for (int b = 0; b < BLOCKS; b++) {
			for (int k = 0; k < N; k++) {
				((float complex*) in)[k] =
					spectra[(j * BLOCKS + b) * N + k] * loc[k];
			}
			fftwf_execute_dft(inv, in, out);
			const float complex* c = (const float complex*) out;
			for (int k = 0; k < N; k++) {
				const int to = k + shift < N ? k + shift : k + shift - N;
				row[to] += crealf(c[k]) * crealf(c[k]) +
					cimagf(c[k]) * cimagf(c[k]);
			}
		}
	}
}

static double metric_of(const float* row, int epoch)
{
	const int guard = (int) ceil(FS / CHIP_RATE) + 1;
	double noise = 0;
	int count = 0;
	for (int k = 0; k < N; k++) {
		int d = abs(k - epoch);
		d = d < N - d ? d : N - d;
		if (d > guard) {
			noise += row[k];
			count++;
		}
	}
	return row[epoch] / (noise / count);
}

/*
 * The search's work for one Doppler bin, shared by THREADS threads: first the blocks' transforms
 * (block i to thread i mod THREADS), then, after all have joined, the PRNs (PRN j to thread j mod
 * THREADS; each writes only its own PRNs' best cells). Each thread has its own buffers; FFTW's
 * new-array execute is safe across threads.
 */
#define THREADS 4

struct worker {
	int t, phase;
	const struct snap* s;
	int nsnaps;
	const int* prn_index;
	int count;
	double fd;
	struct found* best;
	float* peak;
	fftwf_complex *in, *out;
	float* row;
};

static void* work(void* arg)
{
	struct worker* w = arg;
	if (w->phase == 0) {
		for (int i = w->t; i < w->nsnaps * BLOCKS; i += THREADS) {
			mix_block(w->s, i / BLOCKS, i % BLOCKS, w->fd, w->in, w->out);
		}
		return NULL;
	}
	for (int j = w->t; j < w->count; j += THREADS) {
		search_row(
			&locals[w->prn_index[j] * N],
			w->s,
			w->nsnaps,
			w->fd,
			w->in,
			w->out,
			w->row);
		int k_max = 0;
		for (int k = 1; k < N; k++) {
			if (w->row[k] > w->row[k_max]) {
				k_max = k;
			}
		}
		if (w->row[k_max] > w->peak[j]) {
			w->peak[j] = w->row[k_max];
			w->best[j].epoch = k_max;
			w->best[j].doppler = w->fd;
			w->best[j].metric = metric_of(w->row, k_max);
		}
	}
	return NULL;
}

/* gps_model.acquire for the PRNs prn_index[0..count) over the snapshots: the best cell of each
 * over the Dopplers. */
static void acquire(
	const struct snap* s,
	int nsnaps,
	const int* prn_index,
	int count,
	double lo,
	double hi,
	double step,
	struct found* best)
{
	float peak[NPRN];
	for (int j = 0; j < count; j++) {
		best[j] = (struct found){prns[prn_index[j]], 0, 0, 0};
		peak[j] = -1;
	}
	struct worker w[THREADS];
	for (int t = 0; t < THREADS; t++) {
		w[t] = (struct worker){
			t,
			0,
			s,
			nsnaps,
			prn_index,
			count,
			0,
			best,
			peak,
			NULL,
			NULL,
			NULL};
		w[t].in = fftwf_malloc(N * sizeof(fftwf_complex));
		w[t].out = fftwf_malloc(N * sizeof(fftwf_complex));
		w[t].row = malloc(N * sizeof(float));
	}
	for (double fd = lo; fd <= hi + step / 2; fd += step) {
		for (int phase = 0; phase < 2; phase++) {
			pthread_t id[THREADS];
			for (int t = 0; t < THREADS; t++) {
				w[t].phase = phase;
				w[t].fd = fd;
				pthread_create(&id[t], NULL, work, &w[t]);
			}
			for (int t = 0; t < THREADS; t++) {
				pthread_join(id[t], NULL);
			}
		}
	}
	for (int t = 0; t < THREADS; t++) {
		fftwf_free(w[t].in);
		fftwf_free(w[t].out);
		free(w[t].row);
	}
}

static int by_metric(const void* a, const void* b)
{
	const double d =
		((const struct found*) b)->metric - ((const struct found*) a)->metric;
	return (d > 0) - (d < 0);
}

/* gps_model.detect_metric's threshold for the blocks summed, and its cross-correlation rule. */
static int detected(const struct found* f, int p, int nsnaps)
{
	const double threshold = 1 + 10.0 / sqrt(BLOCKS * nsnaps);
	const double cross = 0.02; /* gps_model.CROSS_SHARE */
	return (f[p].metric >= threshold) &&
		((p == 0) || (f[p].metric - 1 >= cross * (f[0].metric - 1)));
}

static int search_and_report(const struct snap* s, int nsnaps, struct found* f)
{
	fin = fftwf_malloc(N * sizeof(fftwf_complex));
	fout = fftwf_malloc(N * sizeof(fftwf_complex));
	fwd = fftwf_plan_dft_1d(N, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE);
	inv = fftwf_plan_dft_1d(N, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE);
	locals = malloc((size_t) NPRN * N * sizeof(float complex));
	spectra = malloc((size_t) nsnaps * BLOCKS * N * sizeof(float complex));
	int all[NPRN];
	for (int p = 0; p < NPRN; p++) {
		make_local(prns[p], &locals[p * N]);
		all[p] = p;
	}

	const double t0 = now_s();
	acquire(s, nsnaps, all, NPRN, -6000, 6000, 250, f);
	printf("search: %d snapshots, %.1f s\n", nsnaps, now_s() - t0);
	qsort(f, NPRN, sizeof(f[0]), by_metric);
	int n = 0;
	for (int p = 0; p < NPRN; p++) {
		const int d = detected(f, p, nsnaps);
		if (p < 12 || d) {
			printf("PRN %3d: metric %6.2f, Doppler %+6.0f Hz, epoch %4d%s\n",
			       f[p].prn,
			       f[p].metric,
			       f[p].doppler,
			       f[p].epoch,
			       d ? "" : " (not detected)");
		}
		n += d;
	}
	return n;
}

/*
 * The arm message for a detection: its Doppler refined in 25 Hz steps, then the first code epoch
 * at least `ahead` seconds after snapshot 0's first sample, with the code and carrier words for
 * that Doppler. Returns the epoch's sample index (the arm's start is one before it); *metric is
 * the refined search's.
 */
static uint32_t plan_arm(
	const struct snap* s,
	int nsnaps,
	const struct found* coarse,
	double ahead,
	int channel,
	uint8_t* msg,
	double* metric)
{
	int idx = 0;
	while (prns[idx] != coarse->prn) {
		idx++;
	}
	struct found fine;
	acquire(s,
		nsnaps,
		&idx,
		1,
		coarse->doppler - 250,
		coarse->doppler + 250,
		25,
		&fine);
	printf("  refined: PRN %d, metric %.2f, Doppler %+.0f Hz, epoch %d\n",
	       fine.prn,
	       fine.metric,
	       fine.doppler,
	       fine.epoch);
	*metric = fine.metric;
	const double rate = CHIP_RATE * (1 + fine.doppler / L1_HZ);
	const double period = 1023 / rate * FS;
	const double m = ceil(ahead * FS / period);
	const uint32_t start1 = s[0].first + (uint32_t) llround(fine.epoch + m * period);
	uint8_t code[1023];
	uint32_t g1, g2, win;
	code_and_arm(fine.prn, code, &g1, &g2, &win);
	const uint32_t W = (uint32_t) llround(rate / FS * 4294967296.0);
	const uint32_t Wc =
		(uint32_t) (int32_t) llround(fine.doppler / FS * 4294967296.0);
	msg[0] = 'A';
	msg[1] = channel;
	msg[2] = fine.prn;
	put32(&msg[3], W);
	put32(&msg[7], Wc);
	msg[11] = g1;
	msg[12] = g1 >> 8;
	msg[13] = g2;
	msg[14] = g2 >> 8;
	msg[15] = win;
	put32(&msg[16], start1 - 1);
	put32(&msg[20], (uint32_t) (int32_t) llround(fine.doppler * 1000));
	return start1;
}

/* The main loop's rate with receiving on and no channel armed: the base for the load. */
static void measure_base_rate(double* base_rate)
{
	uint8_t r[16];
	const uint8_t l = 'L';
	if ((app_message(&l, 1, r, sizeof(r)) >= 12) && get32(&r[4])) {
		*base_rate = (double) get32(&r[0]) / get32(&r[4]);
		printf("main loop, no channel: %.0f passes a second\n", *base_rate * 204e6);
	}
}

/* ---- main */

static double now_s(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/*
 * K snapshots, the engine's first sample index of each, into s (buffers reused when s[j].x is set);
 * *t0 is when the first was asked for. The capture restarts with the requests. Returns 0, or -1.
 */
static int take_snapshots(int k, struct snap* s, double* t0)
{
	const uint8_t snap[2] = {'S', SNAP_VALUE};
	pthread_mutex_lock(&lock);
	capture_n = 0;
	capturing = 1;
	pthread_mutex_unlock(&lock);
	/* K snapshots, then one more that pushes the last one's final words through the M0's
	 * 32-byte chunks and libhackrf's 256 kB transfers. The engine keeps them apart, and the
	 * app takes a request once the snapshot before it has started. */
	*t0 = now_s();
	for (int j = 0; j <= k; j++) {
		if (j) {
			usleep(SNAP_GAP_MS * 1000);
		}
		if (app_simple(snap, 2, "snapshot")) {
			return -1;
		}
	}
	/* Each trailer normally follows its samples at once; the search skips any stray words
	 * before a snapshot. A sample word equal to the trailer is caught by the indices, which
	 * must grow by at least a snapshot each time. */
	const uint16_t trailer = 0xD7 | (SNAP_VALUE << 8);
	long at[MAX_SNAPS];
	uint32_t first[MAX_SNAPS];
	int got = 0;
	size_t from = SNAP_LEN;
	for (int tries = 0; tries < 100 && got < k && !interrupted; tries++) {
		usleep(50000);
		pthread_mutex_lock(&lock);
		const size_t have = capture_n;
		for (size_t j = from; j + 3 <= have && got < k; j++) {
			if (capture[j] != trailer) {
				continue;
			}
			const uint32_t lo = capture[j + 1], hi = capture[j + 2];
			const uint32_t last =
				(((hi - (lo == 0xFFFF)) & 0xFFFF) << 16) | lo;
			const uint32_t f0 = last - 2 - SNAP_LEN + 1;
			if (got && (uint32_t) (f0 - first[got - 1]) < SNAP_LEN) {
				continue;
			}
			at[got] = (long) j;
			first[got++] = f0;
			from = j + 3 + SNAP_LEN;
			j = from - 1;
		}
		pthread_mutex_unlock(&lock);
	}
	capturing = 0;
	if (got < k) {
		fprintf(stderr,
			"%d of %d snapshots found in the stream (%zu words)\n",
			got,
			k,
			capture_n);
		return -1;
	}
	for (int j = 0; j < k; j++) {
		const uint16_t* w = &capture[at[j] - SNAP_LEN];
		float complex* x = s[j].x ? (float complex*) s[j].x :
					    malloc(SNAP_LEN * sizeof(float complex));
		for (int i = 0; i < SNAP_LEN; i++) {
			const int8_t re = w[i] & 0xFF, im = w[i] >> 8;
			x[i] = re + I * im;
		}
		s[j] = (struct snap){x, first[j]};
		printf("snapshot %d: %d samples from engine sample %u\n",
		       j,
		       SNAP_LEN,
		       first[j]);
	}

	return 0;
}

/* The mean position of a .fix file's fixes (report_fix's lines); 0, or -1 with none. */
static int mean_fix(const char* path, double x[3])
{
	FILE* f = fopen(path, "r");
	if (!f) {
		return -1;
	}
	double sum[3] = {0, 0, 0}, v[3], t;
	unsigned n_r, sats;
	float gdop;
	int n = 0;
	while (fscanf(f, "%u %lf %lf %lf %lf %u %f", &n_r, &t, &v[0], &v[1], &v[2], &sats, &gdop) ==
	       7) {
		for (int i = 0; i < 3; i++) {
			sum[i] += v[i];
		}
		n++;
	}
	fclose(f);
	if (n == 0) {
		return -1;
	}
	for (int i = 0; i < 3; i++) {
		x[i] = sum[i] / n;
	}
	printf("held position: the mean of %d fixes from %s\n", n, path);
	return 0;
}

/*
 * A Maidenhead locator of 4 to 10 characters (an amateur's grid square, such as FN20 or FN20xr),
 * with an optional height in meters after a comma, as the ECEF position (WGS 84) of its square's
 * center: 0, or -1 if it is not one. A 6-character square is within about 4 km, about 2 ppb on
 * one satellite's carrier; the survey that follows refines it.
 */
static int grid_position(const char* arg, double x[3], double* height)
{
	char g[11];
	double h = 0;
	size_t n = strcspn(arg, ",");
	if ((n < 4) || (n > 10) || (n % 2) || (arg[n] && (sscanf(&arg[n + 1], "%lf", &h) != 1))) {
		return -1;
	}
	memcpy(g, arg, n);
	g[n] = 0;
	double lon = -180, lat = -90, dlon = 20, dlat = 10;
	for (size_t i = 0; i < n; i += 2) {
		/* Pairs: letters A-R, digits, letters A-X, digits, letters A-X. */
		const int digits = (i / 2) % 2, base = digits ? 10 : (i ? 24 : 18);
		const int a = digits ? g[i] - '0' : toupper((unsigned char) g[i]) - 'A';
		const int b = digits ? g[i + 1] - '0' : toupper((unsigned char) g[i + 1]) - 'A';
		if ((a < 0) || (a >= base) || (b < 0) || (b >= base)) {
			return -1;
		}
		if (i) {
			dlon /= base;
			dlat /= base;
		}
		lon += a * dlon;
		lat += b * dlat;
	}
	const double phi = (lat + dlat / 2) * M_PI / 180, lam = (lon + dlon / 2) * M_PI / 180;
	const double f = 1 / 298.257223563, e2 = f * (2 - f);
	const double r = 6378137.0 / sqrt(1 - e2 * sin(phi) * sin(phi));
	x[0] = (r + h) * cos(phi) * cos(lam);
	x[1] = (r + h) * cos(phi) * sin(lam);
	x[2] = (r * (1 - e2) + h) * sin(phi);
	if (arg[n]) {
		*height = h;
	}
	printf("held position: a %zu-character grid square's center, surveyed from there\n", n);
	return 0;
}

/* The height above the WGS 84 ellipsoid of an ECEF position (Bowring's closed form). */
static double ellipsoid_height(const double x[3])
{
	const double a = 6378137.0, f = 1 / 298.257223563, b = a * (1 - f);
	const double e2 = f * (2 - f), ep2 = e2 / (1 - e2);
	const double p = hypot(x[0], x[1]), th = atan2(x[2] * a, p * b);
	const double phi = atan2(x[2] + ep2 * b * pow(sin(th), 3), p - e2 * a * pow(cos(th), 3));
	const double n = a / sqrt(1 - e2 * sin(phi) * sin(phi));
	return (fabs(cos(phi)) > 1e-6) ? p / cos(phi) - n : fabs(x[2]) - b;
}

/*
 * The warm start: what the app learned (its almanac, the held position once surveyed or given,
 * and the carrier offset a tracked satellite shows: the TCXO, the tuner's offset and the quarter
 * shift together), kept in a small text file between runs and handed back at the start with the
 * GPS time from this computer's clock, so the first search predicts the visible satellites at
 * three Dopplers each instead of searching every PRN at every Doppler (gps_app.c, message 'W').
 * The file holds the antenna's position: it is written readable by its owner only.
 */
#define WARM_ALM_MAX 64
/* An almanac serves acquisition for weeks; its entries are moved to the current week (alm_to_week). */
#define WARM_MAX_AGE (28 * 86400)
#define GPS_UTC_LEAP 18 /* GPS time less UTC since 2017 */

struct warm {
	time_t saved;
	int have_pos, have_bias, alm_len;
	double pos[3], bias;
	int have_alm[32];
	uint8_t alm[32][WARM_ALM_MAX];
};

static int warm_load(const char* path, struct warm* w)
{
	memset(w, 0, sizeof(*w));
	FILE* f = fopen(path, "r");
	if (!f) {
		return -1;
	}
	char line[512], hex[2 * WARM_ALM_MAX + 2];
	int i;
	while (fgets(line, sizeof(line), f)) {
		long long t;
		if (sscanf(line, "saved %lld", &t) == 1) {
			w->saved = (time_t) t;
		} else if (sscanf(line, "position %lf %lf %lf", &w->pos[0], &w->pos[1], &w->pos[2]) == 3) {
			w->have_pos = 1;
		} else if (sscanf(line, "bias %lf", &w->bias) == 1) {
			w->have_bias = 1;
		} else if ((sscanf(line, "alm %d %129s", &i, hex) == 2) && (i >= 0) && (i < 32)) {
			const int n = (int) strlen(hex) / 2;
			if ((n < 1) || (n > WARM_ALM_MAX) || (w->alm_len && (n != w->alm_len))) {
				continue;
			}
			int bad = 0;
			for (int k = 0; k < n; k++) {
				unsigned v = 0;
				bad |= !isxdigit((unsigned char) hex[2 * k]) ||
					!isxdigit((unsigned char) hex[2 * k + 1]) ||
					(sscanf(&hex[2 * k], "%2x", &v) != 1);
				w->alm[i][k] = (uint8_t) v;
			}
			if (bad) {
				w->have_alm[i] = 0;
				continue;
			}
			w->alm_len = n;
			w->have_alm[i] = 1;
		}
	}
	fclose(f);
	return 0;
}

/* GPS time of week now, from this computer's clock (kept by NTP: within milliseconds). */
static double gps_time_of(double unix_s)
{
	return unix_s - 315964800.0 + GPS_UTC_LEAP;
}

static double gps_tow_now(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	return fmod(gps_time_of((double) ts.tv_sec + ts.tv_nsec * 1e-9), GPS_WEEK);
}

static double wrap_pi(double a)
{
	return a - 2 * M_PI * floor((a + M_PI) / (2 * M_PI));
}

/*
 * An almanac entry (struct gps_alm's bytes) moved whole weeks on, to the week nearest now. Its
 * reference time is in seconds of a week, and the app takes it within half a week of now; an entry
 * saved earlier than that keeps its week seconds and moves its mean anomaly by the orbit's mean
 * motion and its node by the node's drift less the Earth's turn, over the weeks between (a week of
 * the Earth's turn is not whole turns). The same orbit, so the app's arithmetic is unchanged. The
 * reference time's week is the one that puts it nearest the time the entry was saved. Returns the
 * weeks moved.
 */
static int alm_to_week(uint8_t* entry, time_t saved)
{
	float f[8]; /* toa, e, i0, omegadot, sqrta, omega0, omega, m0 */
	memcpy(f, entry, sizeof(f));
	const double g_saved = gps_time_of((double) saved);
	double toa = g_saved - fmod(g_saved, GPS_WEEK) + f[0];
	toa += (toa - g_saved > GPS_WEEK / 2) ? -GPS_WEEK : (g_saved - toa > GPS_WEEK / 2) ? GPS_WEEK : 0;
	const int weeks = (int) floor((gps_time_of((double) time(NULL)) - toa) / GPS_WEEK + 0.5);
	if (weeks == 0) {
		return 0;
	}
	const double dt = weeks * GPS_WEEK, a = (double) f[4] * f[4];
	f[7] = (float) wrap_pi(f[7] + sqrt(3.986005e14 / (a * a * a)) * dt);
	f[5] = (float) wrap_pi(f[5] + (f[3] - 7.2921151467e-5) * dt);
	memcpy(entry, f, sizeof(f));
	return weeks;
}

/* The almanac, then the time and offset: the app's first search is then a predicted one. Returns
 * the PRNs sent when the offset went too, else 0. */
static int warm_send(const struct warm* w)
{
	int n = 0, moved = 0;
	for (int i = 0; i < 32; i++) {
		uint8_t msg[2 + WARM_ALM_MAX] = {'W', (uint8_t) i};
		uint8_t r[128];
		if (w->have_alm[i]) {
			memcpy(&msg[2], w->alm[i], w->alm_len);
			moved |= alm_to_week(&msg[2], w->saved);
			n += app_message(msg, 2 + w->alm_len, r, sizeof(r)) == w->alm_len;
		}
	}
	if (n && w->have_bias) {
		uint8_t msg[17] = {'W'};
		const double t = gps_tow_now();
		memcpy(&msg[1], &t, 8);
		memcpy(&msg[9], &w->bias, 8);
		app_simple(msg, 17, "warm start");
	}
	printf("warm start: an almanac of %d PRNs%s, saved %.0f minutes ago%s%s\n",
	       n,
	       moved ? " (moved to this week)" : "",
	       difftime(time(NULL), w->saved) / 60,
	       w->have_pos ? ", the held position" : "",
	       w->have_bias ? ", the carrier offset" : "");
	return w->have_bias ? n : 0;
}

/* What the app learned, over what the file held, written readable by its owner only. */
static void warm_save(const char* path, const struct warm* old)
{
	struct warm w = *old;
	uint8_t r[128];
	const uint8_t q = 'W';
	if (app_message(&q, 1, r, sizeof(r)) != 34) {
		return;
	}
	if (r[24] >= 2) { /* surveyed or given */
		memcpy(w.pos, r, 24);
		w.have_pos = 1;
	}
	if (r[25]) {
		memcpy(&w.bias, &r[26], 8);
		w.have_bias = 1;
	}
	int n = 0;
	for (int i = 0; i < 32; i++) {
		const uint8_t m[2] = {'W', (uint8_t) i};
		const int len = app_message(m, 2, r, sizeof(r));
		/* struct gps_alm: eight floats, health, then have */
		if ((len >= 34) && (len <= WARM_ALM_MAX) && r[33]) {
			memcpy(w.alm[i], r, len);
			w.alm_len = len;
			w.have_alm[i] = 1;
		}
		n += w.have_alm[i];
	}
	char tmp[1200];
	snprintf(tmp, sizeof(tmp), "%s.tmp", path);
	char* slash = strrchr(tmp, '/');
	if (slash && (slash != tmp)) { /* its directory, made if missing (gps/data on a fresh clone) */
		*slash = 0;
		mkdir(tmp, 0700);
		*slash = '/';
	}
	const int fd = open(tmp, O_WRONLY | O_CREAT | O_TRUNC, 0600);
	FILE* f = (fd >= 0) ? fdopen(fd, "w") : NULL;
	if (!f) {
		fprintf(stderr, "%s: cannot write the warm-start state\n", tmp);
		return;
	}
	fprintf(f, "hackrf_gps state 1\nsaved %lld\n", (long long) time(NULL));
	if (w.have_pos) {
		fprintf(f, "position %.3f %.3f %.3f\n", w.pos[0], w.pos[1], w.pos[2]);
	}
	if (w.have_bias) {
		fprintf(f, "bias %.3f\n", w.bias);
	}
	for (int i = 0; i < 32; i++) {
		if (w.have_alm[i]) {
			fprintf(f, "alm %d ", i);
			for (int k = 0; k < w.alm_len; k++) {
				fprintf(f, "%02x", w.alm[i][k]);
			}
			fprintf(f, "\n");
		}
	}
	if ((fclose(f) == 0) && (rename(tmp, path) == 0)) {
		printf("warm-start state saved: an almanac of %d PRNs%s%s\n",
		       n,
		       w.have_pos ? ", the held position" : "",
		       w.have_bias ? ", the carrier offset" : "");
	} else {
		fprintf(stderr, "%s: cannot write the warm-start state\n", path);
		unlink(tmp);
	}
}

/*
 * The board to run on: the one HackRF Pro whose firmware takes apps (its APP_STATUS answers), opened
 * twice: usb for the app's requests, *dev for the stream and the radio. A HackRF One, or a Pro on
 * stock firmware, shares its USB IDs and is passed over: it may run beside it, its CLKIN from the
 * Pro's CLKOUT. Each is asked over EP0, which a board streaming for another program still answers.
 * Two such Pros are refused. 0, or -1 with the reason printed.
 */
static int open_pro(hackrf_device** dev)
{
	const uint8_t in = LIBUSB_ENDPOINT_IN | LIBUSB_REQUEST_TYPE_VENDOR | LIBUSB_RECIPIENT_DEVICE;
	libusb_device** all;
	const ssize_t n = (libusb_init(NULL) == 0) ? libusb_get_device_list(NULL, &all) : -1;
	int pros = 0;
	uint8_t bus = 0, addr = 0;
	for (ssize_t i = 0; i < n; i++) {
		struct libusb_device_descriptor d;
		libusb_device_handle* h;
		if ((libusb_get_device_descriptor(all[i], &d) != 0) || (d.idVendor != 0x1d50) ||
		    (d.idProduct != 0x6089) || (libusb_open(all[i], &h) != 0)) {
			continue;
		}
		uint8_t id = BOARD_ID_UNDETECTED, st[44];
		const int app = (libusb_control_transfer(h, in, REQ_BOARD_ID, 0, 0, &id, 1, 500) == 1) &&
			(id == BOARD_ID_PRALINE) &&
			(libusb_control_transfer(h, in, REQ_APP_STATUS, 0, 0, st, sizeof(st), 500) ==
			 (int) sizeof(st));
		if (app && (pros++ == 0)) {
			usb = h;
			bus = libusb_get_bus_number(all[i]);
			addr = libusb_get_device_address(all[i]);
		} else {
			libusb_close(h);
		}
	}
	if (n >= 0) {
		libusb_free_device_list(all, 1);
	}
	if (pros != 1) {
		fprintf(stderr,
			pros ? "%d HackRF Pros with app support attached: leave one\n"
			     : "no HackRF Pro with app support (the firmware's APP_STATUS request)\n",
			pros);
		return -1;
	}
	/* The same board through libhackrf, by its bus and address. */
	hackrf_device_list_t* list = hackrf_device_list();
	int rc = -1;
	for (int i = 0; list && (i < list->devicecount) && rc; i++) {
		libusb_device* const u = list->usb_devices[list->usb_device_index[i]];
		if ((libusb_get_bus_number(u) == bus) && (libusb_get_device_address(u) == addr)) {
			rc = (hackrf_device_list_open(list, i, dev) == HACKRF_SUCCESS) ? 0 : -1;
		}
	}
	if (list) {
		hackrf_device_list_free(list);
	}
	if (rc) {
		fprintf(stderr, "the HackRF Pro would not open (in use by another program?)\n");
	}
	return rc;
}

/* An option's whole-number argument within lo..hi into *v; -1 with the reason printed if not. */
static int int_arg(int opt, const char* arg, long lo, long hi, int* v)
{
	char* end;
	const long x = strtol(arg, &end, 10);
	if ((end == arg) || *end || (x < lo) || (x > hi)) {
		fprintf(stderr, "-%c: a whole number from %ld to %ld\n", opt, lo, hi);
		return -1;
	}
	*v = (int) x;
	return 0;
}

static void usage(FILE* f)
{
	fputs("Usage: hackrf_gps [options]\n"
	      "       hackrf_gps -f recording.iq [-k snapshots] [-s offset_ms]\n"
	      "\n"
	      "Tracks GPS L1 C/A on a HackRF Pro with the GPS gateware and M4 app, measures\n"
	      "the reference oscillator's frequency error against the satellites' carriers,\n"
	      "and can steer PLL B (CLKOUT and the RF references) to cancel it. It runs on the\n"
	      "one HackRF Pro with app support attached; other HackRFs are left alone.\n"
	      "\n"
	      "Run:\n"
	      "  -t SECONDS        seconds of tracking before stopping; 0 runs until Ctrl-C\n"
	      "                    (default 60)\n"
	      "  -o FILE           record the stream to FILE, the fixes to FILE.fix and the\n"
	      "                    frequency estimates to FILE.freq\n"
	      "  -h, --help        show this help and exit\n"
	      "\n"
	      "Receiver:\n"
	      "  -l DB             LNA gain, 0-40 dB in 8 dB steps (default 32)\n"
	      "  -g DB             VGA gain, 0-62 dB in 2 dB steps (default 30)\n"
	      "  -a 0|1            RF amplifier off or on (default 1)\n"
	      "  -p 0|1            antenna port power (bias tee) off or on (default 0)\n"
	      "\n"
	      "Position (the frequency estimate needs the antenna's position):\n"
	      "  -P X,Y,Z          hold this ECEF position, in meters\n"
	      "  -P RUN.fix        hold the mean of the fixes in a previous run's .fix file\n"
	      "  -P GRID[,HEIGHT]  start from a Maidenhead grid square's center (4 to 10\n"
	      "                    characters, e.g. FN20xr), optionally HEIGHT meters above\n"
	      "                    the WGS 84 ellipsoid; the survey of 60 fixes then\n"
	      "                    replaces it\n"
	      "  -E HEIGHT         antenna height above the WGS 84 ellipsoid (within about\n"
	      "                    100 m of sea level), in meters, held by a three-satellite\n"
	      "                    fix until the first full fix (default 0, or the grid\n"
	      "                    square's or warm start's height)\n"
	      "  -W FILE|none      warm-start state: read at the start (within 4 weeks) and\n"
	      "                    written at the end; none disables it\n"
	      "                    (default gps/data/hackrf_gps.state)\n"
	      "\n"
	      "Clock steering:\n"
	      "  -D SECONDS        turn the steering loop on at this second of tracking\n"
	      "  -B T:A,B,C[;...]  set PLL B's feedback ratio to A + B/C at second T of\n"
	      "                    tracking (nominal 32), for testing steps\n"
	      "\n"
	      "Loading:\n"
	      "  -G FILE           GPS gateware pack (default fpga-build/packs/gps_rx.pack)\n"
	      "  -A FILE           GPS M4 app (default apps/build/gps_app.bin)\n"
	      "  -N                keep the gateware and app already loaded\n"
	      "  -K                keep the GPS image loaded at the end instead of returning\n"
	      "                    to stock image 0\n"
	      "\n"
	      "Search:\n"
	      "  -S SECONDS        print the first search's progress until a pass arms a\n"
	      "                    satellite, for at most this long (default 120)\n"
	      "  -H                search on the host instead of on the box (a reference)\n"
	      "  -k N              snapshots per search with -H or -f (default 6)\n"
	      "  -n N              with -H, fill up to N channels with the next-best PRNs\n"
	      "  -c CHANNEL        with -H, the first channel to arm (default 0)\n"
	      "\n"
	      "Recording (-f):\n"
	      "  -f FILE           search a recording instead of a device: 8-bit I/Q at\n"
	      "                    4 MS/s, as hackrf_transfer writes it\n"
	      "  -s MS             offset of the first snapshot into the recording, in\n"
	      "                    milliseconds (default 100)\n"
	      "\n"
	      "Defaults for -G, -A and -W are relative to this program's repository.\n",
	      f);
}

int main(int argc, char** argv)
{
	const char* out_path = NULL;
	const char* file_path = NULL;
	int offset_ms = 100, nsnaps = 6, fill = 0, first_ch = 0;
	int seconds = 60, lna = 32, vga = 30, amp = 1, bias = 0, opt;
	double held[3];
	int have_held = 0; /* 1 given, 2 a grid square's center to survey from */
	/* -E: the height a fix from three satellites holds, until the first kept fix. */
	double height = 0, grid_h = 0;
	int grid_seed = 0; /* the held position a grid square's center: three satellites may replace it */
	int have_height = 0; /* 1 a grid square's, 2 -E, 3 the warm start's */
	int box_search = 120; /* -S: the first search printed until an arm, at most this long */
	int host_search = 0;  /* -H: the search here instead */
	/* -G and -A: the GPS image and app loaded at the start (the repository's, found from this
	 * program's place in it, by default); -N: whatever is loaded kept; -K: GPS kept at the end. */
	char pack[1100], app[1100], exe[1024];
	const ssize_t el = readlink("/proc/self/exe", exe, sizeof(exe) - 1);
	exe[el > 0 ? el : 0] = 0;
	char* cut = strstr(exe, "/gps/tools/build/");
	if (cut) {
		*cut = 0;
	}
	snprintf(pack, sizeof(pack), "%s/fpga-build/packs/gps_rx.pack", exe);
	snprintf(app, sizeof(app), "%s/apps/build/gps_app.bin", exe);
	/* -W: the warm-start state (gps/data is not kept in the repository); -W none: none. */
	char warm_path[1100];
	snprintf(warm_path, sizeof(warm_path), "%s/gps/data/hackrf_gps.state", exe);
	int load = 1, keep = 0;
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--help")) {
			usage(stdout);
			return 0;
		}
	}
	while ((opt = getopt(argc, argv, "o:t:l:g:a:p:f:s:k:n:c:P:E:B:D:S:HG:A:NKW:h")) != -1) {
		switch (opt) {
		case 'E': {
			char* end;
			height = strtod(optarg, &end);
			if ((end == optarg) || *end || !isfinite(height)) {
				fprintf(stderr, "-E: a height in meters\n");
				return 2;
			}
			have_height = 2;
			break;
		}
		case 'W':
			snprintf(warm_path, sizeof(warm_path), "%s", strcmp(optarg, "none") ? optarg : "");
			break;
		case 'G':
			snprintf(pack, sizeof(pack), "%s", optarg);
			break;
		case 'A':
			snprintf(app, sizeof(app), "%s", optarg);
			break;
		case 'N':
			load = 0;
			break;
		case 'K':
			keep = 1;
			break;
		case 'S':
			if (int_arg(opt, optarg, 0, 86400, &box_search)) {
				return 2;
			}
			break;
		case 'H':
			host_search = 1;
			break;
		case 'D':
			if (int_arg(opt, optarg, 0, INT_MAX, &loop_at)) {
				return 2;
			}
			break;
		case 'B':
			if (parse_steps(optarg) != 0) {
				fprintf(stderr,
					"-B: t:a,b,c;t:a,b,c... (t seconds of tracking from 1, ratio a + b/c)\n");
				return 2;
			}
			break;
		case 'P':
			if (sscanf(optarg, "%lf,%lf,%lf", &held[0], &held[1], &held[2]) == 3) {
				have_held = 1;
			} else if (grid_position(optarg, held, &grid_h) == 0) {
				have_held = 2;
				grid_seed = 1;
				if (strchr(optarg, ',') && (have_height < 2)) {
					height = grid_h; /* -E, given, wins */
					have_height = 1;
				}
			} else if (mean_fix(optarg, held) == 0) {
				have_held = 1;
			} else {
				fprintf(stderr,
					"-P: x,y,z in meters (ECEF), a grid square (FN20xr[,height]), or a .fix "
					"file with fixes\n");
				return 2;
			}
			break;
		case 'f':
			file_path = optarg;
			break;
		case 's':
			if (int_arg(opt, optarg, 0, INT_MAX / N, &offset_ms)) {
				return 2;
			}
			break;
		case 'c':
			if (int_arg(opt, optarg, 0, CHANNELS - 1, &first_ch)) {
				return 2;
			}
			break;
		case 'n':
			if (int_arg(opt, optarg, 0, CHANNELS, &fill)) {
				return 2;
			}
			break;
		case 'k':
			if (int_arg(opt, optarg, 1, MAX_SNAPS, &nsnaps)) {
				return 2;
			}
			break;
		case 'o':
			out_path = optarg;
			break;
		case 't':
			if (int_arg(opt, optarg, 0, INT_MAX, &seconds)) {
				return 2;
			}
			break;
		case 'l':
			if (int_arg(opt, optarg, 0, 40, &lna)) {
				return 2;
			}
			break;
		case 'g':
			if (int_arg(opt, optarg, 0, 62, &vga)) {
				return 2;
			}
			break;
		case 'a':
			if (int_arg(opt, optarg, 0, 1, &amp)) {
				return 2;
			}
			break;
		case 'p':
			if (int_arg(opt, optarg, 0, 1, &bias)) {
				return 2;
			}
			break;
		case 'h':
			usage(stdout);
			return 0;
		default:
			fprintf(stderr, "Try 'hackrf_gps -h' for more information.\n");
			return 2;
		}
	}
	struct snap s[MAX_SNAPS];
	if (file_path) {
		FILE* in = fopen(file_path, "rb");
		int8_t* raw = malloc(2 * SNAP_LEN);
		if (!raw) {
			return 1;
		}
		for (int j = 0; j < nsnaps; j++) {
			const long at = (long) (offset_ms + j * SNAP_GAP_MS) * N;
			if (!in || fseek(in, at * 2, SEEK_SET) ||
			    fread(raw, 2, SNAP_LEN, in) != SNAP_LEN) {
				fprintf(stderr,
					"%s: cannot read %d samples at %ld\n",
					file_path,
					SNAP_LEN,
					at);
				return 1;
			}
			float complex* x = malloc(SNAP_LEN * sizeof(float complex));
			if (!x) {
				return 1;
			}
			for (int i = 0; i < SNAP_LEN; i++) {
				x[i] = raw[2 * i] + I * raw[2 * i + 1];
			}
			s[j] = (struct snap){x, (uint32_t) at};
		}
		struct found f[NPRN];
		const int n = search_and_report(s, nsnaps, f);
		/* The arms as a live run would send them, 2 s ahead, in recording sample indices. */
		for (int p = 0; p < n; p++) {
			uint8_t msg[24];
			double metric;
			const uint32_t start1 =
				plan_arm(s, nsnaps, &f[p], 2.0, p, msg, &metric);
			printf("ARM %d %u %u %u %d\n",
			       f[p].prn,
			       start1,
			       get32(&msg[3]),
			       get32(&msg[7]),
			       (int32_t) get32(&msg[20]));
		}
		return 0;
	}
	/* Without -o the stream and the fix and frequency files are not kept. */
	out_file = fopen(out_path ? out_path : "/dev/null", "wb");
	if (!out_file) {
		perror(out_path);
		return 1;
	}

	hackrf_device* dev = NULL;
	if ((hackrf_init() != HACKRF_SUCCESS) || open_pro(&dev)) {
		hackrf_exit();
		return 1;
	}
	/* No SA_RESTART: a signal ends the sleeps early, and every wait looks at the flag. */
	struct sigaction sa = {0};
	sa.sa_handler = on_signal;
	sigemptyset(&sa.sa_mask);
	sigaction(SIGINT, &sa, NULL);
	sigaction(SIGTERM, &sa, NULL);
	sigaction(SIGHUP, &sa, NULL);
	/* Output piped to a reader that has gone (tee after Ctrl-C) must not end the run before its
	 * cleanup. */
	signal(SIGPIPE, SIG_IGN);
	int status = 0, gps_running = 0, warm_have = 0;
	struct warm warm;
	if (load && switch_to_gps(pack, app)) {
		status = 1;
		goto done;
	}
	int rc = hackrf_set_freq(dev, 1575420000ULL) | hackrf_set_sample_rate(dev, FS) |
		hackrf_set_baseband_filter_bandwidth(
			 dev,
			 hackrf_compute_baseband_filter_bw(3000000)) |
		hackrf_set_lna_gain(dev, lna) | hackrf_set_vga_gain(dev, vga) |
		hackrf_set_amp_enable(dev, amp) | hackrf_set_antenna_enable(dev, bias);
	if (rc != HACKRF_SUCCESS) {
		fprintf(stderr, "radio settings failed\n");
		status = 1;
		goto done;
	}

	capture_max = 1 << 22; /* words; 16 snapshots need about 2.1 million */
	capture = malloc(capture_max * sizeof(uint16_t));
	if (!capture) {
		status = 1;
		goto done;
	}
	/* Records mode before receiving starts, so the stream holds only what is asked for: with
	 * no channel armed, the snapshots. */
	const uint8_t mode[2] = {'M', 1};
	stop_all();
	if (app_simple(mode, 2, "mode")) {
		status = 1;
		goto done;
	}
	usleep(100000);
	if (hackrf_start_rx(dev, rx_callback, NULL) != HACKRF_SUCCESS) {
		fprintf(stderr, "cannot start receiving\n");
		status = 1;
		goto done;
	}

	usleep(200000);
	/* The tuner's offset, which the frequency measurement takes out, and the quarter shift, which
	 * separates PLL B's LO from PLL A's share (gps_app.c, 'F'). */
	{
		const double off = tuner_offset();
		const int32_t uhz = (int32_t) llround(off * 1e6), q = quarter_shift();
		uint8_t msg[9] = {'F'};
		put32(&msg[1], (uint32_t) uhz);
		put32(&msg[5], (uint32_t) q);
		printf("tuner offset from 1575.42 MHz: %+.6f Hz; quarter shift %+d Hz\n", off, q);
		if (app_simple(msg, 9, "tuner offset")) {
			status = 1;
			goto done;
		}
	}
	/* The warm start: a position from it seeds the survey when -P gives none. */
	warm_have = warm_path[0] && !host_search && (warm_load(warm_path, &warm) == 0) &&
		(difftime(time(NULL), warm.saved) < WARM_MAX_AGE);
	if (warm_have && warm.have_pos && !have_held) {
		memcpy(held, warm.pos, 24);
		have_held = 2;
		printf("held position: the warm start's, surveyed from there\n");
		if (have_height < 2) {
			height = ellipsoid_height(held);
			have_height = 3; /* the warm start's: not printed (it comes from the position) */
		}
	}
	/* The height a fix from three satellites holds (message 'H'; the app's default is 0), and
	 * whether such a fix may replace a grid square's center (a warm start's surveyed position
	 * stays). */
	if ((have_height || grid_seed) && (have_held != 1)) {
		uint8_t msg[6] = {'H'};
		const float h = (float) height;
		memcpy(&msg[1], &h, 4);
		msg[5] = (uint8_t) grid_seed;
		if (have_height == 3) {
			printf("three satellites before the first fix: the warm start's height held\n");
		} else {
			printf("three satellites before the first fix: the height held at %.0f m%s\n", height,
			       grid_seed ? ", replacing the grid square's center" : "");
		}
		if (app_simple(msg, 6, "held height")) {
			status = 1;
			goto done;
		}
	}
	/* The frequency estimate's held position: given, or surveyed afresh over 60 kept fixes,
	 * from a grid square's center if one is given (held until the first kept fix). */
	{
		uint8_t msg[29] = {'U'};
		int len = 1;
		if (have_held) {
			memcpy(&msg[1], held, 24);
			len = 25;
		}
		if (have_held != 1) {
			put32(&msg[len], 60);
			len += 4;
		}
		if (app_simple(msg, len, "held position")) {
			status = 1;
			goto done;
		}
		gps_running = 1;
	}
	/* The first search's blocks per PRN and Doppler. A cold pass over 8 finds the strong
	 * satellites, and a pass that arms nothing goes on at twice the blocks. A warm start with a
	 * whole almanac, a position and the carrier offset makes every pass a predicted one (about 8
	 * PRNs at 3 Dopplers). Such a pass is short enough at 128 blocks to find the weak satellites a
	 * cold pass misses. */
	int first_blocks = 8;
	if (warm_have && (warm_send(&warm) == 32) && have_held) { /* 32 only with the offset */
		first_blocks = 128;
	}
	double base_rate = 0;
	int armed = 0;
	if (!host_search) {
		/* The search on the box (gps_app.c, 'C'), progress once a second ('K'). */
		measure_base_rate(&base_rate);
		/* Each candidate is checked over at least 64 blocks. A satellite near the cold pass's
		 * limit passes its check more often over 64 blocks than over 32. */
		uint8_t k[128] = {0};
		const uint8_t kq = 'K';
		if (!load && (app_message(&kq, 1, k, sizeof(k)) >= 21) && k[0]) {
			/* An app kept from an earlier run (-K) searches on: a 'C' now would be refused. */
			printf("search on the box: the kept app's, already under way\n");
		} else {
			uint8_t msg[3] = {'C', (uint8_t) first_blocks, 64};
			if (app_simple(msg, 3, "search on the box")) {
				status = 1;
				goto done;
			}
			printf("search on the box: %d blocks a PRN and Doppler, candidates checked over %d\n",
			       first_blocks,
			       first_blocks > 64 ? first_blocks : 64);
		}
		int searching = 1;
		for (int sec = 1; (sec <= box_search) && searching && !interrupted; sec++) {
			sleep(1);
			const uint8_t q = 'K';
			const int n = app_message(&q, 1, k, sizeof(k));
			if (n < 21) {
				continue;
			}
			searching = k[0];
			printf("%3d s: %s (%d blocks)", sec, search_mode(k), k[19]);
			if (k[0] != 2) {
				printf(", PRN %d, Doppler bin %d", k[1] + 1, k[2]);
			}
			printf("; candidates");
			for (int i = 0; i < k[3]; i++) {
				printf(" %d", k[4 + i] + 1);
			}
			printf("; armed");
			for (int i = 0; i < k[11] && i < CHANNELS; i++) {
				printf(" %d", k[12 + i]);
			}
			if (!searching && (n >= 29)) {
				const uint32_t work = get32(&k[21]), wait = get32(&k[25]);
				printf("; %u ticks working, %u waiting for a snapshot (%.0f%%)", work, wait,
				       100.0 * wait / (work + wait ? work + wait : 1));
			}
			printf("\n");
			fflush(stdout);
			/* Tracking from the first pass's arms: the box goes on searching while fewer
			 * than four channels track. */
			if (k[11] && (k[0] != 2)) {
				break;
			}
		}
		goto track;
	}
	double t_snap0;
	for (int j = 0; j < MAX_SNAPS; j++) {
		s[j].x = NULL;
	}
	if (interrupted || take_snapshots(nsnaps, s, &t_snap0)) {
		status = 1;
		goto done; /* stop receiving: the transceiver would stay on */
	}
	struct found f[NPRN];
	search_and_report(s, nsnaps, f);
	if (interrupted) {
		goto done;
	}
	int count = 0;
	while ((count < NPRN) && (count < CHANNELS) && detected(f, count, nsnaps)) {
		count++;
	}
	const int found = count;
	/* -n: the next-best PRNs too, up to n channels; they track noise, which keeps their loops
	 * writing every millisecond as in pull-in (a load test). */
	if (fill > count) {
		count = fill < CHANNELS ? fill : CHANNELS;
	}
	measure_base_rate(&base_rate);

	/*
	 * The search takes seconds per snapshot, and the code drifts by the Doppler error over that
	 * time: about 100 Hz, 0.065 chip a second. Fresh snapshots, searched only near each
	 * detection, place the arms about a second after them instead; as many as the search had,
	 * so a weak satellite stands out as well.
	 */
	const int fresh = nsnaps;
	if (count && (interrupted || take_snapshots(fresh, s, &t_snap0))) {
		status = 1;
		goto done;
	}
	for (int p = 0; (p < count) && !interrupted; p++) {
		uint8_t msg[24];
		/* snapshot 0 began just after its request */
		const double ahead = (now_s() - t_snap0) + 1.0;
		double metric;
		const uint32_t start1 = plan_arm(
			s,
			fresh,
			&f[p],
			ahead,
			(first_ch + armed) % CHANNELS,
			msg,
			&metric);
		/* A detection gone from the fresh snapshots would be armed on a noise peak. */
		if ((p < found) && (metric < 1 + 10.0 / sqrt(BLOCKS * fresh))) {
			printf("  PRN %d not armed: under the threshold in the fresh snapshots\n",
			       f[p].prn);
			continue;
		}
		if (app_simple(msg, 24, "arm") == 0) {
			printf("  channel %d: PRN %d, Doppler %+.3f Hz, first sample %u\n",
			       (first_ch + armed) % CHANNELS,
			       f[p].prn,
			       (int32_t) get32(&msg[20]) * 1e-3,
			       start1);
			armed++;
		}
	}

track:;
	/* Track, with the app's status and fix once a second. */
	const int dropped0 = fpga_read(REG_DROPPED);
	char fix_path[1024];
	snprintf(fix_path, sizeof(fix_path), out_path ? "%s.fix" : "/dev/null", out_path);
	FILE* fix_file = fopen(fix_path, "w");
	if (!fix_file) {
		perror(fix_path);
		status = 1;
		goto done;
	}
	/* An app left loaded (-N) keeps its last estimate, fix and PLL B writes from an earlier run:
	 * print only newer ones. */
	{
		uint8_t r[128];
		uint8_t q = 'R';
		if (app_message(&q, 1, r, sizeof(r)) == 64) {
			last_estimates = get32(&r[0]);
		}
		q = 'P';
		if (app_message(&q, 1, r, sizeof(r)) == 128) {
			last_fix.count = get32(&r[0]);
		}
		q = 'J';
		if (app_message(&q, 1, r, sizeof(r)) == 116) {
			last_writes = get32(&r[0]);
		}
	}
	char freq_path[1024];
	snprintf(freq_path, sizeof(freq_path), out_path ? "%s.freq" : "/dev/null", out_path);
	FILE* freq_file = fopen(freq_path, "w");
	if (!freq_file) {
		perror(freq_path);
		fclose(fix_file);
		status = 1;
		goto done;
	}
	for (int sec = 0; ((seconds == 0) || (sec < seconds)) && !interrupted; sec++) {
		sleep(1);
		if (interrupted) {
			break;
		}
		for (int i = 0; i < nsteps; i++) {
			if (steps[i].sec == sec + 1) {
				ask_ratio(steps[i].a, steps[i].b, steps[i].c);
			}
		}
		if (loop_at && (sec + 1 == loop_at)) {
			const uint8_t on[2] = {'D', 1};
			if (app_simple(on, 2, "steering loop on") == 0) {
				loop_on = 1;
				printf("  steering loop on\n");
			}
		}
		if (nsteps || loop_on) {
			report_writes();
		}
		uint8_t r[128];
		const uint8_t q = 'Q';
		const int n = app_message(&q, 1, r, sizeof(r));
		if (n != 128) {
			printf("status: no reply\n");
			continue;
		}
		uint8_t st[16 * CHANNELS];
		memcpy(st, r, sizeof(st));
		printf("%3d s:", sec + 1);
		/* Every channel tracking, whoever armed it: the search, the manager or this tool. */
		for (int c = 0; c < CHANNELS; c++) {
			const uint8_t* e = &st[16 * c];
			if (!e[0]) {
				continue;
			}
			printf(" [%d PRN %u: %u rec, edge %d, %+.1f Hz, P %u, gaps %u, lost %u]",
			       c,
			       e[0],
			       get32(&e[4]),
			       e[2] == 255 ? -1 : e[2],
			       (int32_t) get32(&e[8]) * 1e-3,
			       get32(&e[12]),
			       e[1],
			       e[3]);
		}
		/* A search the manager started (the almanac's predictions first, then the PRNs it
		 * placed above the horizon at three Dopplers each, those without an almanac at all). */
		{
			uint8_t k[128];
			const uint8_t kq = 'K';
			if ((app_message(&kq, 1, k, sizeof(k)) >= 21) && k[0]) {
				printf(" [search: %s (%d blocks)", search_mode(k), k[19]);
				if (k[0] != 2) {
					printf(", PRN %d, Doppler bin %d]", k[1] + 1, k[2]);
				} else {
					/* Refining: the PRN index is past the pass's last. */
					printf(", candidates");
					for (int i = 0; i < k[3] && i < CHANNELS; i++) {
						printf(" %d", k[4 + i] + 1);
					}
					printf("]");
				}
			}
		}
		printf(" overruns %u, garbage %u",
		       r[120] | (r[121] << 8),
		       get32(&r[116]));
		const int d = fpga_read(REG_DROPPED);
		if ((d >= 0) && (dropped0 >= 0)) {
			printf(", dropped %d", (d - dropped0) & 0xFF);
		}
		const uint8_t l = 'L';
		const int ln = app_message(&l, 1, r, sizeof(r));
		if ((base_rate > 0) && (ln >= 12) && get32(&r[4])) {
			const double rate = (double) get32(&r[0]) / get32(&r[4]);
			printf(", M4 load %.0f%%, longest pass %.0f us",
			       100 * (1 - rate / base_rate),
			       get32(&r[8]) / 204.0);
		}
		if (ln >= 16) {
			printf(", stack %u of 4096 bytes", get32(&r[12]));
		}
		printf("\n");
		report_channels(st);
		report_fix(fix_file);
		report_freq(freq_file);
		fflush(stdout);
	}
	if (fclose(fix_file) | fclose(freq_file)) {
		fprintf(stderr, "%s, %s: a write failed\n", fix_path, freq_path);
		status = 1;
	}

done:
	if (gps_running && warm_path[0]) {
		warm_save(warm_path, warm_have ? &warm : &(struct warm){0});
	}
	if (interrupted) {
		printf("interrupted: stopping\n");
		status = 130; /* as the shell reports a SIGINT */
	}
	restore_ratio();
	stop_all();
	/* Off before receiving stops, and again by the stop itself (the idle settings). */
	if (hackrf_set_antenna_enable(dev, 0) != HACKRF_SUCCESS) {
		fprintf(stderr, "antenna power off failed\n");
		status = 1;
	}
	hackrf_stop_rx(dev);
	if (!keep) {
		const double t0 = now_s();
		app_stop_wait();
		const double t1 = now_s();
		if (hackrf_set_fpga_bitstream(dev, 0) == HACKRF_SUCCESS) {
			const double t2 = now_s();
			printf("switched to stock image 0 in %.0f ms: app stopped %.0f ms, image %.0f ms\n",
			       (t2 - t0) * 1e3,
			       (t1 - t0) * 1e3,
			       (t2 - t1) * 1e3);
		} else {
			fprintf(stderr, "stock image 0 did not load\n");
			status = 1;
		}
	}
	hackrf_close(dev);
	hackrf_exit();
	libusb_close(usb);
	libusb_exit(NULL);
	const int closed = fclose(out_file);
	pthread_mutex_lock(&lock);
	if (out_path) {
		printf("%.1f MB written to %s\n", stream_bytes / 1e6, out_path);
		if (stream_failed || closed) {
			fprintf(stderr, "%s: a write failed: the recording is short\n", out_path);
			status = 1;
		}
	}
	pthread_mutex_unlock(&lock);
	return status;
}
