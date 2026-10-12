/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Test apps for the runtime-app sandbox, one per value of TEST. Each line of TESTS in the Makefile
 * gives the cause the firmware must report when the app ends.
 */

#include "app.h"

static volatile uint32_t ticks;
static volatile float phase;
#if TEST == 12
static volatile uint32_t mark, stage; /* the cycle counter at a write's next tick; the step */
	#define GAP_CYCLES (110u * 1000u * 204u) /* past the firmware's 100 ms between PLL B writes */
#endif

#if TEST == 5
static uint32_t recurse(uint32_t n)
{
	volatile uint8_t pad[256];
	pad[0] = (uint8_t) n;
	return (n < 100000 ? recurse(n + 1) : 0) + pad[0]; /* 25 MB deep: never returns */
}
#endif

uint32_t app_main(uint32_t op, uint32_t arg)
{
	(void) arg; /* test 12's only */
	if (op == APP_OP_INIT) {
#if TEST == 8
		return 1; /* init refuses */
#elif TEST == 9
		/* Both read-only regions must read. */
		ticks = *(volatile uint32_t*) 0x20007000 + *(volatile uint32_t*) 0x20008000;
#endif
		return 0;
	}
#if TEST == 0
	if (op == APP_OP_MESSAGE) {
		/* Echo the message reversed, with its length first. */
		APP_REPLY_OUT[0] = (uint8_t) arg;
		for (uint32_t i = 0; i < arg; i++) {
			APP_REPLY_OUT[1 + i] = APP_MESSAGE_IN[arg - 1 - i];
		}
		return arg + 1;
	}
#endif
#if TEST == 12
	if ((op == APP_OP_STOP) && (stage != 5)) {
		__asm__ volatile("udf #11"); /* stopped before every step was taken */
	}
#endif
	if (op != APP_OP_TICK) {
		return 0;
	}
	ticks++;
#if TEST == 0
	/* Runs: floating point, the data region. */
	phase += 0.001f;
	APP_DATA[ticks % 40960] = (uint8_t) ticks;
#elif TEST == 1
	*(volatile uint32_t*) 0x40050000 = 0; /* a peripheral */
#elif TEST == 2
	*(volatile uint32_t*) 0x10080000 = 0; /* the firmware's RAM */
#elif TEST == 3
	for (;;) {
	}
#elif TEST == 4
	__asm__ volatile("svc #5");
#elif TEST == 5
	ticks += recurse(0); /* runs off the bottom of the stack */
#elif TEST == 6
	((void (*)(void))(0x10088000 | 1))(); /* data is not executable */
#elif TEST == 7
	__asm__ volatile("udf #0");
#elif TEST == 9
	*(volatile uint32_t*) 0x20008000 = 0; /* the sample buffer is read only */
#elif TEST == 10
	*(volatile uint32_t*) 0xE000E100 = 1; /* the NVIC */
#elif TEST == 11
	/* The register service: register 6 is queued, register 1 refused. A wrong
	 * answer ends the app with a usage fault, which the check reports. */
	if ((app_fpga_write(6, 0) != 0) || (app_fpga_write(1, 0) != 2)) {
		__asm__ volatile("udf #1");
	}
#elif TEST == 12
	/* The services' arguments, fuzzed: a wrong answer ends the app with a usage fault. Tick 1:
	 * every register but 6 refused, the write queue filled until it says full, PLL B ratios on
	 * and just past each bound (5 ppm of 32), one accepted at the bound and a second refused while
	 * it waits. Tick 2: that write made, and another refused (under 100 ms since it). 110 ms on:
	 * one accepted; 110 ms after that, back to 32. */
	uint32_t st[4];
	if (ticks == 1) {
		for (uint32_t reg = 0; reg < 256; reg++) {
			if ((reg != 6) && (app_fpga_write(reg, 0) != 2)) {
				__asm__ volatile("udf #2");
			}
		}
		if ((app_fpga_write(0x106, 0) != 2) || (app_fpga_write(0xFFFFFFFF, 0) != 2)) {
			__asm__ volatile("udf #3");
		}
		uint32_t n = 0;
		while ((app_fpga_write(6, 0) == 0) && (n < 1000)) {
			n++;
		}
		if ((n == 0) || (n >= 1000)) {
			__asm__ volatile("udf #4");
		}
		static const uint32_t bad[][3] = {
			{32, 0, 0},
			{32, 1, 1},
			{32, 5, 4},
			{32, 0, 1048576},
			{32, 1048574, 1048575},
			{31, 0, 1},
			{33, 0, 1},
			{0, 0, 1},
			{0xFFFFFFFF, 0, 1},
			{32, 161, 1000000},
			{31, 999839, 1000000},
		};
		for (uint32_t i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
			if (app_pll_b_set(bad[i][0], bad[i][1], bad[i][2]) != 2) {
				__asm__ volatile("udf #5");
			}
		}
		app_pll_b_status(st);
		if (st[3] != 0) {
			__asm__ volatile("udf #6");
		}
		if ((app_pll_b_set(31, 999840, 1000000) != 0) || (app_pll_b_set(32, 160, 1000000) != 1)) {
			__asm__ volatile("udf #7");
		}
		app_pll_b_status(st);
		if (st[3] != 1) {
			__asm__ volatile("udf #8");
		}
	} else if (ticks == 2) {
		app_pll_b_status(st);
		if ((st[0] == 0) || (st[3] != 0) || (app_pll_b_set(32, 160, 1000000) != 1)) {
			__asm__ volatile("udf #9");
		}
		mark = arg;
		stage = 1;
	} else if ((stage == 1 || stage == 3) && (arg - mark >= GAP_CYCLES)) {
		if (app_pll_b_set(32, (stage == 1) ? 160 : 0, (stage == 1) ? 1000000 : 1) != 0) {
			__asm__ volatile("udf #10");
		}
		stage++;
	} else if (stage == 2 || stage == 4) {
		mark = arg; /* the tick after the ask: the write is made */
		stage++;
	}
#endif
	return 0;
}
