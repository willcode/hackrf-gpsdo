/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * What a runtime app implements: one function, called by the firmware's
 * main loop (firmware/common/app_sandbox.h has the contract). Init returning
 * nonzero refuses the start.
 */

#ifndef APP_H
#define APP_H

#include <stdint.h>

enum {
	APP_OP_INIT = 0,
	APP_OP_TICK = 1, /* arg: the CPU cycle counter, 204 MHz */
	APP_OP_STOP = 2,
	APP_OP_MESSAGE = 3, /* arg: the message's length at APP_MESSAGE_IN; return the
	                       reply's length at APP_REPLY_OUT (at most 128) */
};

#ifdef APP_HOST
/* Built for a host test: the test supplies the memory and the service. */
extern uint8_t app_host_data[0xA000];
extern uint8_t app_host_samples[0x8000];
extern uint32_t app_host_m0_count;
uint32_t app_fpga_write(uint32_t reg, uint32_t value);
uint32_t app_pll_b_set(uint32_t a, uint32_t b, uint32_t c);
void app_pll_b_status(uint32_t out[4]);
	#define APP_DATA         ((volatile uint8_t*) app_host_data)
	#define APP_DATA_USABLE  (0xA000 - 256)
	#define APP_MESSAGE_IN   ((const volatile uint8_t*) &app_host_data[0xA000 - 256])
	#define APP_REPLY_OUT    ((volatile uint8_t*) &app_host_data[0xA000 - 128])
	#define APP_M0_COUNT     app_host_m0_count
	#define APP_SAMPLES      ((const volatile uint8_t*) app_host_samples)
	#define APP_SAMPLES_SIZE 0x8000
#else

/* The app's 40 kB data region, zeroed at start; its top 256 bytes carry host
 * messages (up to 64 bytes) and the app's replies. */
#define APP_DATA         ((volatile uint8_t*) 0x10088000)
#define APP_DATA_USABLE  (0xA000 - 256)
#define APP_MESSAGE_IN   ((const volatile uint8_t*) 0x10091F00)
#define APP_REPLY_OUT    ((volatile uint8_t*) 0x10091F80)

/* The M0's state and the sample buffer it fills, read only. */
#define APP_M0_STATE     ((const volatile uint32_t*) 0x20007000)
#define APP_M0_COUNT     (APP_M0_STATE[2]) /* bytes written to the sample buffer */
#define APP_SAMPLES      ((const volatile uint8_t*) 0x20008000)
#define APP_SAMPLES_SIZE 0x8000

/* Queue a write of FPGA register reg (6 only); the firmware makes it after
 * this call. 0 queued, 1 queue full (try next tick), 2 not allowed. */
static inline uint32_t app_fpga_write(uint32_t reg, uint32_t value)
{
	register uint32_t r0 __asm__("r0") = reg;
	register uint32_t r1 __asm__("r1") = value;
	__asm__ volatile("svc #1" : "+r"(r0) : "r"(r1) : "memory");
	return r0;
}

/* Ask for PLL B's feedback ratio a + b/c (nominal 32, within 5 ppm of it); the
 * firmware sets it after this call's register writes. 0 accepted, 1 one
 * already waits or the last was under 100 ms ago (try again), 2 out of
 * bounds. When the app ends, PLL B goes back to 32. */
static inline uint32_t app_pll_b_set(uint32_t a, uint32_t b, uint32_t c)
{
	register uint32_t r0 __asm__("r0") = a;
	register uint32_t r1 __asm__("r1") = b;
	register uint32_t r2 __asm__("r2") = c;
	__asm__ volatile("svc #2" : "+r"(r0) : "r"(r1), "r"(r2) : "memory");
	return r0;
}

/* The ratio writes made: out[0] the count, out[1] and out[2] the cycle counter
 * before and after the last, out[3] 1 while one waits. */
static inline void app_pll_b_status(uint32_t out[4])
{
	register uint32_t r0 __asm__("r0");
	register uint32_t r1 __asm__("r1");
	register uint32_t r2 __asm__("r2");
	register uint32_t r3 __asm__("r3");
	__asm__ volatile("svc #3" : "=r"(r0), "=r"(r1), "=r"(r2), "=r"(r3) : : "memory");
	out[0] = r0;
	out[1] = r1;
	out[2] = r2;
	out[3] = r3;
}
#endif

/* State kept in the data region rather than after the image (app.ld): zeroed at
 * start like .bss, and not executable. A .bss section name, so an initializer
 * is a compile error ("only zero initializers are allowed") rather than dropped
 * on the device and kept on the host. */
#ifdef APP_HOST
	#define APP_STATE
#else
	#define APP_STATE __attribute__((section(".bss.appdata")))
#endif

uint32_t app_main(uint32_t op, uint32_t arg);

#endif
