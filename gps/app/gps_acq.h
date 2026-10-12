/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * GPS L1 C/A acquisition for the M4, from the GPS image's snapshots: the C form of the host's
 * search (gps/tools/hackrf_gps.c, gps_model.acquire_all) cut to the app's memory. Portable between
 * the GPS app and the host, where test_gps_acq.py checks it against numpy and the recordings.
 *
 * A block is one code period of samples (4,000 at 4 Msps) box-car summed into GPS_ACQ_N bins of
 * about a chip each, mixed by the Doppler under test. Its FFT times the code's conjugate spectrum,
 * transformed back, is the circular correlation at every bin; |.|^2 is summed over blocks, each
 * block's row shifted by its start's offset from the first block's modulo the code period, so the
 * cells line up. The metric is the best cell over the mean of the cells away from it.
 *
 * Single-precision float: the M4 has a single-precision FPU. No C library.
 */

#ifndef GPS_ACQ_H
#define GPS_ACQ_H

#include <stdint.h>

/* Where the app keeps the FFT's twiddle steps: its data region (apps/Makefile), not after its image,
 * whose region is the tighter. Empty on the host. */
#ifndef GPS_ACQ_STATE
#define GPS_ACQ_STATE
#endif

#define GPS_ACQ_N     1024   /* bins per block: one per 3.90625 samples, about a chip */
#define GPS_ACQ_LOG2N 10
#define GPS_ACQ_BLOCK 4000   /* samples per block at 4 Msps: one nominal code period */
#define GPS_ACQ_FS    4e6
#define GPS_ACQ_CHIP  1.023e6
#define GPS_ACQ_L1    1575.42e6

/* The C/A code of a GPS PRN (1-32) or SBAS PRN (120-138) as chips 0/1, and the engine's arm
 * values: G1 and G2 after the first two chips, and the window of chips 1022, 0 and 1. Returns 0,
 * or -1 for an unknown PRN. */
int gps_acq_code(int prn, uint8_t code[1023], uint32_t* g1, uint32_t* g2, uint32_t* win);

/* cos and sin of x in single precision, |x| <= 0.5 (within 1e-7): the double routines are library
 * calls on the M4. */
void gps_acq_cis(float x, float* c, float* s);

/* The FFT's stage constants; once before any transform. */
void gps_acq_init(void);

/* In-place FFT of GPS_ACQ_N complex values (re, im interleaved) in natural order; inverse
 * unscaled. */
void gps_acq_fft(float* x, int inverse);

/* Its radix-4 stages alone, s from 0 to GPS_ACQ_STAGES - 1 in turn, on values already in the
 * stages' order (base-4 digit-reversed), as gps_acq_block, gps_acq_product and gps_acq_code_chips
 * write them: for a caller that spreads a transform over calls, about 0.1 ms each on the M4. */
#define GPS_ACQ_STAGES (GPS_ACQ_LOG2N / 2)
void gps_acq_fft_stage(float* x, int inverse, int s);

/*
 * The code's conjugate spectrum at the block's bins, chip phase 0 at the block's first sample (bin
 * j holds the chip at the middle of its samples), rounded to integers: 4 KB rather than 8. scratch
 * holds 2 N floats.
 */
void gps_acq_code_spectrum(const uint8_t code[1023], float* scratch, int16_t* out);

/* The same in parts: the chips at the bins into out (2 N floats, in the stages' order), then each
 * gps_acq_fft_stage of it, then the conjugate rounded into out. */
void gps_acq_code_chips(const uint8_t code[1023], float* out);
void gps_acq_code_round(const float* spectrum, int16_t* out);

/*
 * One block for the FFT's stages: GPS_ACQ_BLOCK samples from byte first (even) of a ring of size
 * bytes (a power of two, word aligned), each sample a 16-bit word with I in its low byte and Q in
 * its high (signed), summed into the bins, mixed down by doppler (the phase at each bin's middle
 * sample; |doppler| up to 80 kHz), written to out (2 N floats) in the stages' order.
 */
void gps_acq_block(const void* ring, uint32_t size, uint32_t first, float doppler, float* out);

/*
 * A block's spectrum (gps_acq_block, then the stages) times a code's conjugate spectrum, transformed
 * back into work, |.|^2 added into the accumulators: with half 0 the correlation at whole bins,
 * with half 1 at half bins (the code half a bin later, by a phase ramp across the spectrum), one
 * FFT each so each fits a call. acc[j] holds bin j, acc1[j] the cell between bins j and j + 1;
 * shift (half bins, from gps_acq_shift) aligns the block with the first, an odd one trading the
 * rows. With acc1 0, whole bins alone: an odd shift goes to the next bin (gps_acq_skip makes it
 * even). At one bin a chip, a code between cells loses up to 6 dB (2.3 dB on average); with the
 * half bins, up to 2.5 dB (1.1 dB on average), against the host's search at 4 Msps.
 */
void gps_acq_correlate(
	const float* spectrum,
	const int16_t* code,
	float* work,
	float* acc,
	float* acc1,
	int shift,
	int half);

/* gps_acq_correlate in parts: the product into work (in the stages' order), then each
 * gps_acq_fft_stage (inverse) of work, then the accumulation. The product takes the spectrum rot
 * bins on: the block's at its Doppler plus rot kHz, as a block's bins are a code period, 1 ms. */
void gps_acq_product(const float* spectrum, const int16_t* code, float* work, int half, int rot);
void gps_acq_accumulate(const float* work, float* acc, float* acc1, int shift, int half);

/* The shift that aligns a block starting since samples after the first with it, at doppler: since
 * modulo the code period at that Doppler, in half bins, rounded. */
int gps_acq_shift(uint32_t since, float doppler);

/* The samples, 0 to 3, a block starting since samples after the first skips to start within half
 * a sample of one of the first's bins: its shift then even, its cells the first's whole bins. A
 * snapshot's 4,096 samples hold a block of 4,000 from any of them. */
uint32_t gps_acq_skip(uint32_t since, float doppler);

/*
 * The best cell of acc, and of acc1 when given, and its metric: the cell over the mean of acc's
 * cells more than two bins away. *cell counts half bins: 2j for acc[j], 2j + 1 for acc1[j].
 */
float gps_acq_peak(const float* acc, const float* acc1, int* cell);

#endif
