/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * GPS L1 C/A acquisition for the M4 (gps_acq.h).
 */

#include "gps_acq.h"

#include "gps_nav.h"

#define N GPS_ACQ_N

/* gps_model.G2_DELAY */
static const uint16_t g2_gps[33] = {0,   5,   6,   7,   8,   17,  18,  139, 140, 141, 251,
				    252, 254, 255, 256, 257, 258, 469, 470, 471, 472, 473,
				    474, 509, 512, 513, 514, 515, 516, 859, 860, 861, 862};
static const uint16_t g2_sbas[19] = {145, 175, 52,  21,   237, 235, 886, 657, 634, 762,
				     355, 1012, 176, 603, 130, 359, 595, 68,  386};

/* The next output (stage 10) of a 10-stage generator, stage s in bit s - 1, taps a mask of stages,
 * stage 1 in bit 0; then the shift (gps_model.lfsr_out, one chip). */
static uint32_t lfsr_step(uint32_t* reg, uint32_t taps)
{
	const uint32_t out = (*reg >> 9) & 1;
	*reg = ((*reg << 1) | (uint32_t) __builtin_parity(*reg & taps)) & 0x3FF;
	return out;
}

#define G1_TAPS ((1u << 2) | (1u << 9))                                             /* 3, 10 */
#define G2_TAPS ((1u << 1) | (1u << 2) | (1u << 5) | (1u << 7) | (1u << 8) | (1u << 9)) /* 2,3,6,8,9,10 */

int gps_acq_code(int prn, uint8_t code[1023], uint32_t* g1v, uint32_t* g2v, uint32_t* win)
{
	int d;
	if ((prn >= 1) && (prn <= 32)) {
		d = g2_gps[prn];
	} else if ((prn >= 120) && (prn <= 138)) {
		d = g2_sbas[prn - 120];
	} else {
		return -1;
	}
	/* In code itself (the M4 app's stack is 4 KB): G2's output delayed by d, then G1's folded in.
	 * The generators after two shifts from chip 0 take the first ten delayed G2 chips. */
	uint32_t reg = 0x3FF;
	for (int k = 0; k < 1023; k++) {
		code[(k + d) % 1023] = (uint8_t) lfsr_step(&reg, G2_TAPS);
	}
	uint32_t g2 = 0;
	for (int i = 0; i < 10; i++) {
		g2 |= (uint32_t) code[9 - i] << i;
	}
	reg = 0x3FF;
	for (int k = 0; k < 1023; k++) {
		code[k] ^= (uint8_t) lfsr_step(&reg, G1_TAPS);
	}
	uint32_t g1 = 0x3FF;
	for (int s = 0; s < 2; s++) {
		const uint32_t f1 = ((g1 >> 2) ^ (g1 >> 9)) & 1;
		const uint32_t f2 = ((g2 >> 1) ^ (g2 >> 2) ^ (g2 >> 5) ^ (g2 >> 7) ^ (g2 >> 8) ^ (g2 >> 9)) & 1;
		g1 = ((g1 << 1) | f1) & 0x3FF;
		g2 = ((g2 << 1) | f2) & 0x3FF;
	}
	*g1v = g1;
	*g2v = g2;
	*win = code[1022] | (code[0] << 1) | (code[1] << 2);
	return 0;
}

/*
 * Radix 4, decimation in time: N = 4^5, five passes over memory where radix 2 made ten (radix 2
 * took 49 cycles a butterfly on the M4, 1.2 ms a transform). The stages take their input in base-4
 * digit-reversed order, which gps_acq_block, gps_acq_product and gps_acq_code_chips write directly
 * (on the M4 a separate reordering took 257 us a transform). Each stage's twiddles W^k, W^2k and
 * W^3k by rotation (one step constant a stage, exp(-2 pi i / len)): the app has no room for a table
 * (4 KB), and in single precision the rotation stays within about 1e-4 of the exact values.
 */
static float step_r[GPS_ACQ_STAGES] GPS_ACQ_STATE, step_i[GPS_ACQ_STAGES] GPS_ACQ_STATE;

void gps_acq_init(void)
{
	for (int s = 0; s < GPS_ACQ_STAGES; s++) {
		const double a = -2.0 * 3.14159265358979323846 / (double) (4 << (2 * s));
		step_r[s] = (float) gps_nav_cos(a);
		step_i[s] = (float) gps_nav_sin(a);
	}
}

/* The index after r in base-4 digit-reversed order, its top digit's unit top: from the reverse of
 * i, the reverse of i + 1 (a carry runs down from the top digit). */
static uint32_t rev_next(uint32_t r, uint32_t top)
{
	while (top && ((r & (3 * top)) == 3 * top)) {
		r -= 3 * top;
		top >>= 2;
	}
	return r + top;
}

/* Bin j = 4m + d lands at rev(m) + d N/4: rev(m), m's four digits reversed, steps by rev_next once
 * every four bins. */
#define REV_TOP (N >> 4)
#define REV_D   (N >> 2)

void gps_acq_cis(float x, float* c, float* s)
{
	const float x2 = x * x;
	*c = 1.0f - x2 * (0.5f - x2 * (1.0f / 24 - x2 * (1.0f / 720 - x2 * (1.0f / 40320))));
	*s = x * (1.0f - x2 * (1.0f / 6 - x2 * (1.0f / 120 - x2 * (1.0f / 5040))));
}

/* A rotator's magnitude back to one, every 64 steps. */
static void renorm(float* zr, float* zi)
{
	const float g = 1.5f - 0.5f * (*zr * *zr + *zi * *zi);
	*zr *= g;
	*zi *= g;
}

/* At -O2, the stages alone: the app builds at -Os, where they took 47 cycles a butterfly on the M4;
 * the whole file at -O2 cost 312 bytes more and gained nothing elsewhere, nor did gps_acq_block,
 * gps_acq_product and the accumulation alone (250 bytes). */
__attribute__((optimize("O2"))) void gps_acq_fft_stage(float* restrict x, int inverse, int s)
{
	/* Forward: (a1 - a3) times -i; inverse: conjugate twiddles and times +i. */
	const float sgn = inverse ? -1.0f : 1.0f;
	const int len = 4 << (2 * s), q = len >> 2;
	const float sr = step_r[s], si = sgn * step_i[s];
	float w1r = 1, w1i = 0;
	for (int k = 0; k < q; k++) {
		const float w2r = w1r * w1r - w1i * w1i, w2i = 2 * w1r * w1i;
		const float w3r = w2r * w1r - w2i * w1i, w3i = w2r * w1i + w2i * w1r;
		for (int i = k; i < N; i += len) {
			float* restrict p0 = &x[2 * i];
			float* restrict p1 = &x[2 * (i + q)];
			float* restrict p2 = &x[2 * (i + 2 * q)];
			float* restrict p3 = &x[2 * (i + 3 * q)];
			const float a0r = p0[0], a0i = p0[1];
			const float a1r = p1[0] * w1r - p1[1] * w1i, a1i = p1[0] * w1i + p1[1] * w1r;
			const float a2r = p2[0] * w2r - p2[1] * w2i, a2i = p2[0] * w2i + p2[1] * w2r;
			const float a3r = p3[0] * w3r - p3[1] * w3i, a3i = p3[0] * w3i + p3[1] * w3r;
			const float b0r = a0r + a2r, b0i = a0i + a2i;
			const float b1r = a0r - a2r, b1i = a0i - a2i;
			const float b2r = a1r + a3r, b2i = a1i + a3i;
			const float b3r = sgn * (a1i - a3i), b3i = -sgn * (a1r - a3r);
			p0[0] = b0r + b2r;
			p0[1] = b0i + b2i;
			p1[0] = b1r + b3r;
			p1[1] = b1i + b3i;
			p2[0] = b0r - b2r;
			p2[1] = b0i - b2i;
			p3[0] = b1r - b3r;
			p3[1] = b1i - b3i;
		}
		const float nr = w1r * sr - w1i * si;
		w1i = w1r * si + w1i * sr;
		w1r = nr;
	}
}

/* Index i with its base-4 digits reversed. */
static uint32_t digit_reverse(uint32_t i)
{
	uint32_t r = 0;
	for (int d = 0; d < GPS_ACQ_STAGES; d++) {
		r = (r << 2) | (i & 3);
		i >>= 2;
	}
	return r;
}

void gps_acq_fft(float* x, int inverse)
{
	for (uint32_t i = 0; i < N; i++) {
		const uint32_t j = digit_reverse(i);
		if (i < j) {
			const float r = x[2 * i], m = x[2 * i + 1];
			x[2 * i] = x[2 * j];
			x[2 * i + 1] = x[2 * j + 1];
			x[2 * j] = r;
			x[2 * j + 1] = m;
		}
	}
	for (int s = 0; s < GPS_ACQ_STAGES; s++) {
		gps_acq_fft_stage(x, inverse, s);
	}
}

void gps_acq_code_chips(const uint8_t code[1023], float* out)
{
	uint32_t r = 0;
	for (int j = 0; j < N; j += 4) {
		for (int d = 0; d < 4; d++) {
			/* the chip at the bin's middle: (j + 1/2) (BLOCK / N) (CHIP / FS) = (2j + 1) 1023 / 2048 */
			const int chip = (int) (((uint32_t) (2 * (j + d) + 1) * 1023) / 2048) % 1023;
			out[2 * (r + d * REV_D)] = code[chip] ? -1.0f : 1.0f;
			out[2 * (r + d * REV_D) + 1] = 0;
		}
		r = rev_next(r, REV_TOP);
	}
}

void gps_acq_code_round(const float* spectrum, int16_t* out)
{
	for (int j = 0; j < 2 * N; j++) { /* conjugate, rounded; |values| stay under 1,024 */
		const float v = (j & 1) ? -spectrum[j] : spectrum[j];
		out[j] = (int16_t) (v < 0 ? v - 0.5f : v + 0.5f);
	}
}

void gps_acq_code_spectrum(const uint8_t code[1023], float* scratch, int16_t* out)
{
	gps_acq_code_chips(code, scratch);
	for (int s = 0; s < GPS_ACQ_STAGES; s++) {
		gps_acq_fft_stage(scratch, 0, s);
	}
	gps_acq_code_round(scratch, out);
}

void gps_acq_block(
	const void* ring,
	uint32_t size,
	uint32_t first,
	float doppler,
	float* restrict out)
{
	/* The samples a word (two) at a time, I in each half's low byte. A word split by a bin's end
	 * keeps its second sample for the next bin; a block starting at a word's middle starts with
	 * that sample. */
	const uint32_t* w = (const uint32_t*) ring;
	const uint32_t mask = size / 4 - 1;
	uint32_t at = first / 4, k = 0, held = 0;
	int32_t hi = 0, hq = 0;
	if (first & 2) {
		const uint32_t v = w[at++ & mask];
		hi = (int8_t) (v >> 16);
		hq = (int8_t) (v >> 24);
		held = 1;
	}
	/* Mixed down by exp(-2 pi i doppler t) at each bin's middle sample, by a rotator. */
	float cr, ci, zr, zi;
	const float step = doppler * (float) (-2.0 * 3.14159265358979323846 * GPS_ACQ_BLOCK / N / GPS_ACQ_FS);
	gps_acq_cis(step, &cr, &ci);
	gps_acq_cis(0.5f * step, &zr, &zi);
	uint32_t r = 0;
	for (uint32_t j = 0; j < N; j += 4) {
		for (uint32_t d = 0; d < 4; d++) {
			/* Bin j + d takes samples ceil(j BLOCK / N) up to ceil((j + 1) BLOCK / N): 3 or 4. */
			const uint32_t end = ((j + d + 1) * GPS_ACQ_BLOCK + N - 1) / N;
			int32_t si = hi, sq = hq;
			k += held;
			held = 0;
			for (; k + 2 <= end; k += 2) {
				const uint32_t v = w[at++ & mask];
				si += (int8_t) v + (int8_t) (v >> 16);
				sq += (int8_t) (v >> 8) + (int8_t) (v >> 24);
			}
			hi = hq = 0;
			if (k < end) {
				const uint32_t v = w[at++ & mask];
				si += (int8_t) v;
				sq += (int8_t) (v >> 8);
				hi = (int8_t) (v >> 16);
				hq = (int8_t) (v >> 24);
				held = 1;
				k++;
			}
			const float fi = (float) si, fq = (float) sq;
			float* restrict o = &out[2 * (r + d * REV_D)];
			o[0] = fi * zr - fq * zi;
			o[1] = fi * zi + fq * zr;
			const float nr = zr * cr - zi * ci;
			zi = zr * ci + zi * cr;
			zr = nr;
		}
		if ((j & 63) == 60) {
			renorm(&zr, &zi);
		}
		r = rev_next(r, REV_TOP);
	}
}

/* cos and sin of pi / N: the half-bin row's rotator step. */
#define HALF_COS 0.99999529380957619f
#define HALF_SIN 0.0030679567629659761f

void gps_acq_product(
	const float* restrict spectrum,
	const int16_t* restrict code,
	float* restrict work,
	int half,
	int rot)
{
	/* r(j + 1/2) for the half-bin row: the product times exp(i pi k / N) for k < N/2 and, for the
	 * negative frequencies k - N, the same less a half turn; a rotator across k, its magnitude held
	 * at one. Written in the stages' order. */
	float zr = 1, zi = 0;
	uint32_t r = 0;
	for (uint32_t j = 0; j < N; j += 4) {
		if (half && (j == N / 2)) {
			zr = -zr;
			zi = -zi;
		}
		for (uint32_t d = 0; d < 4; d++) {
			const uint32_t k = (j + d + (uint32_t) rot) & (N - 1);
			const float ar = spectrum[2 * k], ai = spectrum[2 * k + 1];
			const float br = (float) code[2 * (j + d)], bi = (float) code[2 * (j + d) + 1];
			const float pr = ar * br - ai * bi, pi = ar * bi + ai * br;
			float* restrict o = &work[2 * (r + d * REV_D)];
			if (!half) {
				o[0] = pr;
				o[1] = pi;
				continue;
			}
			o[0] = pr * zr - pi * zi;
			o[1] = pr * zi + pi * zr;
			const float nr = zr * HALF_COS - zi * HALF_SIN;
			zi = zr * HALF_SIN + zi * HALF_COS;
			zr = nr;
		}
		if (half && ((j & 63) == 60)) {
			renorm(&zr, &zi);
		}
		r = rev_next(r, REV_TOP);
	}
}

/* work's |.|^2 into acc, row j to acc[j + shift] (0 <= shift <= N, modulo N), in two runs. */
static void accumulate(const float* restrict work, float* restrict acc, int shift)
{
	const int s = (shift == N) ? 0 : shift;
	float* restrict a = &acc[s];
	for (int j = 0; j < N - s; j++) {
		a[j] += work[2 * j] * work[2 * j] + work[2 * j + 1] * work[2 * j + 1];
	}
	const float* restrict v = &work[2 * (N - s)];
	for (int j = 0; j < s; j++) {
		acc[j] += v[2 * j] * v[2 * j] + v[2 * j + 1] * v[2 * j + 1];
	}
}

void gps_acq_accumulate(const float* work, float* acc, float* acc1, int shift, int half)
{
	/* Local cell 2j goes to 2j + shift: acc[j + s] for shift 2s, acc1[j + s] for 2s + 1; local
	 * cell 2j + 1 (the half-bin row) to acc1[j + s] or acc[j + s + 1]. */
	const int s = shift >> 1, odd = shift & 1;
	if (!acc1) {
		accumulate(work, acc, s + odd);
	} else if (!half) {
		accumulate(work, odd ? acc1 : acc, s);
	} else {
		accumulate(work, odd ? acc : acc1, s + odd);
	}
}

void gps_acq_correlate(
	const float* spectrum,
	const int16_t* code,
	float* work,
	float* acc,
	float* acc1,
	int shift,
	int half)
{
	gps_acq_product(spectrum, code, work, half, 0);
	for (int s = 0; s < GPS_ACQ_STAGES; s++) {
		gps_acq_fft_stage(work, 1, s);
	}
	gps_acq_accumulate(work, acc, acc1, shift, half);
}

/* since modulo the code period at doppler, in samples. The period is BLOCK / (1 + x), x = doppler /
 * L1 (1,023 chips at 4 Msps are BLOCK samples): whole nominal periods counted off in integers, then
 * their k * BLOCK * x / (1 + x) shortfall (at most about 120 samples) in single precision, which
 * keeps a sample's hundredth. */
static float offset(uint32_t since, float doppler)
{
	const float x = doppler * (float) (1.0 / GPS_ACQ_L1);
	const uint32_t k = since / GPS_ACQ_BLOCK;
	const float p = (float) GPS_ACQ_BLOCK / (1.0f + x);
	float r = (float) (since - k * GPS_ACQ_BLOCK) + (float) (k * GPS_ACQ_BLOCK) * x / (1.0f + x);
	while (r >= p) {
		r -= p;
	}
	while (r < 0) {
		r += p;
	}
	return r;
}

int gps_acq_shift(uint32_t since, float doppler)
{
	const float r = offset(since, doppler);
	const int h = (int) (r * (float) (2 * N) / (float) GPS_ACQ_BLOCK + 0.5f);
	return h % (2 * N);
}

uint32_t gps_acq_skip(uint32_t since, float doppler)
{
	/* A bin is BLOCK / N = 3.90625 samples: of the four starts, the nearest a bin's edge. */
	const float b = (float) GPS_ACQ_BLOCK / N;
	float u = offset(since, doppler);
	u -= (float) (int) (u / b) * b;
	uint32_t best = 0;
	float least = b;
	for (uint32_t d = 0; d < 4; d++) {
		float v = u + (float) d;
		v = (v >= b) ? v - b : v;
		v = (v > b - v) ? b - v : v;
		if (v < least) {
			least = v;
			best = d;
		}
	}
	return best;
}

float gps_acq_peak(const float* acc, const float* acc1, int* cell)
{
	int best = 0;
	float top = acc[0];
	for (int h = 1; h < 2 * N; h++) {
		if ((h & 1) && !acc1) {
			continue;
		}
		const float v = (h & 1) ? acc1[h >> 1] : acc[h >> 1];
		if (v > top) {
			top = v;
			best = h;
		}
	}
	const int bin = best >> 1;
	float noise = 0; /* single precision: on the M4 a double add is a library call */
	int count = 0;
	for (int j = 0; j < N; j++) {
		int d = j - bin;
		d = d < 0 ? -d : d;
		d = (d < N - d) ? d : N - d;
		if (d > 2) {
			noise += acc[j];
			count++;
		}
	}
	*cell = best;
	return top / (noise / (float) count);
}
