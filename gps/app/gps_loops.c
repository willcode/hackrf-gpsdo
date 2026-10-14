/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * See gps_loops.h. Each function names the gps_loop_int.py function it reproduces.
 */

#include "gps_loops.h"

#include <stddef.h>


#define ANG_BITS     16
#define FULL         (1 << ANG_BITS)
#define HALF         (FULL >> 1)
#define QUARTER      (FULL >> 2)
#define CORDIC_ITERS 16
#define PRE_SHIFT    2
#define Q            16
#define AID          10894 /* round(2^24 / 1540) */
#define CHIP_RATE    1.023e6
#define ASSIST_MIN_MS 10

/* round(atan(2^-i) / (2 pi) * 2^16) */
static const int32_t ATAN[CORDIC_ITERS] =
	{8192, 4836, 2555, 1297, 651, 326, 163, 81, 41, 20, 10, 5, 3, 1, 1, 0};

const struct gps_bandwidths gps_bandwidths_default = {
	.fll = 8.0,
	.sync = 3.0,
	.assist = 1.0,
	.dll = 1.0,
	.pll = 15.0,
	.symbol_bt = 0.16,
};

/* wrap() */
static int32_t wrap(int64_t a)
{
	return (int32_t) (((a + HALF) & (FULL - 1)) - HALF);
}

/* fold() */
static int32_t fold(int32_t a)
{
	if (a > QUARTER) {
		return a - HALF;
	}
	if (a <= -QUARTER) {
		return a + HALF;
	}
	return a;
}

/* cordic(): the angle, and the length times the CORDIC gain in *mag. */
static int32_t cordic(int64_t x, int64_t y, int64_t* mag)
{
	int64_t a = 0;
	x >>= PRE_SHIFT;
	y >>= PRE_SHIFT;
	if (x < 0) {
		x = -x;
		y = -y;
		a = HALF;
	}
	for (int i = 0; i < CORDIC_ITERS; i++) {
		const int64_t xs = x >> i, ys = y >> i;
		if (y > 0) {
			x += ys;
			y -= xs;
			a += ATAN[i];
		} else {
			x -= ys;
			y += xs;
			a -= ATAN[i];
		}
	}
	if (mag) {
		*mag = x;
	}
	return wrap(a);
}

static int64_t abs64(int64_t v)
{
	return v < 0 ? -v : v;
}

/*
 * div_q14(), for |n| <= |d| (the DLL's early less late over early plus late): the quotient is at
 * most 2^14, so 15 steps of long division give it exactly, with no 64-bit division (on the M4 a
 * library routine of about 960 bytes).
 */
static int64_t div_q14(int64_t n, int64_t d)
{
	if (d == 0) {
		return 0;
	}
	uint64_t r = (uint64_t) abs64(n);
	const uint64_t b = (uint64_t) abs64(d);
	int64_t q = 0;
	for (int i = 0; i < 15; i++) {
		q <<= 1;
		if (r >= b) {
			r -= b;
			q |= 1;
		}
		r <<= 1;
	}
	return ((n >= 0) == (d > 0)) ? q : -q;
}

/*
 * Round to nearest, half to even, as Python's round does; |v| < 2^51. Adding 1.5 * 2^52 makes the
 * addition's own rounding do it and leaves 2^51 + the result in the low 52 bits: one library call on
 * the M4, where a truncation, a difference and the compares took seven.
 */
static int64_t round_even(double v)
{
	const union {
		double d;
		uint64_t u;
	} x = {.d = v + 6755399441055744.0};

	return (int64_t) (x.u & 0xFFFFFFFFFFFFFull) - ((int64_t) 1 << 51);
}

/* gains() */
int gps_gains(struct gps_gains* k, double fs, int symbol_ms, const struct gps_bandwidths* bw)
{
	const double w = 4294967296.0 / fs;
	const double T = symbol_ms * 1e-3;
	const double lim = bw->symbol_bt / T;
	const double w0 = (bw->pll < lim ? bw->pll : lim) / 0.53;
	const int64_t v[6] = {
		round_even(4 * bw->fll * w),
		round_even(4 * bw->sync * w),
		round_even(w0 * w0 * T * w),
		round_even(1.414 * w0 * w),
		round_even(8 * bw->assist * w),
		round_even(4 * bw->dll * w),
	};
	for (int i = 0; i < 6; i++) {
		if (v[i] < 0 || v[i] >= (1 << 16)) {
			return 0;
		}
	}
	k->fll = v[0];
	k->sync = v[1];
	k->pll_i = v[2];
	k->pll_p = v[3];
	k->assist = v[4];
	k->dll = v[5];
	return 1;
}

int gps_loops_init(
	struct gps_loops* l,
	double fs,
	double doppler,
	int symbol_ms,
	int fll_ms,
	int sync_ms,
	const struct gps_bandwidths* bw)
{
	if (symbol_ms < 1 || symbol_ms > GPS_LOOPS_MAX_SYMBOL_MS || fs <= 0) {
		return -1;
	}
	*l = (struct gps_loops) {0};
	if (!gps_gains(&l->k, fs, symbol_ms, bw ? bw : &gps_bandwidths_default)) {
		return -1;
	}
	l->fs = fs;
	l->w_nom = round_even(CHIP_RATE / fs * 4294967296.0);
	l->freq = round_even(doppler / fs * 4294967296.0) * ((int64_t) 1 << Q);
	l->integ = l->freq;
	l->symbol_ms = symbol_ms;
	l->fll_ms = fll_ms;
	l->sync_ms = sync_ms;
	l->edge = -1;
	return 0;
}

double gps_loops_doppler(const struct gps_loops* l)
{
	return (double) (l->freq >> Q) * l->fs / 4294967296.0;
}

/* IntLoops._dll() */
static void dll(struct gps_loops* l, const struct gps_dump* d)
{
	int64_t me, ml;
	cordic(d->ie, d->qe, &me);
	cordic(d->il, d->ql, &ml);
	const int64_t tau = div_q14(me - ml, me + ml);
	l->corr = ((int64_t) l->k.dll * tau) >> 15;
}

/* IntLoops._words() */
static void words(struct gps_loops* l)
{
	const int64_t carrier = l->freq >> Q;
	l->code = l->w_nom + ((carrier * AID) >> 24) + l->corr;
	l->carrier = carrier;
}

void gps_loops_skip(struct gps_loops* l, uint32_t n)
{
	l->count += n;
	l->have_prev = 0;
	l->have_acc = 0;
}

/* IntLoops.update() */
void gps_loops_update(struct gps_loops* l, const struct gps_dump* d, int64_t* code, int64_t* carrier)
{
	const uint32_t idx = l->count++;
	if (l->edge < 0) {
		const int32_t angle = cordic(d->ip, d->qp, NULL);
		if (l->have_prev) {
			const int32_t diff = wrap((int64_t) angle - l->prev_angle);
			const int64_t gain = idx < (uint32_t) l->fll_ms ? l->k.fll : l->k.sync;
			l->freq += gain * fold(diff);
			l->integ = l->freq;
			if (idx >= (uint32_t) l->fll_ms && !(-QUARTER <= diff && diff <= QUARTER)) {
				l->flips[idx % l->symbol_ms]++;
			}
		}
		l->have_prev = 1;
		l->prev_angle = angle;
		if (idx + 1 >= (uint32_t) (l->fll_ms + l->sync_ms)) {
			/* The most flips; the earliest on a tie. */
			int best = 0;
			for (int i = 1; i < l->symbol_ms; i++) {
				if (l->flips[i] > l->flips[best]) {
					best = i;
				}
			}
			l->edge = best;
		}
		dll(l, d);
		words(l);
	} else {
		const int pos = (int) ((idx - (uint32_t) l->edge) % (uint32_t) l->symbol_ms); /* idx > edge */
		if (pos == 0) {
			l->acc = (struct gps_dump) {0};
			l->half[0][0] = l->half[0][1] = l->half[1][0] = l->half[1][1] = 0;
			l->have_acc = 1;
		}
		if (l->have_acc) {
			l->acc.ie += d->ie;
			l->acc.qe += d->qe;
			l->acc.ip += d->ip;
			l->acc.qp += d->qp;
			l->acc.il += d->il;
			l->acc.ql += d->ql;
			int64_t* h = l->half[pos >= l->symbol_ms / 2];
			h[0] += d->ip;
			h[1] += d->qp;
			if (pos == l->symbol_ms - 1) {
				if (l->symbol_ms >= ASSIST_MIN_MS) {
					const int32_t a0 = cordic(l->half[0][0], l->half[0][1], NULL);
					const int32_t a1 = cordic(l->half[1][0], l->half[1][1], NULL);
					l->integ += (int64_t) l->k.assist * wrap((int64_t) a1 - a0);
				}
				const int32_t e = fold(cordic(l->acc.ip, l->acc.qp, NULL));
				l->integ += (int64_t) l->k.pll_i * e;
				l->freq = l->integ + (int64_t) l->k.pll_p * e;
				dll(l, &l->acc);
				l->have_acc = 0;
				words(l);
			}
		}
	}
	*code = l->code;
	*carrier = l->carrier;
}
