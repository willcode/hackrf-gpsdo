/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * The oscillator's frequency offset from the satellites' carriers (gps_freq.h).
 */

#include "gps_freq.h"

#define OMEGA_E 7.2921151467e-5
#define TWO32   4294967296.0

void gps_freq_open(struct gps_freq_window* w, uint32_t n0, uint32_t wc)
{
	w->n_first = n0;
	w->n_last = n0;
	w->wc_last = wc;
	w->cycles = 0;
	w->records = 1;
	w->open = 1;
	w->bad = 0;
}

void gps_freq_record(struct gps_freq_window* w, uint32_t n0, uint32_t wc)
{
	const uint32_t dn = n0 - w->n_last;
	if (dn > GPS_FREQ_GAP) {
		w->bad = 1;
	}
	/* Under 2^14 samples times 2^31 (a longer gap marks the window bad): 45 bits a record, under
	 * 56 bits over a second. */
	w->cycles += (int64_t) dn * (int32_t) w->wc_last;
	w->n_last = n0;
	w->wc_last = wc;
	w->records++;
}

double gps_freq_measured(const struct gps_freq_window* w, double lo_word, double fs)
{
	const double span = (double) (w->n_last - w->n_first);
	return ((double) w->cycles + span * lo_word) / span * fs / TWO32;
}

void gps_freq_end(
	const struct gps_eph* eph,
	double t,
	const double rx[3],
	int end,
	double* off,
	double* clk)
{
	const double s = end ? GPS_FREQ_HALF : -GPS_FREQ_HALF;
	double pos[3];
	gps_nav_sat_position(eph, t + s, pos, clk);
	double d[3] = {pos[0] - rx[0], pos[1] - rx[1], pos[2] - rx[2]};
	const double rot = OMEGA_E * gps_nav_sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) / GPS_NAV_C;
	const double c = 1 - 0.5 * rot * rot, sn = rot; /* rot under 7e-6 rad: to second order */
	d[0] = pos[0] * c + pos[1] * sn - rx[0];
	d[1] = -pos[0] * sn + pos[1] * c - rx[1];
	*off = s + gps_nav_sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) / GPS_NAV_C;
}

double gps_freq_expected(const double off[2], const double clk[2])
{
	return GPS_FREQ_L1 * ((2 * GPS_FREQ_HALF + clk[1] - clk[0]) / (off[1] - off[0]) - 1.0);
}

float gps_freq_predict(const float* t, const float* e, int n, float at)
{
	if (n < 3) {
		return e[n - 1];
	}
	float st = 0, se = 0;
	for (int i = 0; i < n; i++) {
		st += t[i];
		se += e[i];
	}
	const float tm = st / (float) n, em = se / (float) n;
	float stt = 0, ste = 0;
	for (int i = 0; i < n; i++) {
		stt += (t[i] - tm) * (t[i] - tm);
		ste += (t[i] - tm) * (e[i] - em);
	}
	return em + ste / stt * (at - tm);
}

void gps_freq_fraction(double x, uint32_t max_c, uint32_t* b, uint32_t* c)
{
	/* Convergents p/q of x's continued fraction until q passes max_c; then the best of the
	 * last convergent and the largest semiconvergent that fits. */
	double p0 = 0, q0 = 1, p1 = 1, q1 = 0, v = x;
	for (int i = 0; i < 64; i++) {
		const double a = (double) (uint64_t) v; /* v >= 0 */
		const double p2 = a * p1 + p0, q2 = a * q1 + q0;
		if (q2 > max_c) {
			const double k = (double) (uint64_t) ((max_c - q0) / q1);
			const double ps = k * p1 + p0, qs = k * q1 + q0;
			const double es = x - ps / qs, e1 = x - p1 / q1;
			if ((qs > 0) && ((es < 0 ? -es : es) < (e1 < 0 ? -e1 : e1))) {
				p1 = ps;
				q1 = qs;
			}
			break;
		}
		p0 = p1;
		q0 = q1;
		p1 = p2;
		q1 = q2;
		const double r = v - a;
		if (r < 1e-300) {
			break;
		}
		v = 1.0 / r;
		if (v > 4e18) {
			break;
		}
	}
	*b = (uint32_t) p1;
	*c = (uint32_t) q1; /* 1/1 when x rounds to 1: the caller takes the next integer */
}

/* Neither NaN nor infinite: an exponent not all ones (bits, not a library call, on the M4). */
static int finite(double x)
{
	const union {
		double d;
		uint64_t u;
	} v = {.d = x};

	return ((v.u >> 52) & 0x7FF) != 0x7FF;
}

static double weight(double var)
{
	const double v = (var == 0) ? GPS_FREQ_PRIOR_VAR : var;
	return 1.0 / ((v < GPS_FREQ_FLOOR * GPS_FREQ_FLOOR) ? GPS_FREQ_FLOOR * GPS_FREQ_FLOOR : v);
}

void gps_freq_var_update(double* var, double r, double reject)
{
	if (!finite(r)) {
		return; /* a NaN would stay in the variance, and its weight make every estimate NaN */
	}
	const double a = __builtin_fabs(r), c = (a > reject) ? reject : a; /* only c * c counts */
	const double v = (*var == 0) ? GPS_FREQ_PRIOR_VAR : *var;
	*var = v + (c * c - v) / GPS_FREQ_VAR_K;
}

uint32_t gps_freq_combine(const double* ppb, const double* var, int n, double reject, double* est)
{
	double v[32];
	int m = 0;
	if (n > 32) {
		return 0;
	}
	for (int i = 0; i < n; i++) { /* insertion sort of the finite offsets, for the median */
		if (!finite(ppb[i])) {
			continue; /* never within reject of the median below either */
		}
		int j = m++;
		for (; (j > 0) && (v[j - 1] > ppb[i]); j--) {
			v[j] = v[j - 1];
		}
		v[j] = ppb[i];
	}
	if (m < 1) {
		return 0;
	}
	const double med = (m & 1) ? v[m / 2] : 0.5 * (v[m / 2 - 1] + v[m / 2]);
	uint32_t used = 0;
	double sum = 0, wsum = 0;
	for (int i = 0; i < n; i++) {
		if (__builtin_fabs(ppb[i] - med) <= reject) {
			const double w = weight(var[i]);
			used |= 1u << i;
			sum += w * ppb[i];
			wsum += w;
		}
	}
	if (used) {
		*est = sum / wsum;
	}
	return used;
}
