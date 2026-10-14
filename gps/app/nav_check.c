/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Host driver for gps_nav.c, for test_gps_nav.py. One command per line on stdin, one answer line per
 * command on stdout, numbers with 17 significant digits:
 *   math x y                    sin x, cos x, atan2(y, x), sqrt |x|
 *   parity word d29 d30         ok and the 24 data bits (word and data in hex)
 *   sub prn w1 ... w10          a subframe's 24-bit data words (hex) into prn's ephemeris: id and TOW
 *   eph prn                     complete, week, health, iodc, iode2, iode3, then the 21 doubles
 *   alm w1 ... w10              a subframe 4 or 5 page: the PRN (0: none), then its entry: health,
 *                               toa, e, i0, omegadot, sqrta, omega0, omega, m0
 *   almpos prn t                the PRN's satellite position (m) at t from its entry
 *   solve t_rx n prn t_tx ...   iterations, x, y, z, bias (m), the residuals, then the five DOPs
 *                               at the fix for the satellites' positions at t_tx
 *   freq ppb x y z t n m e ...  a frequency window as gps_app.c's freq_finish takes it (gps_freq.c),
 *                               channel i's carrier from ephemeris m and expected from e, the
 *                               oscillator ppb off: the estimate (nan for none), the mask used, then
 *                               each channel's variance, kept from window to window
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gps_freq.h"
#include "gps_nav.h"

#define MAX_PRN 64

static struct gps_eph ephs[MAX_PRN];
static struct gps_alm alm[32];
static double fvar[GPS_NAV_MAX_SV];

/* The expected carrier frequency (Hz) from eph at transmit time t, received at rx. */
static double expected(const struct gps_eph* eph, double t, const double rx[3])
{
	double off[2], clk[2];
	for (int end = 0; end < 2; end++) {
		gps_freq_end(eph, t, rx, end, &off[end], &clk[end]);
	}
	return gps_freq_expected(off, clk);
}

/* The line's next word, "" past its end. */
static const char* next(char** save)
{
	const char* t = strtok_r(NULL, " \n", save);
	return t ? t : "";
}

/* A whole number from lo to hi, or -1. */
static int index_in(const char* s, int lo, int hi)
{
	char* end;
	const long v = strtol(s, &end, 10);
	return ((end == s) || (v < lo) || (v > hi)) ? -1 : (int) v;
}

int main(void)
{
	char line[4096];
	while (fgets(line, sizeof(line), stdin)) {
		char* save = NULL;
		const char* cmd = strtok_r(line, " \n", &save);
		if (!cmd) {
			continue;
		}
#define NEXT() next(&save)
		if (!strcmp(cmd, "math")) {
			const double x = strtod(NEXT(), NULL);
			const double y = strtod(NEXT(), NULL);
			printf("%.17g %.17g %.17g %.17g\n",
			       gps_nav_sin(x),
			       gps_nav_cos(x),
			       gps_nav_atan2(y, x),
			       gps_nav_sqrt(x < 0 ? -x : x));
		} else if (!strcmp(cmd, "parity")) {
			const uint32_t w = (uint32_t) strtoul(NEXT(), NULL, 16);
			const int d29 = atoi(NEXT());
			const int d30 = atoi(NEXT());
			uint32_t d = 0;
			const int ok = gps_nav_parity(w, d29, d30, &d);
			printf("%d %06x\n", ok, (unsigned) d);
		} else if (!strcmp(cmd, "sub")) {
			const int prn = index_in(NEXT(), 0, MAX_PRN - 1);
			uint32_t ws[10], tow;
			if (prn < 0) {
				printf("?\n");
				fflush(stdout);
				continue;
			}
			for (int i = 0; i < 10; i++) {
				ws[i] = (uint32_t) strtoul(NEXT(), NULL, 16);
			}
			const int id = gps_nav_subframe(ws, &ephs[prn], &tow);
			printf("%d %u\n", id, (unsigned) tow);
		} else if (!strcmp(cmd, "alm")) {
			uint32_t ws[10];
			for (int i = 0; i < 10; i++) {
				ws[i] = (uint32_t) strtoul(NEXT(), NULL, 16);
			}
			const int prn = gps_nav_almanac(ws, alm);
			printf("%d", prn);
			if (prn) {
				const struct gps_alm* a = &alm[prn - 1];
				printf(" %u %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g", a->health, a->toa, a->e,
				       a->i0, a->omegadot, a->sqrta, a->omega0, a->omega, a->m0);
			}
			printf("\n");
		} else if (!strcmp(cmd, "almpos")) {
			const int prn = index_in(NEXT(), 1, 32);
			const double t = strtod(NEXT(), NULL);
			struct gps_eph e;
			if (prn < 0) {
				printf("?\n");
				fflush(stdout);
				continue;
			}
			double pos[3], clk;
			gps_nav_alm_eph(&alm[prn - 1], &e);
			gps_nav_sat_position(&e, t, pos, &clk);
			printf("%.17g %.17g %.17g\n", pos[0], pos[1], pos[2]);
		} else if (!strcmp(cmd, "eph")) {
			const int prn = index_in(NEXT(), 0, MAX_PRN - 1);
			if (prn < 0) {
				printf("?\n");
				fflush(stdout);
				continue;
			}
			const struct gps_eph* e = &ephs[prn];
			printf("%d %u %u %u %u %u",
			       gps_nav_eph_complete(e),
			       e->week,
			       e->health,
			       e->iodc,
			       e->iode2,
			       e->iode3);
			const double v[] = {e->tgd, e->toc, e->af0,    e->af1,
					    e->af2, e->crs, e->dn,     e->m0,
					    e->cuc, e->e,   e->cus,    e->sqrta,
					    e->toe, e->cic, e->omega0, e->cis,
					    e->i0,  e->crc, e->omega,  e->omegadot,
					    e->idot};
			for (unsigned i = 0; i < sizeof(v) / sizeof(v[0]); i++) {
				printf(" %.17g", v[i]);
			}
			printf("\n");
		} else if (!strcmp(cmd, "solve")) {
			const double t_rx = strtod(NEXT(), NULL);
			const int n = index_in(NEXT(), 1, GPS_NAV_MAX_SV);
			double t_tx[GPS_NAV_MAX_SV];
			const struct gps_eph* e[GPS_NAV_MAX_SV];
			int bad = n < 0;
			for (int k = 0; k < n; k++) {
				const int prn = index_in(NEXT(), 0, MAX_PRN - 1);
				bad |= prn < 0;
				e[k] = &ephs[prn < 0 ? 0 : prn];
				t_tx[k] = strtod(NEXT(), NULL);
			}
			if (bad) {
				printf("?\n");
				fflush(stdout);
				continue;
			}
			struct gps_solve s;
			gps_nav_solve_init(&s);
			while (!s.done) {
				gps_nav_solve_iteration(&s, n, t_tx, e, t_rx);
			}
			printf("%d", s.iteration);
			for (int i = 0; i < 4; i++) {
				printf(" %.17g", s.x[i]);
			}
			for (int k = 0; k < n; k++) {
				printf(" %.17g", s.res[k]);
			}
			double dop[5] = {0};
			gps_nav_dop(&s, dop);
			for (int i = 0; i < 5; i++) {
				printf(" %.17g", dop[i]);
			}
			printf("\n");
		} else if (!strcmp(cmd, "freq")) {
			const double ppb_true = strtod(NEXT(), NULL);
			double rx[3];
			for (int i = 0; i < 3; i++) {
				rx[i] = strtod(NEXT(), NULL);
			}
			const double t = strtod(NEXT(), NULL);
			const int n = index_in(NEXT(), 1, GPS_NAV_MAX_SV);
			double ppb[GPS_NAV_MAX_SV];
			int bad = n < 0;
			for (int k = 0; k < n; k++) {
				const int m = index_in(NEXT(), 0, MAX_PRN - 1);
				const int e = index_in(NEXT(), 0, MAX_PRN - 1);
				bad |= (m < 0) || (e < 0);
				if (!bad) {
					const double f = expected(&ephs[m], t, rx) - GPS_FREQ_L1 * ppb_true * 1e-9;
					ppb[k] = -(f - expected(&ephs[e], t, rx)) / GPS_FREQ_L1 * 1e9;
				}
			}
			if (bad) {
				printf("?\n");
				fflush(stdout);
				continue;
			}
			double est = 0;
			const uint32_t used = gps_freq_combine(ppb, fvar, n, 1.0, &est);
			if (!used || ((n > 1) && !(used & (used - 1)))) {
				printf("nan 0");
			} else {
				for (int k = 0; k < n; k++) {
					gps_freq_var_update(&fvar[k], ppb[k] - est, 1.0);
				}
				printf("%.17g %u", est, (unsigned) used);
			}
			for (int k = 0; k < n; k++) {
				printf(" %.17g", fvar[k]);
			}
			printf("\n");
		} else {
			printf("?\n");
		}
		fflush(stdout);
	}
	return 0;
}
