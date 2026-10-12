/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Feeds a stream of dumps through gps_loops for test_gps_loops.py.
 *
 * Input, one line each: "fs doppler symbol_ms fll_ms sync_ms", then "IE QE IP QP IL QL" per dump.
 * Output: "fll sync pll_i pll_p assist dll" (the gains), then "code carrier edge" per dump.
 */

#include "gps_loops.h"

#include <inttypes.h>
#include <stdio.h>

int main(void)
{
	double fs, doppler;
	int symbol_ms, fll_ms, sync_ms;
	if (scanf("%lf %lf %d %d %d", &fs, &doppler, &symbol_ms, &fll_ms, &sync_ms) != 5) {
		return 2;
	}
	struct gps_loops l;
	if (gps_loops_init(&l, fs, doppler, symbol_ms, fll_ms, sync_ms, NULL) != 0) {
		return 3;
	}
	printf("%d %d %d %d %d %d\n", l.k.fll, l.k.sync, l.k.pll_i, l.k.pll_p, l.k.assist, l.k.dll);
	struct gps_dump d;
	while (scanf("%" SCNd64 " %" SCNd64 " %" SCNd64 " %" SCNd64 " %" SCNd64 " %" SCNd64,
		     &d.ie, &d.qe, &d.ip, &d.qp, &d.il, &d.ql) == 6) {
		int64_t code, carrier;
		gps_loops_update(&l, &d, &code, &carrier);
		printf("%" PRId64 " %" PRId64 " %d\n", code, carrier, l.edge);
	}
	return 0;
}
