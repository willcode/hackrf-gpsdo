/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * The GPS/SBAS tracking loops for one correlator channel, in integer arithmetic: the C form of
 * fpga-build/dsp/gps_loop_int.py (IntLoops), word for word. Portable between the GPS app's M4 code
 * and the host, where test_gps_loops.py checks it against the Python model.
 *
 * One call per dump record (one code period): the six sums in, the channel's next code and carrier
 * words out. The words go to the correlator as a pending request; they take effect at the period
 * after next.
 *
 * Requires arithmetic right shifts of negative values (as GCC gives on the host and the M4): the model
 * floors.
 */

#ifndef GPS_LOOPS_H
#define GPS_LOOPS_H

#include <stdint.h>

#define GPS_LOOPS_MAX_SYMBOL_MS 20

struct gps_dump {
	int64_t ie, qe, ip, qp, il, ql;
};

/* Loop bandwidths in hertz, and the symbol-sum PLL's limit on bandwidth times symbol time. */
struct gps_bandwidths {
	double fll, sync, assist, dll, pll, symbol_bt;
};

/* The defaults of gps_loop_int.gains(). */
extern const struct gps_bandwidths gps_bandwidths_default;

struct gps_gains {
	int32_t fll, sync, pll_i, pll_p, assist, dll;
};

struct gps_loops {
	struct gps_gains k;
	double fs;
	int64_t w_nom;
	int64_t freq;  /* carrier word, Q16 */
	int64_t integ; /* PLL integrator, carrier word, Q16 */
	int64_t corr;  /* DLL code-word correction */
	int symbol_ms, fll_ms, sync_ms;
	uint32_t count;
	int have_prev;
	int32_t prev_angle;
	uint32_t flips[GPS_LOOPS_MAX_SYMBOL_MS];
	int edge; /* -1 until found */
	int have_acc;
	struct gps_dump acc;
	int64_t half[2][2];
	int64_t code, carrier; /* the words last returned */
};

/* The integer gains for a sample rate and symbol length; 0 when one does not fit 16 bits. */
int gps_gains(struct gps_gains* k, double fs, int symbol_ms, const struct gps_bandwidths* bw);

/* Start a channel at the acquisition's Doppler. symbol_ms is 20 for GPS, 2 for SBAS; fll_ms and
 * sync_ms are the stages' lengths (300 and 1000 in the model). Returns 0, or -1 for bad arguments. */
int gps_loops_init(
	struct gps_loops* l,
	double fs,
	double doppler,
	int symbol_ms,
	int fll_ms,
	int sync_ms,
	const struct gps_bandwidths* bw);

/* One dump in; the words to send next in *code and *carrier. */
void gps_loops_update(struct gps_loops* l, const struct gps_dump* d, int64_t* code, int64_t* carrier);

/*
 * n records missed (a snapshot pauses the stream; not in IntLoops, which never sees a gap): their
 * periods counted, so the symbol's position stays on the record count, the symbol in progress
 * dropped, and no frequency step taken across the gap.
 */
void gps_loops_skip(struct gps_loops* l, uint32_t n);

/* The loops' Doppler estimate in hertz. */
double gps_loops_doppler(const struct gps_loops* l);

#endif
