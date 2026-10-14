/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * The oscillator's frequency offset from each satellite's carrier, against a held position: the C
 * form of gps/freq.py's expected() and combine(), with the measurement's accumulator. Portable
 * between the GPS app's M4 code and the host, where test_gps_app.py's replays of device streams check
 * the app's offsets against freq.py's.
 *
 * Per window of records and satellite: the carrier cycles the engine applied (each record's
 * carrier word times the samples to the next record's start, summed exactly in 64 bits) plus the
 * tuner's offset give the received carrier frequency. The ephemeris and the held position give the
 * frequency a perfect oscillator would see. Measured less expected is -L1 times the oscillator's
 * fractional error, common to every satellite while the tuner and the sample clock share it.
 */

#ifndef GPS_FREQ_H
#define GPS_FREQ_H

#include <stdint.h>

#include "gps_nav.h"

#define GPS_FREQ_L1 1575.42e6
#define GPS_FREQ_HALF 0.5 /* half the transmit-time interval of the expected frequency's difference, s */
/*
 * Samples between record starts beyond which a window is not used: two records missing (a period
 * is about 4,000). A snapshot pauses the stream for a record or two; the cycles across the gap are
 * the last carrier word's, which the loops move by millihertz in that time.
 */
#define GPS_FREQ_GAP 12500u

/* One satellite's window of records: from the record opening it to the latest. */
struct gps_freq_window {
	uint32_t n_first, n_last; /* period starts of the first and latest record */
	uint32_t wc_last;         /* the latest record's carrier word, applied until the next */
	int64_t cycles;           /* carrier word times samples, summed: cycles * 2^32 */
	uint32_t records;
	uint8_t open, bad; /* bad: a gap over GPS_FREQ_GAP */
};

void gps_freq_open(struct gps_freq_window* w, uint32_t n0, uint32_t wc);

/* A record starting at n0 with carrier word wc. */
void gps_freq_record(struct gps_freq_window* w, uint32_t n0, uint32_t wc);

/* The carrier frequency (Hz) over the window, the tuner's offset lo_word (carrier word units)
 * added; fs the sample rate. */
double gps_freq_measured(const struct gps_freq_window* w, double lo_word, double fs);

/*
 * One end of the expected frequency's difference (end 0 at t - H, 1 at t + H, t a transmit time):
 * the receive time less t (s), from the satellite's range to rx with the Earth's rotation during
 * the flight, and the satellite clock correction (s). One satellite position: on the M4 a few
 * hundred microseconds, so a caller takes one end per tick.
 */
void gps_freq_end(
	const struct gps_eph* eph,
	double t,
	const double rx[3],
	int end,
	double* off,
	double* clk);

/* The expected carrier frequency (Hz) from the two ends. */
double gps_freq_expected(const double off[2], const double clk[2]);

/*
 * Each satellite's weight is the inverse of its offset's recent variance about the estimate
 * (ppb^2), which gps_freq_var_update keeps: GPS_FREQ_PRIOR_VAR before any history, never under
 * GPS_FREQ_FLOOR^2, so one satellite's weight is bounded. A satellite near the lock gate can sit
 * 0.5 ppb from the others while they agree to 0.01-0.03 ppb; weighted, it moves the estimate by
 * thousandths.
 */
#define GPS_FREQ_PRIOR_VAR (0.1 * 0.1)
#define GPS_FREQ_FLOOR     0.01
#define GPS_FREQ_VAR_K     16 /* the variance follows the squared residual with weight 1/16 */

/*
 * The common offset over n satellites' offsets (ppb): the weighted mean of those within reject of
 * their median, var[i] each one's variance (0 for none yet). Returns the mask of those used (bit i
 * for ppb[i]); n at most 32. A NaN or infinite offset (an ephemeris that gives no expected
 * frequency) takes no part. None is within reject of the median only when they disagree (two
 * more than twice reject apart) or none is finite: 0 then, and *est is left alone.
 */
uint32_t gps_freq_combine(const double* ppb, const double* var, int n, double reject, double* est);

/* A satellite's variance after a window whose offset lay r from the estimate (r limited to
 * +-reject, as for one left out); unchanged when r is NaN or infinite. */
void gps_freq_var_update(double* var, double r, double reject);

/*
 * The steering loop's arithmetic (the app holds its state; freq.py's Servo is the reference).
 *
 * gps_freq_predict: the least-squares line through n estimates e[i] at times t[i] (seconds),
 * evaluated at time at; the latest estimate when n < 3. Single precision: keep the times small
 * (the app counts them back from the latest estimate).
 *
 * gps_freq_fraction: the fraction b/c nearest x (0 <= x < 1) with c at most max_c, by continued
 * fractions with the last step's semiconvergents: 0/1 for x nearest 0, 1/1 for x nearest 1.
 */
float gps_freq_predict(const float* t, const float* e, int n, float at);
void gps_freq_fraction(double x, uint32_t max_c, uint32_t* b, uint32_t* c);

#endif
