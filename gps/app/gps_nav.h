/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * GPS navigation in C: subframe parity and decoding, the broadcast ephemeris, satellite position and
 * clock, and the least-squares fix (IS-GPS-200 20.3). The C form of gps/nav.py's decode_subframe,
 * sat_position, solve and dop, computed the same way. Portable between the GPS app's M4 code and the
 * host, where test_gps_nav.py checks it against nav.py.
 *
 * No C library: the M4 app links none, so the square root, sine, cosine and arctangent are here, in
 * double precision. The M4's FPU is single precision, so on the M4 every double operation is a
 * library call: the solve is split into iterations a caller can spread over ticks.
 */

#ifndef GPS_NAV_H
#define GPS_NAV_H

#include <stdint.h>

#define GPS_NAV_C      299792458.0
#define GPS_NAV_MAX_SV 12

/* A satellite's broadcast clock and orbit, scaled to seconds, meters and radians. */
struct gps_eph {
	uint8_t have;  /* bit n-1: subframe n seen */
	uint16_t week; /* mod 1024 */
	uint8_t health;
	uint16_t iodc;
	uint8_t iode2, iode3;
	double tgd, toc, af0, af1, af2;
	double crs, dn, m0, cuc, e, cus, sqrta, toe;
	double cic, omega0, cis, i0, crc, omega, omegadot, idot;
};

/* The fix in progress: x, y, z (m, ECEF) and the clock bias (m). */
struct gps_solve {
	double x[4];
	double res[GPS_NAV_MAX_SV]; /* each satellite's residual (m) in the last iteration */
	int iteration;
	int done;
	int converged; /* done because the step fell below 0.1 mm, not at the tenth iteration */
	/* Within an iteration: the step, the satellite, and what the steps carry between them. */
	int phase, k, tropo;
	double up[3];
	double sv[GPS_NAV_MAX_SV][4]; /* each satellite's position (m, ECEF) and clock (s): the first
				       * iteration's, held for the rest */
	double ata[4][4], atr[4];
	double q[4][4]; /* the last iteration's normal matrix (with the held height's row), for the DOP */
	/* Set after gps_nav_solve_init for a fix from three satellites: the height (m, above the WGS 84
	 * ellipsoid) held as a fourth equation, and the start in x[0..2] (zero: the Earth's center). */
	int held;
	double height;
};

double gps_nav_sqrt(double x);
double gps_nav_sin(double x);
double gps_nav_cos(double x);
double gps_nav_atan2(double y, double x);

/*
 * One 30-bit word, D1 in bit 29, with the last two bits of the word before it. Returns 1 and the 24
 * data bits (D1 in bit 23) with the D30* inversion removed when the six parity bits match, else 0.
 */
int gps_nav_parity(uint32_t word, int d29, int d30, uint32_t* data);

/*
 * A subframe's ten words of 24 data bits into the ephemeris. Returns the subframe ID (1-5) and its
 * TOW count (the HOW's 17 bits: the GPS time of the next subframe's start, in 6 s units).
 */
int gps_nav_subframe(const uint32_t words[10], struct gps_eph* eph, uint32_t* tow);

/*
 * The orbit a GPS satellite can have: its semi-major axis's square root (m^0.5; 5,153.6 nominal)
 * and its eccentricity. Every field the subframes carry is bounded by its width; within these the
 * orbit's radius is 24,000-29,000 km and every position and frequency finite. A spoofed ephemeris
 * with sqrt(A) = 0 made the expected frequency NaN.
 */
#define GPS_NAV_SQRTA_MIN 5000.0
#define GPS_NAV_SQRTA_MAX 5300.0
#define GPS_NAV_E_MAX     0.03

/* Subframes 1-3 present, from one issue of data, and their orbit within the bounds above. */
int gps_nav_eph_complete(const struct gps_eph* eph);

/*
 * An almanac entry's orbit (IS-GPS-200 20.3.3.5.1.2) in single precision: an almanac places a
 * satellite within a few kilometers, far coarser than a float's resolution at orbit radius
 * (meters). i0 includes the 0.3 semicircles the broadcast offset is from. Its clock terms are left
 * out: the entry serves the Doppler prediction, which af1 moves by about 0.01 Hz.
 */
struct gps_alm {
	float toa, e, i0, omegadot, sqrta, omega0, omega, m0;
	uint8_t health;
	uint8_t have; /* 1 from a page here; the GPS app's warm start also writes 2 and 3 (gps_app.c) */
};

/* A subframe 4 or 5 page's almanac into alm[prn - 1]. Returns the PRN (1-32), or 0 for a page that
 * holds none (data ID not 1, or an SV ID outside 1-32) or whose orbit is outside the ephemeris's
 * bounds (alm[prn - 1] then unchanged). */
int gps_nav_almanac(const uint32_t words[10], struct gps_alm alm[32]);

/* An almanac entry as an ephemeris for gps_nav_sat_position: the terms it lacks 0 (the clock
 * terms among them), its reference time for toe and toc. */
void gps_nav_alm_eph(const struct gps_alm* alm, struct gps_eph* eph);

/* ECEF position (m) at GPS time of week t and the satellite clock correction (s). */
void gps_nav_sat_position(
	const struct gps_eph* eph,
	double t,
	double pos[3],
	double* clock);

/*
 * The fix from transmit times at one receive instant: t_tx[i] is the GPS time of week at which the
 * signal received at t_rx left satellite i; t_rx is the receiver's time of that instant, any guess
 * within tens of milliseconds. Start with gps_nav_solve_init, then call gps_nav_solve_step until
 * s->done (at most ten iterations). gps_nav_solve_init leaves s->sv, the last solve's satellite
 * positions, as they were: the app's held-height start reads them. The bias is s->x[3] / C, true
 * time = t_rx - bias.
 * With s->held, three satellites suffice: the height above the ellipsoid is held at s->height (along
 * the radius, the ellipsoid's own within about 0.2 degree). Three spheres and the ellipsoid meet
 * twice; the solve finds the one nearer its start, and from the Earth's center the wrong one about
 * one time in twenty (a nominal sky, satellites above 10 degrees).
 *
 * A step is one satellite position (the first iteration's two a satellite; later iterations hold
 * them), a satellite's row, the local vertical, or an iteration's solution: on the M4 a position
 * costs a few hundred microseconds, the first iteration with seven satellites some milliseconds,
 * a later one far less. gps_nav_solve_iteration takes the steps of one iteration.
 */
void gps_nav_solve_init(struct gps_solve* s);

/* The step taken, or 0 when the fix is done. */
enum {
	GPS_NAV_STEP_UP = 1,     /* the local vertical, for the tropospheric delay */
	GPS_NAV_STEP_FIRST = 2,  /* a satellite's position at its transmit time */
	GPS_NAV_STEP_SECOND = 3, /* its row (the first iteration: its position again, less its clock) */
	GPS_NAV_STEP_SOLVE = 4,  /* the iteration's solution */
};

int gps_nav_solve_step(
	struct gps_solve* s,
	int n,
	const double* t_tx,
	const struct gps_eph* const* eph,
	double t_rx);
void gps_nav_solve_iteration(
	struct gps_solve* s,
	int n,
	const double* t_tx,
	const struct gps_eph* const* eph,
	double t_rx);

/* G, P and T dilution of precision of a finished solve, from its last iteration's normal matrix
 * (with a held height's row), in dop[0], [1] and [4] (H and V, [2] and [3], are 0). Returns 0 if
 * singular. */
int gps_nav_dop(const struct gps_solve* s, double dop[5]);

#endif
