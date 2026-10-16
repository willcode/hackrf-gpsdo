/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * GPS navigation in C; see gps_nav.h.
 */

#include "gps_nav.h"

#include <stddef.h>

#define MU      3.986005e14
#define OMEGA_E 7.2921151467e-5
#define GPS_PI  3.1415926535898 /* the semicircle of IS-GPS-200, not pi */
#define REL_F   -4.442807633e-10
#define WEEK_S  604800.0

#define PI       3.14159265358979323846
#define PIO2_HI  1.57079632673412561417e+00 /* pi/2 in 33 bits: k * PIO2_HI is exact */
#define PIO2_LO  6.07710050650619224932e-11 /* pi/2 - PIO2_HI */
#define TAN_PI_8 0.41421356237309504880

/* ---- arithmetic */

static double floor_(double x)
{
	const double t = (double) (int64_t) x;
	return (t > x) ? t - 1.0 : t;
}

static double round_(double x)
{ return floor_(x + 0.5); }

double gps_nav_sqrt(double x)
{
	if (x <= 0) {
		return 0; /* a NaN goes on through: a solve on garbage then never reads as converged */
	}

	/*
	 * The FPU's single-precision root as the first guess (within 6e-8; the two conversions are soft
	 * but cheap), then Newton: each step doubles the bits, so two reach the last one (within 1 ulp
	 * of libm from 1e-12 to 1e30 on the host). Below 1e-45 or above 3e38 the guess is 0 or
	 * infinite and so is the result; nothing here takes such a root.
	 */
	double y = __builtin_sqrtf((float) x);
	for (int i = 0; i < 2; i++) {
		y = 0.5 * (y + x / y);
	}
	return y;
}

/*
 * sin and cos of r in [-pi/4, pi/4]: Taylor series to r^15 and r^16, error below 5e-17 (under the
 * last bit of a double near 0.7; measured within 1.1e-16 of the exact over the range), nested so
 * each term is the last times -r^2 / (n (n - 1)). The reciprocals are constants: on the M4 a division
 * costs several multiplications.
 */
static const double sin_k[7] = {
	1.0 / (15 * 14),
	1.0 / (13 * 12),
	1.0 / (11 * 10),
	1.0 / (9 * 8),
	1.0 / (7 * 6),
	1.0 / (5 * 4),
	1.0 / (3 * 2),
};
static const double cos_k[8] = {
	1.0 / (16 * 15),
	1.0 / (14 * 13),
	1.0 / (12 * 11),
	1.0 / (10 * 9),
	1.0 / (8 * 7),
	1.0 / (6 * 5),
	1.0 / (4 * 3),
	1.0 / (2 * 1),
};

static double sin_kernel(double r)
{
	const double r2 = r * r;
	double s = 0;
	for (int i = 0; i < 7; i++) {
		s = -(s + 1.0) * r2 * sin_k[i];
	}
	return r * (1.0 + s);
}

static double cos_kernel(double r)
{
	const double r2 = r * r;
	double c = 0;
	for (int i = 0; i < 8; i++) {
		c = -(c + 1.0) * r2 * cos_k[i];
	}
	return 1.0 + c;
}

/* x = k * pi/2 + r with |r| <= pi/4; returns r and k mod 4. Exact for |x| up to about 1.6e6 (k
 * within 20 bits); every argument here is under 100. */
static double reduce(double x, int* quadrant)
{
	const double k = round_(x * (2 / PI));
	*quadrant = (int) ((int64_t) k & 3);
	return (x - k * PIO2_HI) - k * PIO2_LO;
}

double gps_nav_sin(double x)
{
	int q;
	const double r = reduce(x, &q);
	switch (q) {
	case 0:
		return sin_kernel(r);
	case 1:
		return cos_kernel(r);
	case 2:
		return -sin_kernel(r);
	default:
		return -cos_kernel(r);
	}
}

double gps_nav_cos(double x)
{
	int q;
	const double r = reduce(x, &q);
	switch (q) {
	case 0:
		return cos_kernel(r);
	case 1:
		return -sin_kernel(r);
	case 2:
		return -cos_kernel(r);
	default:
		return sin_kernel(r);
	}
}

/* atan of t in [0, 1]. */
static double atan_unit(double t)
{
	double base = 0;
	if (t > TAN_PI_8) {
		/* atan t = pi/4 + atan((t - 1) / (t + 1)), and |(t - 1) / (t + 1)| <= tan(pi/8). */
		base = PI / 4;
		t = (t - 1) / (t + 1);
	}
	/* atan u = 2 atan(u / (1 + sqrt(1 + u^2))) halves the angle: |u| <= 0.199, then the series to
	 * u^25, error below 1e-19. */
	static const double k[12] = {
		1.0 / 25,
		1.0 / 23,
		1.0 / 21,
		1.0 / 19,
		1.0 / 17,
		1.0 / 15,
		1.0 / 13,
		1.0 / 11,
		1.0 / 9,
		1.0 / 7,
		1.0 / 5,
		1.0 / 3,
	};
	const double u = t / (1 + gps_nav_sqrt(1 + t * t));
	const double u2 = u * u;
	double s = 0;
	for (int i = 0; i < 12; i++) {
		s = u2 * (k[i] - s);
	}
	return base + 2 * u * (1 - s);
}

double gps_nav_atan2(double y, double x)
{
	const double ax = x < 0 ? -x : x;
	const double ay = y < 0 ? -y : y;
	if (ax == 0 && ay == 0) {
		return 0;
	}
	double a = (ay <= ax) ? atan_unit(ay / ax) : PI / 2 - atan_unit(ax / ay);
	if (x < 0) {
		a = PI - a;
	}
	return y < 0 ? -a : a;
}

/* ---- navigation data */

/* IS-GPS-200 Table 20-XIV: the data bits (D1 = bit 1) each parity bit D25-D30 covers. */
static const uint32_t parity_mask[6] = {
	0xEC7CD2, /* 1 2 3 5 6 10 11 12 13 14 17 18 20 23 */
	0x763E69, /* 2 3 4 6 7 11 12 13 14 15 18 19 21 24 */
	0xBB1F34, /* 1 3 4 5 7 8 12 13 14 15 16 19 20 22 */
	0x5D8F9A, /* 2 4 5 6 8 9 13 14 15 16 17 20 21 23 */
	0xAEC7CD, /* 1 3 5 6 7 9 10 14 15 16 17 18 21 22 24 */
	0x2DEA27, /* 3 5 6 8 9 10 11 13 15 19 22 23 24 */
};

static const uint8_t parity_prev[6] = {29, 30, 29, 30, 30, 29};

static int odd(uint32_t v)
{
	v ^= v >> 16;
	v ^= v >> 8;
	v ^= v >> 4;
	v ^= v >> 2;
	v ^= v >> 1;
	return (int) (v & 1);
}

int gps_nav_parity(uint32_t word, int d29, int d30, uint32_t* data)
{
	const uint32_t d = ((word >> 6) & 0xFFFFFF) ^ (d30 ? 0xFFFFFF : 0);
	for (int i = 0; i < 6; i++) {
		const int p =
			odd(d & parity_mask[i]) ^ (parity_prev[i] == 29 ? d29 : d30);
		if (p != (int) ((word >> (5 - i)) & 1)) {
			return 0;
		}
	}
	*data = d;
	return 1;
}

/* Bits a to b (from 1, most significant first) of 24-bit word w (from 1). */
static uint32_t bits(const uint32_t* ws, int w, int a, int b)
{ return (ws[w - 1] >> (24 - b)) & ((UINT32_C(1) << (b - a + 1)) - 1); }

static double sgn(uint32_t v, int n)
{ return (v >> (n - 1)) ? (double) v - (double) (UINT64_C(1) << n) : (double) v; }

/* 2^-n as a double, exactly. */
static double p2(int n)
{
	union {
		double d;
		uint64_t u;
	} v = {.u = (uint64_t) (1023 - n) << 52};

	return v.d;
}

/*
 * The ephemeris's double fields (IS-GPS-200 Table 20-III): subframe, index from tgd in struct
 * gps_eph, word, first bit, width (32: bits 17-24 of the word, then all of the next), and the
 * scale 2^-p, signed (S) and in semicircles (PI) as flagged. A table and one loop: as code, each
 * field's scaling was its own library calls.
 */
#define EPH_S  1
#define EPH_PI 2
_Static_assert(offsetof(struct gps_eph, idot) - offsetof(struct gps_eph, tgd) ==
		       20 * sizeof(double),
	       "the 21 doubles from tgd to idot stand in a row");
static const struct {
	uint8_t sf, idx, w, a, n, flags;
	int8_t p;
} eph_fields[] = {
	{1, 0, 7, 17, 8, EPH_S, 31},             /* tgd */
	{1, 1, 8, 9, 16, 0, -4},                 /* toc */
	{1, 2, 10, 1, 22, EPH_S, 31},            /* af0 */
	{1, 3, 9, 9, 16, EPH_S, 43},             /* af1 */
	{1, 4, 9, 1, 8, EPH_S, 55},              /* af2 */
	{2, 5, 3, 9, 16, EPH_S, 5},              /* crs */
	{2, 6, 4, 1, 16, EPH_S | EPH_PI, 43},    /* dn */
	{2, 7, 4, 17, 32, EPH_S | EPH_PI, 31},   /* m0 */
	{2, 8, 6, 1, 16, EPH_S, 29},             /* cuc */
	{2, 9, 6, 17, 32, 0, 33},                /* e */
	{2, 10, 8, 1, 16, EPH_S, 29},            /* cus */
	{2, 11, 8, 17, 32, 0, 19},               /* sqrta */
	{2, 12, 10, 1, 16, 0, -4},               /* toe */
	{3, 13, 3, 1, 16, EPH_S, 29},            /* cic */
	{3, 14, 3, 17, 32, EPH_S | EPH_PI, 31},  /* omega0 */
	{3, 15, 5, 1, 16, EPH_S, 29},            /* cis */
	{3, 16, 5, 17, 32, EPH_S | EPH_PI, 31},  /* i0 */
	{3, 17, 7, 1, 16, EPH_S, 5},             /* crc */
	{3, 18, 7, 17, 32, EPH_S | EPH_PI, 31},  /* omega */
	{3, 19, 9, 1, 24, EPH_S | EPH_PI, 43},   /* omegadot */
	{3, 20, 10, 9, 14, EPH_S | EPH_PI, 43},  /* idot */
};

int gps_nav_subframe(const uint32_t ws[10], struct gps_eph* e, uint32_t* tow)
{
	const int id = (int) bits(ws, 2, 20, 22);
	*tow = bits(ws, 2, 1, 17);
	if (id == 1) {
		e->week = (uint16_t) bits(ws, 3, 1, 10);
		e->health = (uint8_t) bits(ws, 3, 17, 22);
		e->iodc = (uint16_t) ((bits(ws, 3, 23, 24) << 8) | bits(ws, 8, 1, 8));
	} else if (id == 2) {
		e->iode2 = (uint8_t) bits(ws, 3, 1, 8);
	} else if (id == 3) {
		e->iode3 = (uint8_t) bits(ws, 10, 1, 8);
	} else {
		return id;
	}
	double* d = &e->tgd;
	for (unsigned i = 0; i < sizeof(eph_fields) / sizeof(eph_fields[0]); i++) {
		const unsigned w = eph_fields[i].w, a = eph_fields[i].a, n = eph_fields[i].n;
		if (eph_fields[i].sf != id) {
			continue;
		}
		const uint32_t v = (n == 32) ? (bits(ws, w, 17, 24) << 24) | bits(ws, w + 1, 1, 24)
					     : bits(ws, w, a, a + n - 1);
		double x = ((eph_fields[i].flags & EPH_S) ? sgn(v, n) : (double) v) * p2(eph_fields[i].p);
		if (eph_fields[i].flags & EPH_PI) {
			x *= GPS_PI;
		}
		d[eph_fields[i].idx] = x;
	}
	e->have |= 1u << (id - 1);
	return id;
}

int gps_nav_eph_complete(const struct gps_eph* e)
{
	return ((e->have & 7) == 7) && (e->iode2 == e->iode3) &&
		(e->iode2 == (e->iodc & 0xFF));
}

/* An almanac field: bits a to b of word w, signed when n is (n its width), times scale, in single
 * precision (the M4 has no double-precision hardware). */
static float alm_field(const uint32_t* ws, int w, int a, int b, int is_signed, float scale)
{
	const int n = b - a + 1;
	const uint32_t v = bits(ws, w, a, b);
	const int32_t x = is_signed ? (int32_t) (v << (32 - n)) >> (32 - n) : (int32_t) v;
	return (float) x * scale;
}

int gps_nav_almanac(const uint32_t ws[10], struct gps_alm alm[32])
{
	const int id = (int) bits(ws, 2, 20, 22), prn = (int) bits(ws, 3, 3, 8);
	if (((id != 4) && (id != 5)) || (bits(ws, 3, 1, 2) != 1) || (prn < 1) || (prn > 32)) {
		return 0;
	}
	struct gps_alm* a = &alm[prn - 1];
	const float sc23 = (float) (GPS_PI / 8388608.0); /* 2^-23 semicircles */
	a->e = alm_field(ws, 3, 9, 24, 0, (float) (1.0 / 2097152.0));
	a->toa = alm_field(ws, 4, 1, 8, 0, 4096.0f);
	a->i0 = (float) (0.3 * GPS_PI) + alm_field(ws, 4, 9, 24, 1, (float) (GPS_PI / 524288.0));
	a->omegadot = alm_field(ws, 5, 1, 16, 1, (float) (GPS_PI / 274877906944.0));
	a->health = (uint8_t) bits(ws, 5, 17, 24);
	a->sqrta = alm_field(ws, 6, 1, 24, 0, 1.0f / 2048.0f);
	a->omega0 = alm_field(ws, 7, 1, 24, 1, sc23);
	a->omega = alm_field(ws, 8, 1, 24, 1, sc23);
	a->m0 = alm_field(ws, 9, 1, 24, 1, sc23);
	a->have = 1;
	return prn;
}

void gps_nav_alm_eph(const struct gps_alm* a, struct gps_eph* e)
{
	*e = (struct gps_eph){0};
	e->toe = e->toc = a->toa;
	e->e = a->e;
	e->i0 = a->i0;
	e->omegadot = a->omegadot;
	e->sqrta = a->sqrta;
	e->omega0 = a->omega0;
	e->omega = a->omega;
	e->m0 = a->m0;
	e->health = a->health;
}

/* ---- orbits and the fix */

/* Into half a week either way; every use is within a week and a half of it. */
static double wrap(double t)
{
	return (t > WEEK_S / 2) ? t - WEEK_S : (t < -WEEK_S / 2) ? t + WEEK_S : t;
}

void gps_nav_sat_position(const struct gps_eph* e, double t, double pos[3], double* clock)
{
	const double a = e->sqrta * e->sqrta;
	const double tk = wrap(t - e->toe);
	const double n = gps_nav_sqrt(MU / (a * a * a)) + e->dn;
	const double m = e->m0 + n * tk;
	double ea = m;
	/* Until it stops changing, at most ten times (e below 0.03 for GPS: about ten for 1e-16). */
	for (int i = 0; i < 10; i++) {
		const double next = m + e->e * gps_nav_sin(ea);
		if (next == ea) {
			break;
		}
		ea = next;
	}
	const double se = gps_nav_sin(ea), ce = gps_nav_cos(ea);
	/* The true anomaly's sine and cosine directly, then the argument of latitude's by angle sums,
	 * its harmonic correction (under 1e-4 rad) to second order: no arctangent and four fewer
	 * kernels (within 2e-8 m of the arctangent's). */
	const double den = 1 - e->e * ce;
	const double sv = gps_nav_sqrt(1 - e->e * e->e) * se / den, cv = (ce - e->e) / den;
	const double sw = gps_nav_sin(e->omega), cw = gps_nav_cos(e->omega);
	const double sp = sv * cw + cv * sw, cp = cv * cw - sv * sw;
	const double s2 = 2 * sp * cp, c2 = (cp - sp) * (cp + sp);
	const double du = e->cus * s2 + e->cuc * c2;
	const double sdu = du * (1 - du * du * (1.0 / 6)), cdu = 1 - 0.5 * du * du;
	const double r = a * den + e->crs * s2 + e->crc * c2;
	const double i = e->i0 + e->idot * tk + e->cis * s2 + e->cic * c2;
	const double xp = r * (cp * cdu - sp * sdu), yp = r * (sp * cdu + cp * sdu);
	const double om = e->omega0 + (e->omegadot - OMEGA_E) * tk - OMEGA_E * e->toe;
	const double so = gps_nav_sin(om), co = gps_nav_cos(om);
	const double si = gps_nav_sin(i), ci = gps_nav_cos(i);
	pos[0] = xp * co - yp * ci * so;
	pos[1] = xp * so + yp * ci * co;
	pos[2] = yp * si;
	const double dt = wrap(t - e->toc);
	*clock = e->af0 + e->af1 * dt + e->af2 * dt * dt + REL_F * e->e * e->sqrta * se -
		e->tgd;
}

static double norm3(const double* v)
{ return gps_nav_sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }

/*
 * Sine and cosine of the WGS 84 latitude and longitude at p, by Bowring's closed form: square roots
 * and divisions, no trigonometry, within about 1e-10 rad of the iterated latitude near the ground.
 * On the M4 the usual iteration took over a millisecond.
 */
static void local_frame(
	const double p[3],
	double* sla,
	double* cla,
	double* slo,
	double* clo)
{
	const double a = 6378137.0, f = 1 / 298.257223563;
	const double b = a * (1 - f), e2 = f * (2 - f), ep2 = e2 / (1 - e2);
	const double rho = gps_nav_sqrt(p[0] * p[0] + p[1] * p[1]);
	*clo = rho > 0 ? p[0] / rho : 1;
	*slo = rho > 0 ? p[1] / rho : 0;
	const double u = p[2] * a, v = rho * b;
	const double t = gps_nav_sqrt(u * u + v * v);
	const double st = u / t, ct = v / t;
	const double nn = p[2] + ep2 * b * st * st * st;
	const double dd = rho - e2 * a * ct * ct * ct;
	const double r = gps_nav_sqrt(nn * nn + dd * dd);
	*sla = nn / r;
	*cla = dd / r;
}

/* The local up vector at rx. */
static void up_at(const double rx[3], double up[3])
{
	double sla, cla, slo, clo;
	local_frame(rx, &sla, &cla, &slo, &clo);
	up[0] = cla * clo;
	up[1] = cla * slo;
	up[2] = sla;
}

/* Solve the n x n system a x = b in place (Gaussian elimination, partial pivoting). 0 if singular. */
static int linsolve(int n, double a[4][4], double b[4])
{
	for (int c = 0; c < n; c++) {
		int p = c;
		for (int r = c + 1; r < n; r++) {
			const double v = a[r][c] < 0 ? -a[r][c] : a[r][c];
			const double w = a[p][c] < 0 ? -a[p][c] : a[p][c];
			if (v > w) {
				p = r;
			}
		}
		if (a[p][c] == 0) {
			return 0;
		}
		for (int k = 0; k < n; k++) {
			const double t = a[c][k];
			a[c][k] = a[p][k];
			a[p][k] = t;
		}
		const double t = b[c];
		b[c] = b[p];
		b[p] = t;
		for (int r = c + 1; r < n; r++) {
			const double f = a[r][c] / a[c][c];
			for (int k = c; k < n; k++) {
				a[r][k] -= f * a[c][k];
			}
			b[r] -= f * b[c];
		}
	}
	for (int c = n - 1; c >= 0; c--) {
		for (int k = c + 1; k < n; k++) {
			b[c] -= a[c][k] * b[k];
		}
		b[c] /= a[c][c];
	}
	return 1;
}

void gps_nav_solve_init(struct gps_solve* s)
{
	for (int i = 0; i < 4; i++) {
		s->x[i] = 0;
	}
	s->iteration = 0;
	s->done = 0;
	s->converged = 0;
	s->phase = 0;
	s->held = 0;
	s->height = 0;
}

enum { START, UP, FIRST, SECOND, SOLVE };

/* A row h with residual r into the normal equations (one copy: in the app, room counts). */
__attribute__((noinline)) static void add_row(struct gps_solve* s, const double h[4], double r)
{
	for (int a = 0; a < 4; a++) {
		s->atr[a] += h[a] * r;
		for (int b = 0; b < 4; b++) {
			s->ata[a][b] += h[a] * h[b];
		}
	}
}

int gps_nav_solve_step(
	struct gps_solve* s,
	int n,
	const double* t_tx,
	const struct gps_eph* const* eph,
	double t_rx)
{
	if (s->done) {
		return 0;
	}
	if (s->phase == START) {
		/* The normal equations of the linearized ranges, from zero. */
		s->tropo = (s->iteration > 2) && (norm3(s->x) > 6e6);
		for (int a = 0; a < 4; a++) {
			s->atr[a] = 0;
			for (int b = 0; b < 4; b++) {
				s->ata[a][b] = 0;
			}
		}
		s->k = 0;
		s->phase = s->tropo ? UP : s->iteration ? SECOND : FIRST;
	}
	switch (s->phase) {
	case UP:
		up_at(s->x, s->up);
		s->phase = SECOND; /* the troposphere from the third iteration: the positions are held */
		return GPS_NAV_STEP_UP;
	case FIRST: /* the first iteration only: the positions and clocks are held for the fix */
		gps_nav_sat_position(eph[s->k], t_tx[s->k], s->sv[s->k], &s->sv[s->k][3]);
		s->phase = SECOND;
		return GPS_NAV_STEP_FIRST;
	case SECOND: {
		const int k = s->k;
		double* const pos = s->sv[k];
		if (s->iteration == 0) {
			gps_nav_sat_position(eph[k], t_tx[k] - pos[3], pos, &pos[3]);
		}
		/* wrap: at the week's end one satellite's time can have rolled over and another's not */
		const double rho = GPS_NAV_C * wrap(t_rx - (t_tx[k] - pos[3]));
		double d[3];
		for (int j = 0; j < 3; j++) {
			d[j] = pos[j] - s->x[j];
		}
		/* The Earth turns during the flight (under 7e-6 rad: its sine and cosine to second order). */
		const double rot = OMEGA_E * norm3(d) / GPS_NAV_C;
		const double cr = 1 - 0.5 * rot * rot, sr = rot;
		const double sat[3] = {
			pos[0] * cr + pos[1] * sr,
			-pos[0] * sr + pos[1] * cr,
			pos[2]};
		for (int j = 0; j < 3; j++) {
			d[j] = sat[j] - s->x[j];
		}
		const double rng = norm3(d);
		double tropo = 0;
		if (s->tropo) {
			const double z =
				(d[0] * s->up[0] + d[1] * s->up[1] + d[2] * s->up[2]) /
				rng;
			tropo = 2.47 / (z + 0.0121); /* z is the sine of the elevation */
		}
		const double r = rho - (rng + s->x[3] + tropo);
		const double h[4] = {-d[0] / rng, -d[1] / rng, -d[2] / rng, 1.0};
		s->res[k] = r;
		add_row(s, h, r);
		s->k = k + 1;
		s->phase = (s->k >= n) ? SOLVE : s->iteration ? SECOND : FIRST;
		return GPS_NAV_STEP_SECOND;
	}
	default:
		break;
	}
	/* SOLVE */
	const int it = s->iteration;
	s->phase = START;
	if (s->held) {
		/* The height's row: the radius along u, the vertical from the Earth's center (at the
		 * center itself, the satellites' mean direction: the rows' sum is -ata[.][3]). In double:
		 * single precision's rounding (half a meter at the Earth's radius) kept the steps from
		 * settling, and one such fix in eight ended unconverged. */
		double u[3];
		const int away = norm3(s->x) > 1e6;
		for (int j = 0; j < 3; j++) {
			u[j] = away ? s->x[j] : -s->ata[j][3];
		}
		const double l = norm3(u);
		for (int j = 0; j < 3; j++) {
			u[j] /= l;
		}
		/* The ellipsoid's radius along u: a / sqrt(1 + e'^2 sin^2), e'^2 = 1/(1 - f)^2 - 1. */
		const double rad = 6378137.0 / gps_nav_sqrt(1 + 0.00673949674228 * u[2] * u[2]);
		const double h[4] = {u[0], u[1], u[2], 0};
		add_row(s, h, rad + s->height - (u[0] * s->x[0] + u[1] * s->x[1] + u[2] * s->x[2]));
	}
	for (int a = 0; a < 16; a++) {
		s->q[a / 4][a % 4] = s->ata[a / 4][a % 4]; /* for the DOP: the solve destroys ata */
	}
	if (!linsolve(4, s->ata, s->atr)) {
		s->done = 1;
		return GPS_NAV_STEP_SOLVE;
	}
	const double* dx = s->atr;
	for (int a = 0; a < 4; a++) {
		s->x[a] += dx[a];
	}
	s->iteration = it + 1;
	const double step = gps_nav_sqrt(
		dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2] + dx[3] * dx[3]);
	if ((it > 4) && (step < 1e-4)) {
		s->done = 1;
		s->converged = 1;
	} else if (s->iteration >= 10) {
		s->done = 1;
	}
	return GPS_NAV_STEP_SOLVE;
}

void gps_nav_solve_iteration(
	struct gps_solve* s,
	int n,
	const double* t_tx,
	const struct gps_eph* const* eph,
	double t_rx)
{
	const int it = s->iteration;
	while (!s->done && (s->iteration == it)) {
		gps_nav_solve_step(s, n, t_tx, eph, t_rx);
	}
}

int gps_nav_dop(const struct gps_solve* s, double dop[5])
{
	/* The inverse's diagonal, a column at a time. */
	double inv[4];
	for (int c = 0; c < 4; c++) {
		double a[4][4], b[4] = {0, 0, 0, 0};
		for (int i = 0; i < 16; i++) {
			a[i / 4][i % 4] = s->q[i / 4][i % 4];
		}
		b[c] = 1;
		if (!linsolve(4, a, b) || !(b[c] > 0)) {
			return 0; /* singular, or so near it that rounding gave a nonpositive variance */
		}
		inv[c] = b[c];
	}
	const double tr3 = inv[0] + inv[1] + inv[2];
	dop[0] = gps_nav_sqrt(tr3 + inv[3]);
	dop[1] = gps_nav_sqrt(tr3);
	dop[2] = 0; /* H and V need the local frame's rotation: left out for room in the app */
	dop[3] = 0;
	dop[4] = gps_nav_sqrt(inv[3]);
	return 1;
}
