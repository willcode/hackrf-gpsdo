/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * The GPS app's M4 side: the tracking loops for the correlator in the GPS
 * image (fpga-build/top/gps_rx.py, dsp/gps_engine.py).
 *
 * Each tick reads the record stream from the sample buffer, where the M0
 * writes it while the image streams records (setting 0 = 1) and a receive
 * session runs. Each record goes through its channel's loops
 * (gps/app/gps_loops.c, the integer model word for word); new code and
 * carrier words go back to the engine as a pending request, written into
 * FPGA register 6 byte by byte through the firmware's register service.
 *
 * The host arms channels and asks for snapshots through messages; it finds
 * satellites in the snapshots. Messages (little endian), each answered with
 * one status byte (0 done, 1 busy, 2 refused) unless noted:
 *   'M' v                setting 0 (output mode)
 *   'S' v                setting 1: a snapshot of (v + 1) * 1024 samples
 *   'A' ch prn W(4) Wc(4) G1(2) G2(2) win start(4) doppler(4, mHz)
 *                        arm a channel; start is the sample before its first
 *   'X' ch               stop a channel
 *   'Q'                  status (see status_reply)
 *   'L'                  the main loop's last full second (see load_reply)
 *   'P'                  navigation: the last fix and each channel's data (see nav_reply)
 *   'F' offset(4) [q(4)] the tuner's offset from 1575.42 MHz, in microhertz (signed); and the
 *                        quarter shift, the digital center less the analog, in hertz (signed)
 *   'B' a(4) b(4) c(4)   set PLL B's feedback ratio to a + b/c (nominal 32; see steer)
 *   'J'                  the last PLL B writes (see steer_reply)
 *   'D' on(1)            the steering loop on (1) or off (0, the ratio then held; see servo)
 *   'C' blocks refine    search for satellites and arm them (see acq_search); blocks per cold
 *                        (PRN, Doppler) and per candidate's check (at least the cold's); and from then the satellite
 *                        manager (mgr_second); 'C' 0 0 the manager alone
 *   'K'                  the search's state (see acq_reply)
 *   'R'                  the oscillator's frequency offset from the carriers (see freq_reply)
 *   'U' N(4)             survey the held position over N kept fixes (default 60)
 *   'U' x(8) y(8) z(8)   hold this ECEF position (doubles, meters) instead
 *   'U' x y z N(4)       survey over N kept fixes from this rough position (a grid square), held
 *                        until the first kept fix
 *   'H' h(4) [rough(1)]  the height (float, meters above the WGS 84 ellipsoid) held by a fix from
 *                        three satellites (default 0); rough 1: such a fix may replace the held
 *                        position given (a grid square's center)
 *   'W' ...              the warm start: what a run learned, read out and handed back (see warm)
 *
 * Navigation (gps/app/gps_nav.c): each GPS channel sums its prompt over the loops' symbol edge into
 * bits, finds subframes by their TLM and HOW, and keeps its ephemeris and the record at which its
 * latest subframe began. Once a second of samples the app takes every channel's transmit time at
 * one sample index and solves for position and clock, one iteration a tick.
 *
 * Frequency (gps/app/gps_freq.c): each GPS channel sums the carrier cycles the engine applied over
 * windows of a second of samples on a common grid, and gates each window on its symbol-sum phase
 * lock. Against the held position (the mean of the first kept fixes, or one given) and the
 * ephemeris, each satellite's carrier less its expected frequency gives the oscillator's offset;
 * the satellites' combined offset is the estimate. The steering loop (message 'D', servo) sets
 * PLL B, which clocks CLKOUT, against it.
 */

#include <stdint.h>

#include "app.h"
#include "gps_acq.h"
#include "gps_freq.h"
#include "gps_loops.h"
#include "gps_nav.h"

#define CHANNELS     7
#define RECORD_WORDS 21
#define MAGIC        0xD5
#define FS           4e6
#define MAX_WORDS    1024 /* per tick */
#define OUT_SIZE     512  /* command bytes waiting for the register service */
#define CHIP_RATE    1.023e6
#define TWO32        4294967296.0
#define FIX_SAMPLES  4000000u /* a fix each second of samples */
#define MAX_GDOP     50.0     /* a fix with a worse geometry is not kept */
#define SURVEY_FIXES 60       /* kept fixes averaged into the held position */
/*
 * A satellite whose offset is this far from the window's median is left out: locked satellites
 * agree to 0.02-0.07 ppb a second on the recordings, so 1 ppb is many times their scatter and far
 * under what a lost lock gives. gps/freq.py's REJECT_PPB.
 */
#define REJECT_PPB 1.0
#define MIN_LOCK   0.5 /* symbol-sum phase lock under which a window is not used (locked 0.90-0.98) */

enum {
	OP_ARM = 0,
	OP_PENDING = 1,
	OP_STOP = 2,
	OP_SETTING = 3,
};

/* A GPS channel's navigation data. */
struct chnav {
	uint64_t sr; /* the last 64 bits, the newest in bit 0 */
	uint32_t bits;
	int64_t acc, accq; /* the prompt summed over the bit so far */
	uint8_t have_acc;
	uint8_t resume; /* a gap crossed a bit's end: the next record starts the bit after it */
	uint8_t nwords; /* words of the subframe being read; 0 while searching */
	uint8_t wbits;  /* bits of the word being read */
	uint8_t anchored;
	uint32_t words[10];
	uint32_t sf_record; /* the record that opens the subframe's first period */
	uint32_t anchor_record;
	double anchor_t; /* GPS time of week at which that period began at the satellite */
	uint16_t subframes;
	uint32_t tow; /* the latest HOW's count */
	struct gps_eph eph;
	/*
	 * The latest two records: index, and their periods' start index, code phase and code word. A
	 * transmit time is taken within the period that holds its sample, with the code word the
	 * engine used there: the DLL moves the word at each period, and extrapolating across a
	 * boundary is off by up to 0.2 m.
	 */
	struct {
		uint32_t rec, n0, F, W;
	} last[2]; /* [0] the latest; both there once the channel's records count two */
};

struct channel {
	uint8_t active;
	uint8_t prn;
	uint8_t send;        /* words changed and not yet queued */
	uint8_t outstanding; /* a pending request sent, its words not yet seen in a record */
	uint8_t waited;      /* records since it was sent */
	uint32_t sent_W, sent_Wc;
	int32_t sent_steer;  /* the steering word in the request outstanding */
	int32_t steer_in_use; /* the steering word in the words the records show */
	double fvar;          /* the frequency offset's recent variance, ppb^2 (gps_freq.h) */
	uint32_t lost; /* pending requests never seen applied */
	uint32_t records;
	uint32_t gaps; /* records missing by their start indices */
	uint32_t last_n;
	int32_t prompt; /* |IP| + |QP| of the last record */
	uint16_t age;        /* seconds since the arm, up to MGR_SETTLE (the manager's) */
	uint8_t unlocked;    /* seconds in a row under the manager's lock */
	/* The loops' Doppler in the last second it was locked (the arm's until then): unlocked, the
	 * channel is kept while its loops stay near it. */
	float fd_lock;
	uint8_t odd;         /* records in a row whose period start jumped too far to be one */
	int64_t code, carrier;
	struct gps_loops loops;
	struct chnav nav;
	/* The frequency window open since the grid point before next_b, the next grid point. */
	struct gps_freq_window fw;
	uint32_t next_b;
	uint8_t have_b;
	uint8_t fw_hold; /* a PLL B write: the next window opened is not used either */
	double lock_i, lock_q; /* the window's whole symbols' prompt sums, squared and summed */
	struct {
		uint32_t end; /* the grid point the window ended at */
		uint8_t valid;
		float lock;
		uint32_t span;
		double f;     /* carrier frequency, Hz */
		double t_end; /* transmit time at the closing record */
	} done;
};

static struct {
	uint32_t pos; /* bytes of the sample buffer consumed */
	uint8_t synced;
	uint8_t in_record;
	uint8_t rec_n;
	uint16_t rec[RECORD_WORDS];
	uint8_t snap_pending;
	uint8_t stopping, stopped; /* channels stopped this second and the second before (bits) */
	uint32_t snap_words; /* words of the snapshot asked for */
	uint32_t snap_skip;  /* words of a snapshot still to pass */
	uint32_t snaps, records, garbage, overruns, refused;
	uint8_t toggle;
	uint8_t out[OUT_SIZE];
	uint32_t out_head, out_tail;
	struct channel ch[CHANNELS];
} g APP_STATE;

/* ---- commands to the engine (dsp/gps_cmd.py) ------------------------------------------ */

static uint32_t out_room(void)
{ return OUT_SIZE - (g.out_head - g.out_tail); }

static void out_byte(uint32_t start, uint32_t data)
{
	g.toggle ^= 1;
	g.out[g.out_head % OUT_SIZE] = (g.toggle << 7) | (start << 6) | (data & 0x3F);
	g.out_head++;
}

/* A command: the start byte, then the payload bits six at a time, least significant first. */
static int command(uint32_t op, uint32_t ch, const uint8_t* bits, uint32_t nbytes)
{
	if (out_room() < nbytes + 1) {
		return 0;
	}
	out_byte(1, (op << 4) | ch);
	for (uint32_t i = 0; i < nbytes; i++) {
		uint32_t v = 0;
		for (uint32_t b = 0; b < 6; b++) {
			const uint32_t k = 6 * i + b;
			v |= ((bits[k / 8] >> (k % 8)) & 1u) << b;
		}
		out_byte(0, v);
	}
	return 1;
}

static void put_bits(uint8_t* bits, uint32_t at, uint64_t v, uint32_t width)
{
	for (uint32_t i = 0; i < width; i++) {
		if ((v >> i) & 1) {
			bits[(at + i) / 8] |= 1u << ((at + i) % 8);
		}
	}
}

static int send_pending(uint32_t ch, uint32_t W, uint32_t Wc)
{
	uint8_t bits[9] = {0};
	put_bits(bits, 0, W, 32);
	put_bits(bits, 32, Wc, 32);
	return command(OP_PENDING, ch, bits, 11);
}

static int send_setting(uint32_t which, uint32_t value)
{
	uint8_t bits[2] = {0};
	put_bits(bits, 0, which, 4);
	put_bits(bits, 4, value, 8);
	return command(OP_SETTING, 0, bits, 2);
}

static void flush_commands(void)
{
	while (g.out_tail != g.out_head) {
		const uint32_t r = app_fpga_write(6, g.out[g.out_tail % OUT_SIZE]);
		if (r == 1) {
			return; /* the firmware's queue is full: next tick */
		}
		if (r != 0) {
			g.refused++;
		}
		g.out_tail++;
	}
}

/* ---- records ---------------------------------------------------------------------------- */

static int32_t s32(const uint16_t* w)
{ return (int32_t) ((uint32_t) w[0] | ((uint32_t) w[1] << 16)); }

static int32_t iabs(int32_t v)
{ return v < 0 ? -v : v; }

/* ---- navigation ------------------------------------------------------------------------ */

/* The almanac, from any channel's subframes 4 and 5 (gps_nav_almanac). */
static struct gps_alm alm[32] APP_STATE;

static void nav_subframe(struct chnav* v)
{
	uint32_t tow;
	gps_nav_subframe(v->words, &v->eph, &tow);
	gps_nav_almanac(v->words, alm);
	v->subframes++;
	v->tow = tow;
	v->anchor_record = v->sf_record;
	v->anchor_t = tow * 6.0 - 6.0; /* the HOW counts the next subframe's start */
	v->anchored = 1;
}

/* One data bit; first is the record that opened the bit's first period. */
static void nav_bit(struct chnav* v, uint32_t bit, uint32_t first)
{
	v->sr = (v->sr << 1) | bit;
	v->bits++;
	uint32_t data;
	if (v->nwords) {
		if (++v->wbits < 30) {
			return;
		}
		v->wbits = 0;
		if (!gps_nav_parity(
			    (uint32_t) (v->sr & 0x3FFFFFFF),
			    (int) ((v->sr >> 31) & 1),
			    (int) ((v->sr >> 30) & 1),
			    &data)) {
			v->nwords = 0;
			return;
		}
		v->words[v->nwords++] = data;
		if (v->nwords == 10) {
			nav_subframe(v);
			v->nwords = 0;
		}
		return;
	}
	/* A TLM with the preamble and a HOW with a subframe ID, each passing parity. Parity takes
	 * either polarity, and the preamble is tested after the D30* inversion is removed. */
	if (v->bits < 62) {
		return;
	}
	uint32_t tlm, how;
	if (!gps_nav_parity(
		    (uint32_t) ((v->sr >> 30) & 0x3FFFFFFF),
		    (int) ((v->sr >> 61) & 1),
		    (int) ((v->sr >> 60) & 1),
		    &tlm) ||
	    ((tlm >> 16) != 0x8B)) {
		return;
	}
	if (!gps_nav_parity(
		    (uint32_t) (v->sr & 0x3FFFFFFF),
		    (int) ((v->sr >> 31) & 1),
		    (int) ((v->sr >> 30) & 1),
		    &how)) {
		return;
	}
	const uint32_t id = (how >> 2) & 7;
	if ((id < 1) || (id > 5)) {
		return;
	}
	v->words[0] = tlm;
	v->words[1] = how;
	v->nwords = 2;
	v->wbits = 0;
	v->sf_record = first - 20 * 59;
}

/*
 * The tuner's offset from 1575.42 MHz as a carrier word (message 'F'): the synthesizer lands a few
 * hertz off (9.46 Hz low at this tuning), and that offset is in every carrier word: the frequency
 * measurement takes it out.
 */
static double nav_lo_word APP_STATE;

/*
 * PLL B's steering (message 'B'). PLL B clocks the tuner's LOs and the samples are on PLL A, so a
 * ratio change moves every carrier by -(the analog LO) times the change, and no code rate. The app
 * adds that to every carrier word it sends, outside the loops: the loops, and the code aid they
 * take from their own carrier, never see the step. A large step waits for a tick with no channel's
 * request outstanding, a small one goes at once (steer_tick); every channel's words follow, and the
 * firmware sets the ratio after the call's register writes.
 *
 * The frequency measurement takes from each record's carrier word the steering word its request
 * carried (known once a record shows the request's words): what remains is the carrier an
 * unsteered receiver would track, so the estimate is the TCXO's with no steering term. Between
 * the write and a channel's next period boundary (about a millisecond) the LO has moved and the
 * word not yet: a step over LARGE_STEP leaves out the window open across it and the next.
 */
#define LARGE_STEP 10737 /* carrier word for 10 Hz: 1 ms of it is 0.006 ppb of a second's window */
#define PLL_B_NOMINAL 32
static struct {
	uint8_t want;    /* a ratio waits for the tick */
	uint8_t refused; /* the last one was out of the firmware's bounds */
	uint8_t log_fw;  /* the latest log entry waits for the firmware's cycle counts */
	uint32_t a, b, c;
	double delta;     /* the ratio in effect over nominal, less one */
	int32_t word;     /* the carrier word added for it */
	double lo_hz;     /* the analog LO, on PLL B (message 'F') */
	uint32_t writes;
	struct {
		uint32_t a, b, c;
		uint32_t tick;       /* the cycle counter at the tick that asked */
		uint32_t n;          /* the highest active channel's last period start then */
		uint32_t start, end; /* the cycle counter around the firmware's write */
	} log[4]; /* the last four, the latest at (writes - 1) % 4 */
} steer APP_STATE;

/* The bit whose last period is end: its sum, from the periods received, to the lock and the
 * subframe reader. A record carries the start of the period it opens and the sums of the one
 * before, so a bit whose last sums came in record idx began in the period record idx - 20 opened. */
static void nav_bit_end(struct channel* ch, uint32_t end)
{
	struct chnav* v = &ch->nav;
	v->have_acc = 0;
	ch->lock_i += (double) v->acc * (double) v->acc;
	ch->lock_q += (double) v->accq * (double) v->accq;
	nav_bit(v, v->acc < 0, end - 20);
}

static void nav_record(
	struct channel* ch,
	uint32_t n0,
	uint32_t F,
	uint32_t W,
	const struct gps_dump* d)
{
	struct chnav* v = &ch->nav;
	const uint32_t idx = ch->loops.count - 1;
	v->last[1] = v->last[0];
	v->last[0].rec = idx;
	v->last[0].n0 = n0;
	v->last[0].F = F;
	v->last[0].W = W;
	if ((ch->loops.symbol_ms != 20) || (ch->loops.edge < 0)) {
		return;
	}
	const uint32_t pos = (idx - (uint32_t) ch->loops.edge) % 20;
	if ((pos == 0) || v->resume) {
		v->acc = 0;
		v->accq = 0;
		v->have_acc = 1;
		v->resume = 0;
	}
	if (!v->have_acc) {
		return;
	}
	v->acc += d->ip;
	v->accq += d->qp;
	if (pos == 19) {
		nav_bit_end(ch, idx);
	}
}

static int nav_usable(const struct channel* ch)
{
	const struct chnav* v = &ch->nav;
	return ch->active && (ch->prn < 120) && (ch->records >= 2) && v->anchored &&
		(v->last[0].rec - v->anchor_record < 30000) &&
		gps_nav_eph_complete(&v->eph) && (v->eph.health == 0);
}

/*
 * GPS time of week at which the code the channel received at sample n_r left the satellite; n_r
 * lies in the latest period or the one before.
 */
static double nav_transmit_time(const struct chnav* v, uint32_t n_r)
{
	const int k = (int32_t) (n_r - v->last[0].n0) < 0;
	const double t0 = v->anchor_t +
		(double) (int32_t) (v->last[k].rec - v->anchor_record) * 1e-3 +
		(double) v->last[k].F / TWO32 / CHIP_RATE;
	return t0 +
		(double) (int32_t) (n_r - v->last[k].n0) * ((double) v->last[k].W / TWO32) /
		CHIP_RATE;
}

/*
 * The steering loop (message 'D'; gps/freq.py's Servo is the reference). Each new estimate of the
 * TCXO's error joins the last SERVO_N; the line through them, at the middle of the coming second,
 * is the error to cancel, so PLL B's steering is its negative: the nearest ratio with c up to
 * 2^20 - 1, within SERVO_MAX_STEP of the last asked for (each step's carrier change stays under
 * LARGE_STEP, so no window is lost), and not asked for when it changes the steering less than
 * SERVO_MIN_CHANGE. CLKOUT's error is then the TCXO's less its prediction. No ratio lies within
 * 29.8 ppb of nominal (32 +- 1/1048575): a target in that gap becomes whichever of -29.8, 0 and
 * +29.8 is nearest the line's own (a large step, across the gap if need be), so a TCXO within
 * about 15 ppb of nominal is left uncorrected and one that crosses nominal is followed across.
 */
#define SERVO_N          10
#define SERVO_HORIZON    1.0   /* s past the latest estimate's middle */
#define SERVO_MAX_STEP   5.0   /* ppb a write: 7.9 Hz at the carriers */
#define SERVO_MIN_CHANGE 0.002 /* ppb */
#define SERVO_MAX        4000.0 /* ppb; the firmware's bound is 5 ppm */
#define SERVO_GAP        29.80232 /* ppb: 1e9 / (32 * 1048575), the nearest ratios to nominal */
static struct {
	uint8_t on, have;
	int n;
	/* The estimates' middles (s, from the latest window's end) and values (ppb), oldest first:
	 * single precision serves (0.002 ppb in 4,000 is 5e-7) and saves the double library's calls. */
	float t[SERVO_N], e[SERVO_N];
	uint32_t end_last; /* the latest window's end, in samples */
	double target;                 /* the steering last asked for, ppb */
} servo APP_STATE;

static void mgr_second(void);
static struct {
	uint8_t on;     /* from the first search */
	uint8_t quiet;  /* seconds of the M4's clock since the last call */
	uint16_t idle;  /* seconds since the search went idle with a channel free */
} mgr APP_STATE;

/* The frequency estimate's state (freq_tick). */
static struct {
	/* The held position: 0 none, 1 surveying (the mean so far), 2 surveyed, 3 given. */
	uint8_t pos_state;
	uint32_t survey_n, surveyed;
	double sum[3], rx[3];
	/* Before the first kept fix, three satellites with the height held (message 'H', above the
	 * ellipsoid) replace the held position each second when there is none, it is theirs (from3),
	 * or the host marked it rough ('H', a grid square's center); flip says the satellites
	 * disagreed at one of theirs, so the next such fix starts from its mirror. */
	float height;
	uint8_t from3, rough, flip;
	/*
	 * Windows run between grid points FIX_SAMPLES apart; est_end is the end of the window to
	 * estimate next.
	 */
	uint8_t have_grid;
	uint32_t est_end;
	uint8_t state; /* 0 waiting, 1 the expected frequencies, one end a tick */
	int n, k, end;
	uint8_t chs[CHANNELS];
	double t[CHANNELS], off[CHANNELS][2], clk[CHANNELS][2];
	/* The last window's estimate: see freq_reply. */
	uint32_t estimates, window_end;
	uint8_t valid, used, result;
	float est, rms, ppb[CHANNELS];
	int8_t lock[CHANNELS];
	/* The loops' Doppler less the expected carrier, a used satellite's at the last estimate: what
	 * an almanac's expected carrier is offset by in the loops (the TCXO, the tuner's offset and the
	 * quarter shift together; not PLL B's steering, which the samples show and the loops do not). */
	double bias;
	uint8_t have_bias;
} fr APP_STATE;


static struct {
	uint8_t state; /* 0 waiting, 1 solving, 2 the DOP */
	uint8_t have_next;
	uint32_t next_n;
	int n;
	uint8_t used;
	uint32_t n_r;
	double t_rx;
	double t_tx[CHANNELS];
	const struct gps_eph* eph[CHANNELS];
	struct gps_solve s;
	/* The last fix. */
	uint32_t fixes, fix_n;
	uint8_t fix_sats, fix_iterations, fix_used, tried_sats;
	uint8_t result; /* of the last attempt: see nav_reply */
	uint32_t rejected;
	double fix_t, fix_x[3];
	float fix_dop[5];
} nav APP_STATE;

/* What a tick's navigation did: 0 nothing, 1-4 a solve step (gps_nav.h), then these. */
enum {
	STEP_DOP = 6,   /* the DOP */
	STEP_START = 7, /* the transmit times of a new fix */
	STEP_FREQ_END = 8, /* a satellite position for an expected frequency */
	STEP_FREQ = 9,     /* the frequency estimate */
	STEP_ACQ_CODE = 10, /* a code's spectrum */
	STEP_ACQ_BINS = 11, /* a block's bins, mixed */
	STEP_ACQ_FFT = 12,  /* their FFT's stages */
	STEP_ACQ_CORR = 13, /* a correlation row (whole or half bins) */
};

/*
 * One step of the fix a tick (a satellite position, a row, the local vertical, a solution, the
 * DOP): on the M4 a position costs a few hundred microseconds of the call's 2 ms.
 */
static int nav_tick(void)
{
	if (nav.state == 1) {
		const int step =
			gps_nav_solve_step(&nav.s, nav.n, nav.t_tx, nav.eph, nav.t_rx);
		if (nav.s.done) {
			nav.state = 2;
		}
		return step;
	}
	if (nav.state == 2) {
		double dop[5] = {0, 0, 0, 0, 0};
		const int have_dop = gps_nav_dop(&nav.s, dop);
		nav.state = 0;
		/* Four satellites in a poor geometry can send the solve far off the Earth or leave it
		 * unconverged (a DOP above 1,000 has been seen); such a fix is counted, not kept. Three
		 * with the height held have GDOP 18-67 on good skies, the vertical held: their bound is
		 * looser, their error mostly the height's. */
		const double* x = nav.s.x;
		const double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
		nav.result = !nav.s.converged                          ? 2 :
			((r2 < 6.0e6 * 6.0e6) || (r2 > 6.6e6 * 6.6e6)) ? 3 :
			!have_dop || !(dop[0] <= (nav.s.held ? 4 : 1) * MAX_GDOP) ? 4 :
									 0;
		if (nav.result != 0) {
			nav.rejected += !nav.s.held; /* the fixes that would be kept */
			return STEP_DOP;
		}
		if (nav.s.held) {
			/* Not a kept fix: the held position until the first kept one. */
			for (int i = 0; i < 3; i++) {
				fr.rx[i] = x[i];
			}
			fr.pos_state = 1;
			fr.from3 = 1;
			nav.result = 5;
			return STEP_DOP;
		}
		nav.fixes++;
		nav.fix_n = nav.n_r;
		nav.fix_sats = (uint8_t) nav.n;
		nav.fix_iterations = (uint8_t) nav.s.iteration;
		nav.fix_used = nav.used;
		nav.fix_t = nav.t_rx - nav.s.x[3] / GPS_NAV_C;
		for (int i = 0; i < 3; i++) {
			nav.fix_x[i] = nav.s.x[i];
		}
		for (int i = 0; i < 5; i++) {
			nav.fix_dop[i] = (float) dop[i];
		}
		if ((fr.pos_state < 2) && fr.survey_n) {
			fr.surveyed++;
			for (int i = 0; i < 3; i++) {
				fr.sum[i] += x[i];
				fr.rx[i] = fr.sum[i] / fr.surveyed;
			}
			fr.pos_state = (fr.surveyed >= fr.survey_n) ? 2 : 1;
		}
		return STEP_DOP;
	}
	/* Once a second of samples, at the earliest of the usable channels' latest period starts. */
	int n = 0, first = -1;
	uint32_t n_r = 0;
	for (int c = 0; c < CHANNELS; c++) {
		if (!nav_usable(&g.ch[c])) {
			continue;
		}
		const uint32_t n0 = g.ch[c].nav.last[0].n0;
		if ((first < 0) || ((int32_t) (n0 - n_r) < 0)) {
			n_r = n0;
		}
		if (first < 0) {
			first = c;
		}
		n++;
	}
	if (first < 0) {
		return 0;
	}
	if (!nav.have_next) {
		nav.have_next = 1;
		nav.next_n = n_r;
	}
	if ((int32_t) (n_r - nav.next_n) < 0) {
		return 0;
	}
	nav.next_n += FIX_SAMPLES;
	if ((int32_t) (n_r - nav.next_n) >= 0) {
		nav.next_n = n_r + FIX_SAMPLES; /* after a pause */
	}
	nav.tried_sats = (uint8_t) n;
	/* Three: with no held position, one theirs, or a rough one (a grid square's center: within
	 * kilometers, where a fix with the height held is within about one per km of height error); a
	 * warm start's surveyed position stays until the first kept fix. */
	const int held = (n == 3) && fr.survey_n && !fr.surveyed && (!fr.pos_state || fr.from3 || fr.rough);
	if ((n < 4) && !held) {
		nav.result = 1;
		return 0;
	}
	nav.n = 0;
	nav.used = 0;
	double latest = 0;
	for (int c = 0; c < CHANNELS; c++) {
		if (!nav_usable(&g.ch[c])) {
			continue;
		}
		const double t = nav_transmit_time(&g.ch[c].nav, n_r);
		nav.t_tx[nav.n] = t;
		nav.eph[nav.n] = &g.ch[c].nav.eph;
		/* Later by the week's wrap: in its last seconds one channel's time can have rolled over
		 * to the next week's (a few seconds) while another's has not. */
		const double later = t - latest;
		if ((nav.n == 0) || ((later > 0) ? (later < 302400.0) : (later < -302400.0))) {
			latest = t;
		}
		nav.n++;
		nav.used |= 1u << c;
	}
	nav.n_r = n_r;
	nav.t_rx = latest + 0.075; /* any guess within tens of milliseconds */
	gps_nav_solve_init(&nav.s);
	if (held) {
		nav.s.held = 1;
		nav.s.height = fr.height;
		if (fr.pos_state) {
			/* From the held position, so the solve keeps to its side; after a disagreement,
			 * from its mirror through the mean direction of the last solve's first three
			 * satellites (gps_nav_solve_init leaves them). */
			float c[3], m = 0, cm = 0; /* a start: single precision is ample */
			for (int i = 0; i < 3; i++) {
				c[i] = (float) (nav.s.sv[0][i] + nav.s.sv[1][i] + nav.s.sv[2][i]);
				cm += c[i] * c[i];
				m += c[i] * (float) fr.rx[i];
			}
			for (int i = 0; i < 3; i++) {
				nav.s.x[i] = fr.flip ? (double) (2 * m / cm * c[i] - (float) fr.rx[i]) : fr.rx[i];
			}
			fr.flip = 0;
		}
	}
	nav.state = 1;
	return STEP_START;
}

/* ---- frequency ---------------------------------------------------------------------------- */

static void freq_close(struct channel* ch)
{
	const double tot = ch->lock_i + ch->lock_q;
	ch->done.end = ch->next_b;
	ch->done.lock = (tot > 0) ? (float) ((ch->lock_i - ch->lock_q) / tot) : 0.0f;
	ch->done.span = ch->fw.n_last - ch->fw.n_first;
	ch->done.valid = !ch->fw.bad && (ch->fw.records > 500) && (ch->done.lock >= MIN_LOCK) &&
		nav_usable(ch);
	if (ch->done.valid) {
		ch->done.f = gps_freq_measured(&ch->fw, nav_lo_word, FS);
		ch->done.t_end = nav_transmit_time(&ch->nav, ch->fw.n_last);
	}
}

/*
 * A GPS channel's record, after nav_record: the record that first starts at or past a grid point
 * closes the window begun at the point before and opens the next.
 */
static void freq_record(struct channel* ch, uint32_t n0, uint32_t wc)
{
	if (!fr.have_grid) {
		fr.have_grid = 1;
		fr.est_end = n0 + FIX_SAMPLES;
	}
	if (!ch->have_b) { /* the first grid point after n0 */
		uint32_t b = fr.est_end;
		while ((int32_t) (n0 - b) >= 0) {
			b += FIX_SAMPLES;
		}
		while ((int32_t) (n0 - (b - FIX_SAMPLES)) < 0) {
			b -= FIX_SAMPLES;
		}
		ch->next_b = b;
		ch->have_b = 1;
		ch->fw.open = 0;
		return;
	}
	if (ch->fw.open) {
		gps_freq_record(&ch->fw, n0, wc);
	}
	if ((int32_t) (n0 - ch->next_b) < 0) {
		return;
	}
	const uint32_t past = n0 - ch->next_b;
	if (ch->fw.open && (past < FIX_SAMPLES)) { /* open, it began at the point before next_b */
		freq_close(ch);
	}
	ch->next_b += (past / FIX_SAMPLES + 1) * FIX_SAMPLES;
	gps_freq_open(&ch->fw, n0, wc);
	ch->fw.bad = ch->fw_hold;
	ch->fw_hold = 0;
	ch->lock_i = ch->lock_q = 0;
}

static void servo_estimate(uint32_t end, double est)
{
	if (!servo.on || !((est < SERVO_MAX) && (est > -SERVO_MAX))) {
		return; /* a NaN or a wild estimate never reaches PLL B */
	}
	const int32_t gap = (int32_t) (end - servo.end_last);
	if ((gap <= 0) || (gap > SERVO_N * FS)) {
		/* After a pause longer than the history (every channel lost, say) a line through
		 * estimates minutes old does not predict the TCXO now, and past 537 s the gap wraps:
		 * start over. */
		servo.have = 0;
		servo.n = 0;
	}
	const float step = servo.have ? (float) gap * (float) (1.0 / FS) : 0.0f;
	servo.have = 1;
	servo.end_last = end;
	const int drop = servo.n == SERVO_N;
	for (int i = drop; i < servo.n; i++) {
		servo.t[i - drop] = servo.t[i] - step;
		servo.e[i - drop] = servo.e[i];
	}
	servo.n -= drop;
	servo.t[servo.n] = -0.5f;
	servo.e[servo.n] = (float) est;
	servo.n++;
	const double raw = -gps_freq_predict(servo.t, servo.e, servo.n, (float) (-0.5 + SERVO_HORIZON));
	double target = raw;
	if (target > servo.target + SERVO_MAX_STEP) {
		target = servo.target + SERVO_MAX_STEP;
	} else if (target < servo.target - SERVO_MAX_STEP) {
		target = servo.target - SERVO_MAX_STEP;
	}
	target = (target > SERVO_MAX) ? SERVO_MAX : (target < -SERVO_MAX) ? -SERVO_MAX : target;
	if ((target < SERVO_GAP) && (target > -SERVO_GAP)) {
		/* Clamped into the gap from either side, it never left it: the step limit kept the
		 * nearest ratio at the gap's edge. */
		target = (raw > SERVO_GAP / 2) ? SERVO_GAP : (raw < -SERVO_GAP / 2) ? -SERVO_GAP : 0;
	}
	uint32_t a = (target >= 0) ? PLL_B_NOMINAL : PLL_B_NOMINAL - 1, b, c;
	gps_freq_fraction(
		(target >= 0) ? PLL_B_NOMINAL * target * 1e-9 : 1.0 + PLL_B_NOMINAL * target * 1e-9,
		1048575,
		&b,
		&c);
	if (b == c) {
		a++;
		b = 0;
		c = 1;
	}
	const double got = ((double) ((int32_t) a - PLL_B_NOMINAL) * c + b) / (PLL_B_NOMINAL * (double) c) * 1e9;
	const double d = got - servo.target;
	if ((d < SERVO_MIN_CHANGE) && (d > -SERVO_MIN_CHANGE)) {
		return;
	}
	steer.a = a; /* replaces one not yet written */
	steer.b = b;
	steer.c = c;
	steer.want = 1;
	servo.target = got;
}

static void freq_finish(void)
{
	/*
	 * The windows hold the carriers less their steering words (freq_record), so measured less
	 * expected is -(the digital center) times the TCXO's error. PLL B's own error, CLKOUT's, is
	 * the TCXO's plus the steering.
	 */
	double ppb[CHANNELS], bias[CHANNELS], var[CHANNELS], est = 0;
	for (int i = 0; i < fr.n; i++) {
		const struct channel* ch = &g.ch[fr.chs[i]];
		const double f = gps_freq_expected(fr.off[i], fr.clk[i]);
		ppb[i] = -(ch->done.f - f) / GPS_FREQ_L1 * 1e9;
		bias[i] = gps_loops_doppler(&ch->loops) - f;
		var[i] = ch->fvar;
	}
	const uint32_t used = gps_freq_combine(ppb, var, fr.n, REJECT_PPB, &est);
	if (!used || ((fr.n > 1) && !(used & (used - 1)))) {
		/* Fewer than two agree (two satellites apart, or three scattered): one is wrong (a false
		 * lock, a wrong ephemeris, a held position far off) and nothing says which. No estimate,
		 * so none reaches the steering. */
		fr.result = 4;
		fr.flip = fr.from3 && !fr.surveyed;
		return;
	}
	double ss = 0;
	int k = 0;
	fr.used = 0;
	for (int i = 0; i < fr.n; i++) {
		fr.ppb[fr.chs[i]] = (float) ppb[i];
		gps_freq_var_update(&g.ch[fr.chs[i]].fvar, ppb[i] - est, REJECT_PPB);
		if (used & (1u << i)) {
			fr.bias = bias[i];
			fr.have_bias = 2; /* this run's (1: a warm start's) */
			fr.used |= 1u << fr.chs[i];
			ss += (ppb[i] - est) * (ppb[i] - est);
			k++;
		}
	}
	fr.est = (float) est;
	fr.rms = (float) gps_nav_sqrt(ss / k);
	fr.estimates++;
	fr.result = 0;
	servo_estimate(fr.window_end, est);
}

/*
 * Once every armed GPS channel's records pass a window's end (or one is 100 ms past it, another
 * stalled): each satellite's expected frequency, one satellite position a tick, then the estimate.
 */
static int freq_tick(void)
{
	if (fr.state == 1) {
		const int i = fr.k;
		gps_freq_end(
			&g.ch[fr.chs[i]].nav.eph,
			fr.t[i],
			fr.rx,
			fr.end,
			&fr.off[i][fr.end],
			&fr.clk[i][fr.end]);
		fr.end ^= 1;
		if (fr.end || (++fr.k < fr.n)) {
			return STEP_FREQ_END;
		}
		freq_finish();
		fr.state = 0;
		return STEP_FREQ;
	}
	if (!fr.have_grid) {
		return 0;
	}
	int armed = 0, waiting = 0, late = 0;
	uint32_t latest = fr.est_end;
	for (int c = 0; c < CHANNELS; c++) {
		const struct channel* ch = &g.ch[c];
		if (!ch->active || (ch->prn >= 120) || !ch->have_b) {
			continue;
		}
		const int32_t d = (int32_t) (ch->last_n - fr.est_end);
		armed = 1;
		waiting |= d < 0;
		late |= d >= (int32_t) (FIX_SAMPLES / 10);
		if ((int32_t) (ch->last_n - latest) > 0) {
			latest = ch->last_n;
		}
	}
	if (!armed || (waiting && !late)) {
		return 0;
	}
	const uint32_t end = fr.est_end;
	fr.est_end += FIX_SAMPLES;
	while ((int32_t) (latest - fr.est_end) >= (int32_t) FIX_SAMPLES) {
		fr.est_end += FIX_SAMPLES; /* after a pause */
	}
	fr.window_end = end;
	fr.n = 0;
	fr.valid = 0;
	for (int c = 0; c < CHANNELS; c++) {
		const struct channel* ch = &g.ch[c];
		fr.lock[c] = 0;
		if (!ch->active || (ch->prn >= 120) || (ch->done.end != end)) {
			continue;
		}
		const float l = ch->done.lock * 100.0f;
		fr.lock[c] = (int8_t) (l > 100.0f ? 100 : l < -100.0f ? -100 : l);
		if (!ch->done.valid) {
			continue;
		}
		fr.valid |= 1u << c;
		fr.chs[fr.n] = (uint8_t) c;
		/* The transmit time at the window's middle; within microseconds is enough. */
		fr.t[fr.n] = ch->done.t_end - ch->done.span / (2.0 * FS);
		fr.n++;
	}
	mgr_second();
	if (fr.pos_state == 0) {
		fr.result = 1;
		return 0;
	}
	if (fr.n == 0) {
		fr.result = 2;
		return 0;
	}
	fr.k = 0;
	fr.end = 0;
	fr.state = 1;
	return 0;
}

static void handle_record(void)
{
	const uint32_t c = g.rec[0] >> 8;
	struct channel* ch = &g.ch[c];
	g.records++;
	if (!ch->active) {
		return;
	}
	const uint32_t n0 = (uint32_t) s32(&g.rec[1]);
	/* A record misread (words that look like a record start) carries a random period start: one
	 * more than 6 periods from the last, or behind it, is left out, unless three in a row (a real
	 * resynchronization). */
	if ((ch->records != 0) && (n0 - ch->last_n - 2000 > 22000) && (++ch->odd < 3)) {
		ch->gaps++;
		return;
	}
	ch->odd = 0;
	if ((ch->records != 0) && (n0 - ch->last_n > 6000)) {
		/* A period is about 4,000 samples. Records missed (a snapshot pauses the stream) are
		 * counted by their periods, so the record index, which bit timing and the transmit
		 * time go by, stays on time. A gap shorter than a bit drops no bit (a lost bit shifts the
		 * subframe's words and fails their parity to the next preamble): the missed periods add
		 * nothing to the bit in progress, and a gap across its end closes it from the periods
		 * received, the next record starting the bit after. */
		ch->gaps++;
		const uint32_t missed = (n0 - ch->last_n + 2000) / 4000 - 1;
		const uint32_t last = ch->loops.count - 1;
		const uint32_t pos = (last - (uint32_t) ch->loops.edge) % 20;
		if ((ch->loops.symbol_ms == 20) && (ch->loops.edge >= 0) && (pos + missed >= 19)) {
			if (ch->nav.have_acc) {
				nav_bit_end(ch, last + 19 - pos);
			}
			ch->nav.resume = 1;
		}
		gps_loops_skip(&ch->loops, missed);
	}
	ch->last_n = n0;
	ch->records++;
	const struct gps_dump d = {
		s32(&g.rec[9]),
		s32(&g.rec[11]),
		s32(&g.rec[13]),
		s32(&g.rec[15]),
		s32(&g.rec[17]),
		s32(&g.rec[19]),
	};
	ch->prompt = iabs((int32_t) d.ip) + iabs((int32_t) d.qp);
	/*
	 * A record carries the words in use from its period's start. The engine holds a pending
	 * request until the channel's previous one is applied, and the command port drops whatever
	 * command finishes meanwhile (no flow control reaches the M4); so a channel has one request
	 * outstanding at a time, released when a record shows its words. A request sent after record
	 * k applies at the next boundary, or the one after when the record came late.
	 */
	if (ch->outstanding) {
		if (((uint32_t) s32(&g.rec[5]) == ch->sent_W) &&
		    ((uint32_t) s32(&g.rec[7]) == ch->sent_Wc)) {
			ch->outstanding = 0;
			ch->steer_in_use = ch->sent_steer;
		} else if (++ch->waited >= 3) {
			ch->outstanding = 0;
			ch->lost++;
			ch->send = 1;
		}
	}
	int64_t code, carrier;
	gps_loops_update(&ch->loops, &d, &code, &carrier);
	nav_record(
		ch,
		n0,
		(uint32_t) s32(&g.rec[3]),
		(uint32_t) s32(&g.rec[5]),
		&d);
	if (ch->prn < 120) {
		freq_record(ch, n0, (uint32_t) s32(&g.rec[7]) - (uint32_t) ch->steer_in_use);
	}
	if ((code != ch->code) || (carrier != ch->carrier)) {
		ch->code = code;
		ch->carrier = carrier;
		ch->send = 1;
	}
}

/* ---- acquisition (gps/app/gps_acq.c) ------------------------------------------------------ */

static struct {
	float spec[2 * GPS_ACQ_N]; /* a block's spectrum */
	float work[2 * GPS_ACQ_N];
	int16_t code[2 * GPS_ACQ_N]; /* a code's conjugate spectrum */
	float acc[GPS_ACQ_N], acc1[GPS_ACQ_N];
	/* The snapshot being read: where its samples begin in the ring, its trailer. */
	uint32_t snap_byte;
	uint16_t trailer[3];
	/* Completed snapshots not yet binned: their place in the ring and first sample's index, the
	 * oldest first. A request goes out while one is in flight (with nothing else streaming, a
	 * snapshot's last words reach the ring only when more data follows), so two can complete
	 * before the search takes one; each stays in the ring for the 3 snapshots after it. */
	uint8_t q_n;
	uint32_t q_byte[2], q_first[2];
	uint32_t newest;
	uint32_t ticks_work, ticks_wait; /* the search's ticks: an operation, or none for want of a snapshot */
	/* The search: mode, then the PRN, Doppler and block in hand. */
	uint8_t mode;
	uint8_t blocks, refine_blocks;
	uint8_t prn_i, code_prn, phase, blk, op, code_op;
	uint8_t pair; /* cold: the blocks serve dop_i + 4 (1 kHz on) too, its sums in acc1 */
	int8_t dop_i;
	uint32_t first0, blk_first;
	float doppler;
	struct {
		float metric;
		int16_t dop_i;
	} best[32]; /* the cold pass's best cell per GPS PRN */
	/* Refining: the candidates by metric, the one in hand and its best so far. */
	uint8_t cand[CHANNELS], ncand, cand_i;
	float ref_doppler;
	int ref_cell;
	uint32_t ref_first0;
	uint8_t armed[CHANNELS], narmed; /* PRNs armed by the search */
	uint8_t pass_armed;              /* armed in this pass */
	/* Each GPS PRN's Doppler indices for the cold pass: all, the three around its almanac's
	 * prediction, or none (lo > hi: below the horizon or unhealthy). */
	int8_t lo[32], hi[32];
	/* The prediction: GPS time, then one satellite end a tick. */
	double pred_t;
	uint8_t pred_end;
	uint8_t predicted; /* this pass's PRNs and Dopplers from the almanac */
	/* The host's GPS time (message 'W'), carried on by the M4's clock: a prediction with no
	 * usable channel. */
	double warm_t;
	uint8_t warm;
	/* Each PRN's predicted Doppler less PLL B's steering (Hz), 0 where none: an arm's Doppler. */
	float pred_hz[32];
	/* Locating (acq_locate): the tracked satellites (PRN index, loops' Doppler), the PRNs with a
	 * position and velocity (in work, six floats each), the first of the two positions, the
	 * grid's ring and point (its latitude's and longitude's cos and sin, the longitude's step),
	 * and whether any point was kept. */
	uint8_t loc_n, loc_k[CHANNELS], loc_ring, loc_kept;
	uint16_t loc_i, loc_np;
	float loc_m[CHANNELS];
	float loc_cl, loc_sl, loc_co, loc_so, loc_cs, loc_ss;
	double loc_s0[3];
	uint32_t loc_ok;
} acq APP_STATE;

static uint32_t arm_channel(
	uint32_t c,
	uint32_t prn,
	uint32_t W,
	uint32_t Wc,
	uint32_t g1,
	uint32_t g2,
	uint32_t win,
	uint32_t start,
	int32_t doppler_mhz);

/*
 * The search on the box (message 'C'). Cold: every GPS PRN at every Doppler from -ACQ_DOPPLER to
 * +ACQ_DOPPLER in ACQ_STEP steps, each (PRN, Doppler) over its own blocks (snapshots of 4,096
 * samples, binned as they complete: no room to keep blocks for every PRN), whole bins, the best
 * cell kept per PRN; Dopplers 1 kHz apart in pairs, one block's bins and FFT serving both. Each
 * block starts where its bins fall on the first block's (gps_acq_skip): a block half a bin off
 * would count only in the half-bin cells, which the cold pass does not keep. The candidates (over
 * the detection threshold and not a strong one's cross-correlation) are checked in turn: once at
 * the cold pass's Doppler, whole and half bins, fresh snapshots; one still over the threshold is
 * armed there on a free channel, its code epoch carried from the checking snapshots to ACQ_AHEAD
 * samples past the newest seen. A pass that arms nothing, or leaves fewer than four channels
 * tracking, starts over with twice the blocks, up to ACQ_MAX_BLOCKS, the PRNs already tracked left
 * out. With an almanac a pass predicts first (acq_locate): a PRN below the horizon is left out, a
 * visible one searched at three Dopplers instead of ACQ_NDOP.
 */
#define ACQ_MAX_BLOCKS 64
#define ACQ_DOPPLER 4500.0
#define ACQ_STEP    250.0
#define ACQ_NDOP    37 /* 2 * ACQ_DOPPLER / ACQ_STEP + 1 */
#define ACQ_AHEAD   200000u /* 50 ms */
#define ACQ_SNAP    3 /* setting 1: (3 + 1) * 1,024 samples */
enum { ACQ_IDLE = 0, ACQ_COLD = 1, ACQ_REFINE = 2, ACQ_PREDICT = 3, ACQ_LOCATE = 4 };
enum { PH_WAIT = 1, PH_BLOCK = 2, PH_PEAK = 3 };

/* The first sample of a snapshot begins at ring byte at (parse_word). */
static void acq_snap_start(uint32_t at) { acq.snap_byte = at % APP_SAMPLES_SIZE; }

/* A snapshot's last three words have passed: its first sample's index, and its bins when the
 * search wants a block. */
static void acq_snap_done(void)
{
	if ((acq.trailer[0] & 0xFF) != 0xD7) {
		return;
	}
	const uint32_t count = ((acq.trailer[0] >> 8) + 1u) * 1024;
	const uint32_t lo = acq.trailer[1], hi = acq.trailer[2];
	const uint32_t second = ((((hi - (lo == 0xFFFF)) & 0xFFFF)) << 16) | lo;
	const uint32_t first = second - 2 - count + 1;
	acq.newest = first + count;
	if ((acq.mode == ACQ_IDLE) || (acq.q_n == 2) || (count != (ACQ_SNAP + 1u) * 1024)) {
		/* only the search's own: a host's ('S') can be long, its start long overwritten */
		return;
	}
	acq.q_byte[acq.q_n] = acq.snap_byte;
	acq.q_first[acq.q_n] = first;
	acq.q_n++;
}

static void acq_request_snapshot(void)
{
	if (!g.snap_pending && (acq.q_n < 2) && send_setting(1, ACQ_SNAP)) {
		g.snap_pending = 1;
		g.snap_words = (ACQ_SNAP + 1u) * 1024 + 3;
	}
}

/*
 * Thresholds, k / sqrt(blocks) over 1. On a 5-minute recording the best of 37 Dopplers x 1,024
 * cells for an absent PRN reached 1 + 8.4 and 7.8 / sqrt(blocks) at 8 and 16 blocks (200 and 120
 * trials; 99% under 7.8 and 7.3): the cold pass passes candidates from 8, a noise one now and then
 * (its check costs about 0.3 s). The check of
 * a candidate, one Doppler and fresh snapshots, arms over 8: on a 5-minute recording an absent
 * PRN's check reached 1 + 6.3 / sqrt(blocks) at most in 96 (32 and 64 blocks), and two satellites
 * near the cold pass's limit (about 29-31 dB-Hz) passed 92% and 98% of theirs.
 */
static float acq_threshold(int blocks, float k) { return 1.0f + k / __builtin_sqrtf((float) blocks); }

/* After the cold pass: the candidates, strongest first, over the threshold and not within 2% of
 * the strongest's excess (a cross-correlation of it). */
static void acq_candidates(void)
{
	float top = 0;
	for (int i = 0; i < 32; i++) {
		top = (acq.best[i].metric > top) ? acq.best[i].metric : top;
	}
	acq.ncand = 0;
	const float thr = acq_threshold(acq.blocks, 8.0f);
	/* No more than the free channels: each candidate refined costs seconds, and with no channel
	 * left one would be refined for nothing. */
	int free = 0;
	for (int c = 0; c < CHANNELS; c++) {
		free += !g.ch[c].active;
	}
	for (int n = 0; n < free; n++) {
		int pick = -1;
		for (int i = 0; i < 32; i++) {
			int taken = 0;
			for (int k = 0; k < acq.ncand; k++) {
				taken |= acq.cand[k] == i;
			}
			const float m = acq.best[i].metric;
			if (!taken && (m >= thr) &&
			    ((m == top) || (m - 1 >= 0.02f * (top - 1))) &&
			    ((pick < 0) || (m > acq.best[pick].metric))) {
				pick = i;
			}
		}
		if (pick < 0) {
			break;
		}
		acq.cand[acq.ncand++] = (uint8_t) pick;
	}
}

/* Arm channel c on prn: its code's arm values generated, the carrier word from Doppler fd (Hz, the
 * loops' own: PLL B's steering not in it, arm_channel adds it to the words as to every word), chip 0
 * on sample start1. 0 when the engine took it. */
static uint32_t arm_at(uint32_t c, uint32_t prn, uint32_t W, double fd, uint32_t start1)
{
	uint8_t code[1023];
	uint32_t g1, g2, win;
	gps_acq_code((int) prn, code, &g1, &g2, &win);
	const double wc = fd / FS * TWO32;
	return arm_channel(
		c,
		prn,
		W,
		(uint32_t) (int32_t) (wc < 0 ? wc - 0.5 : wc + 0.5),
		g1,
		g2,
		win,
		start1 - 1,
		(int32_t) (fd * 1000.0));
}

/* The refined candidate on a free channel: the code epoch carried to ACQ_AHEAD past the newest. */
static void acq_arm(int prn)
{
	int c = 0;
	while ((c < CHANNELS) && g.ch[c].active) {
		c++;
	}
	if (c == CHANNELS) {
		return;
	}
	/* With this run's offset (fr.have_bias 2) the almanac's Doppler is within a hertz or two; the
	 * cold pass's bin is good to 125 Hz at best, as over 1 ms blocks the peak is nearly flat (25
	 * Hz steps across it once picked noise, up to 220 Hz off within the bin, and a weak satellite
	 * often peaks in the next bin). A warm start's offset is older. */
	/* In the loops' frame: the samples show PLL B's steering in the carrier and not in the code,
	 * whose rate is the signal's, so neither the rate nor the loops take it (at 105 ppb of
	 * steering a rate with it put the arm 0.2-0.8 chip off by the check's end). */
	double fd = acq.ref_doppler - (float) steer.word * (float) (FS / TWO32);
	if ((fr.have_bias == 2) && (acq.pred_hz[prn - 1] != 0.0f)) {
		fd = acq.pred_hz[prn - 1];
	}
	const double rate = GPS_ACQ_CHIP * (1.0 + fd / GPS_ACQ_L1);
	const double period = 1023.0 / rate * GPS_ACQ_FS;
	/* chip 0 at half bin ref_cell of the first block: ref_cell * BLOCK / (2 N) samples on */
	const double epoch = acq.ref_cell * ((double) GPS_ACQ_BLOCK / (2 * GPS_ACQ_N));
	const double ahead = (double) (int32_t) (acq.newest + ACQ_AHEAD - acq.ref_first0) - epoch;
	const double m = (double) (uint32_t) (ahead / period) + 1.0;
	const double off = epoch + m * period;
	const uint32_t start1 = acq.ref_first0 + (uint32_t) (off + 0.5);
	if (arm_at((uint32_t) c, (uint32_t) prn, (uint32_t) (rate / GPS_ACQ_FS * TWO32 + 0.5), fd, start1) ==
	    0) {
		acq.armed[acq.narmed++ % CHANNELS] = (uint8_t) prn;
		acq.pass_armed++;
	}
}

/* The Doppler of the (PRN, Doppler) in hand. */
static float acq_doppler_of(void)
{
	if (acq.mode == ACQ_COLD) {
		return -ACQ_DOPPLER + acq.dop_i * ACQ_STEP;
	}
	return -ACQ_DOPPLER + acq.best[acq.cand[acq.cand_i]].dop_i * ACQ_STEP;
}

static void acq_pass_done(void);
static void acq_cold_pass(void);
static void acq_begin(void);
static int mgr_locked(int c);

/* Refining over at least the blocks that found the candidate: a satellite found over 64 would
 * not hold its metric over 32. */
static int acq_refine_blocks(void)
{
	return (acq.blocks > acq.refine_blocks) ? acq.blocks : acq.refine_blocks;
}

/* An FFT's stages two a call (about 0.2 ms on the M4): the call's other work is then paid half as
 * often. A row is its product and the inverse stages, the accumulation with the last. */
#define PAIRS    ((GPS_ACQ_STAGES + 1) / 2)
#define CODE_OPS (2 + PAIRS)
#define ROW_OPS  (1 + PAIRS)
#define BLOCK_OPS(rows) (PAIRS + (rows) * ROW_OPS)

static void acq_stages(float* x, int inverse, int pair)
{
	for (int s = 2 * pair; (s < 2 * pair + 2) && (s < GPS_ACQ_STAGES); s++) {
		gps_acq_fft_stage(x, inverse, s);
	}
}

/* The cold pass from PRN index prn_i on: past the PRNs tracked or with no Doppler to search, to
 * its first Doppler; after the last, the candidates refined (or the pass done). */
static void acq_next(void)
{
	for (; acq.prn_i < 32; acq.prn_i++) {
		int tracked = acq.lo[acq.prn_i] > acq.hi[acq.prn_i];
		for (int c = 0; c < CHANNELS; c++) {
			tracked |= g.ch[c].active && (g.ch[c].prn == acq.prn_i + 1);
		}
		if (!tracked) {
			acq.dop_i = acq.lo[acq.prn_i];
			return;
		}
	}
	acq_candidates();
	acq.mode = ACQ_REFINE;
	acq.cand_i = 0;
	acq.dop_i = 0;
	if (!acq.ncand) {
		acq_pass_done();
	}
}

/*
 * The almanac's prediction (ACQ_PREDICT): each PRN at its Doppler seen from the held position,
 * offset as the tracked satellites are (fr.bias), a bin either side, and not at all below about
 * 5 degrees there. With no position (ACQ_LOCATE), where the tracked satellites' Dopplers agree: a
 * grid of points about 300 km apart on a sphere of the Earth's mean radius (the ellipsoid moves a
 * Doppler under 20 Hz), LOC_POINTS a tick. At a point that sees every tracked satellite, their
 * loops' Dopplers less the point's share a mean, the offset; a point within LOC_TOL of it on each
 * is kept. Each other PRN is searched at the Dopplers it shows from the kept points, a bin either
 * side, and not at all where none sees it. Either way a PRN without an almanac entry is searched at
 * every Doppler. On the almanac of a run, 12 skies over a day, every satellite visible within its
 * range: one tracked left about 650 of the 1,184 cells, two 250, three 90 (100 Hz missed some at
 * three).
 */
#define LOC_STEP   0.0470880f /* 300 km over 6,371 km, radians */
#define LOC_RINGS  67         /* pi / LOC_STEP */
#define LOC_MIN_NP 13         /* points a ring at least: a step under 0.5 rad (gps_acq_cis) */
#define LOC_TOL    200.0f
#define LOC_POINTS 32
#define LOC_R      6371e3f
#define LOC_LAMBDA 0.190293673f /* L1's wavelength, m */

/* PRN index p's Doppler (Hz) seen from point x, else 0 when below the horizon there (from a held
 * position, about 3-5 degrees: a range over 25,300 km). */
static int loc_doppler(int p, const float* x, float* f)
{
	const float* sv = &acq.work[6 * p];
	const float d0 = sv[0] - x[0], d1 = sv[1] - x[1], d2 = sv[2] - x[2];
	const float r = __builtin_sqrtf(d0 * d0 + d1 * d1 + d2 * d2);
	if ((d0 * x[0] + d1 * x[1] + d2 * x[2] <= 0) || (!acq.loc_n && (r > 25.3e6f))) {
		return 0;
	}
	*f = -(d0 * sv[3] + d1 * sv[4] + d2 * sv[5]) / (r * LOC_LAMBDA);
	return 1;
}

/* One point, the held position or the grid's: kept, each PRN's range widened to its bin there. */
static void loc_point(const float* x)
{
	float f[CHANNELS], b = (float) fr.bias;
	if (acq.loc_n) {
		b = 0;
		for (int j = 0; j < acq.loc_n; j++) {
			if (!loc_doppler(acq.loc_k[j], x, &f[j])) {
				return;
			}
			b += acq.loc_m[j] - f[j];
		}
		b /= (float) acq.loc_n;
		for (int j = 0; j < acq.loc_n; j++) {
			const float r = acq.loc_m[j] - f[j] - b;
			if ((r > LOC_TOL) || (r < -LOC_TOL)) {
				return;
			}
		}
	}
	acq.loc_kept = 1;
	/* The samples show PLL B's steering; the loops, and so the offset, do not. */
	const float bs = b + (float) steer.word * (float) (FS / TWO32) + (float) ACQ_DOPPLER;
	for (int p = 0; p < 32; p++) {
		float fp;
		if ((acq.loc_ok >> p & 1) && loc_doppler(p, x, &fp)) {
			if (!acq.loc_n) {
				acq.pred_hz[p] = fp + b;
			}
			/* In range before the conversion: a host's offset or almanac entry can be anything.
			 * A bin either side, -1 and ACQ_NDOP included: valid Dopplers too. */
			const float k = (fp + bs) * (1.0f / (float) ACQ_STEP) + 0.5f;
			if ((k >= 0) && (k < ACQ_NDOP)) {
				const int8_t lo = (int8_t) ((int) k - 1), hi = (int8_t) ((int) k + 1);
				acq.lo[p] = (lo < acq.lo[p]) ? lo : acq.lo[p];
				acq.hi[p] = (hi > acq.hi[p]) ? hi : acq.hi[p];
			}
		}
	}
}

/* The pass from the ranges: every Doppler for a PRN without an almanac entry, and for all when no
 * point was kept (a tracked Doppler wrong: a false lock). Deeper only once the almanac is whole: a
 * PRN without one keeps every Doppler, and 128 blocks over all of them took 300 s on the box. */
static void loc_done(void)
{
	for (int p = 0; p < 32; p++) {
		if (!acq.loc_kept || !alm[p].have) {
			acq.lo[p] = 0;
			acq.hi[p] = ACQ_NDOP - 1;
		}
		acq.predicted = acq.predicted && alm[p].have && acq.loc_kept;
	}
	if (!acq.predicted && (acq.blocks > ACQ_MAX_BLOCKS)) {
		acq.blocks = ACQ_MAX_BLOCKS; /* deeper only predicted (acq_pass_done) */
	}
	acq_cold_pass();
}

/* A tick of prediction or locating: a satellite's position and velocity from its almanac entry,
 * half a second either side of the time, a step each; then the held position, or LOC_POINTS of the
 * grid; then the pass. */
static int acq_locate(void)
{
	if (acq.prn_i < 32) {
		const int p = acq.prn_i;
		if (alm[p].have && !alm[p].health) {
			struct gps_eph e;
			double x[3], clk;
			gps_nav_alm_eph(&alm[p], &e);
			gps_nav_sat_position(&e, acq.pred_t + acq.pred_end - 0.5, x, &clk);
			float* sv = &acq.work[6 * p];
			for (int k = 0; k < 3; k++) {
				if (!acq.pred_end) {
					acq.loc_s0[k] = x[k];
				}
				sv[k] = (float) (0.5 * (x[k] + acq.loc_s0[k]));
				sv[3 + k] = (float) (x[k] - acq.loc_s0[k]);
			}
			if ((acq.pred_end ^= 1)) {
				return STEP_FREQ_END;
			}
			acq.loc_ok |= 1u << p;
		}
		acq.prn_i++;
		return STEP_FREQ_END;
	}
	if (!acq.loc_n) {
		const float x[3] = {(float) fr.rx[0], (float) fr.rx[1], (float) fr.rx[2]};
		loc_point(x);
		loc_done();
		return STEP_FREQ_END;
	}
	for (int n = 0; n < LOC_POINTS; n++) {
		if (acq.loc_i == acq.loc_np) {
			if (acq.loc_ring == LOC_RINGS) {
				loc_done();
				return STEP_FREQ_END;
			}
			/* The next ring's latitude, from the south pole half a step at a time. */
			float c, s;
			gps_acq_cis(acq.loc_ring++ ? LOC_STEP : 0.5f * LOC_STEP, &c, &s);
			const float cl = acq.loc_cl * c - acq.loc_sl * s;
			acq.loc_sl = acq.loc_cl * s + acq.loc_sl * c;
			acq.loc_cl = cl;
			const int np = (int) (2 * 3.14159265f * cl / LOC_STEP + 0.5f);
			acq.loc_np = (uint16_t) (np < LOC_MIN_NP ? LOC_MIN_NP : np);
			gps_acq_cis(2 * 3.14159265f / acq.loc_np, &acq.loc_cs, &acq.loc_ss);
			acq.loc_co = 1;
			acq.loc_so = 0;
			acq.loc_i = 0;
		}
		const float x[3] = {LOC_R * acq.loc_cl * acq.loc_co, LOC_R * acq.loc_cl * acq.loc_so,
				    LOC_R * acq.loc_sl};
		loc_point(x);
		const float co = acq.loc_co * acq.loc_cs - acq.loc_so * acq.loc_ss;
		acq.loc_so = acq.loc_co * acq.loc_ss + acq.loc_so * acq.loc_cs;
		acq.loc_co = co;
		acq.loc_i++;
	}
	return STEP_FREQ_END;
}

/*
 * One short operation of the search a tick, none much over 0.45 ms on the M4 (records share the
 * call's 2 ms once channels track): a code's spectrum in CODE_OPS parts (its chips, the FFT's
 * stages, the rounding); then each block as its bins mixed, the FFT's stages, and per row (whole
 * bins, and half bins when refining) the product, the inverse FFT's stages and the accumulation;
 * then the peak.
 */
static int acq_search(void)
{
	if (acq.mode >= ACQ_PREDICT) {
		return acq_locate();
	}
	const int prn = (acq.mode == ACQ_COLD) ? acq.prn_i + 1 : acq.cand[acq.cand_i] + 1;
	if (acq.code_prn != prn) {
		if (acq.code_op == 0) {
			uint8_t code[1023];
			uint32_t g1, g2, win;
			gps_acq_code(prn, code, &g1, &g2, &win);
			gps_acq_code_chips(code, acq.work);
		} else if (acq.code_op < CODE_OPS - 1) {
			acq_stages(acq.work, 0, acq.code_op - 1);
		} else {
			gps_acq_code_round(acq.work, acq.code);
			acq.code_prn = (uint8_t) prn;
		}
		acq.code_op = (acq.code_op + 1) % CODE_OPS;
		return STEP_ACQ_CODE;
	}
	const int cold = acq.mode == ACQ_COLD;
	const int rows = cold ? 1 + acq.pair : 2;
	if (acq.phase == PH_WAIT) {
		acq_request_snapshot();
		if (acq.q_n == 0) {
			return 0;
		}
		/* The oldest completed snapshot binned and mixed into spec, from the sample that
		 * puts its bins on the first block's (at the pair's middle Doppler). */
		acq.doppler = acq_doppler_of();
		uint32_t first = acq.q_first[0];
		if (acq.blk == 0) {
			acq.first0 = first;
			acq.pair = cold && !((acq.dop_i - acq.lo[acq.prn_i]) & 4) &&
				   (acq.dop_i + 4 <= acq.hi[acq.prn_i]);
			for (int j = 0; j < GPS_ACQ_N; j++) {
				acq.acc[j] = 0;
				acq.acc1[j] = 0;
			}
		}
		const uint32_t skip =
			gps_acq_skip(first - acq.first0, acq.doppler + (float) (500 * acq.pair));
		gps_acq_block((const void*) APP_SAMPLES, APP_SAMPLES_SIZE, acq.q_byte[0] + 2 * skip,
			      acq.doppler, acq.spec);
		first += skip;
		acq.q_byte[0] = acq.q_byte[1];
		acq.q_first[0] = acq.q_first[1];
		acq.q_n--;
		acq.blk_first = first;
		acq.phase = PH_BLOCK;
		acq.op = 0;
		return STEP_ACQ_BINS;
	}
	if (acq.phase == PH_BLOCK) {
		const int op = acq.op;
		int kind = STEP_ACQ_CORR;
		if (op < PAIRS) {
			acq_stages(acq.spec, 0, op);
			kind = STEP_ACQ_FFT;
		} else {
			const int r = (op - PAIRS) / ROW_OPS;
			const int k = (op - PAIRS) % ROW_OPS;
			/* Cold, row 1 is the pair's Doppler, 1 kHz on: one bin's rotation of the
			 * spectrum, into acc1 as whole bins. */
			if (k == 0) {
				gps_acq_product(acq.spec, acq.code, acq.work, !cold && r, cold && r);
			} else {
				acq_stages(acq.work, 1, k - 1);
			}
			if (k == PAIRS) {
				gps_acq_accumulate(
					acq.work,
					(cold && r) ? acq.acc1 : acq.acc,
					cold ? 0 : acq.acc1,
					gps_acq_shift(acq.blk_first - acq.first0,
						      acq.doppler + (float) (1000 * r)),
					!cold && r);
			}
		}
		if (++acq.op == BLOCK_OPS(rows)) {
			const int blocks = (acq.mode == ACQ_COLD) ? acq.blocks : acq_refine_blocks();
			acq.phase = (++acq.blk < blocks) ? PH_WAIT : PH_PEAK;
		}
		return kind;
	}
	acq.phase = PH_WAIT;
	acq.blk = 0;
	int cell;
	const float m = gps_acq_peak(acq.acc, cold ? 0 : acq.acc1, &cell);
	if (cold) {
		for (int r = 0; r < rows; r++) {
			const float v = r ? gps_acq_peak(acq.acc1, 0, &cell) : m;
			if (v > acq.best[acq.prn_i].metric) {
				acq.best[acq.prn_i].metric = v;
				acq.best[acq.prn_i].dop_i = acq.dop_i + 4 * r;
			}
		}
		/* The next Doppler not a pair's second: four in a run pair with the four after. */
		acq.dop_i += (((acq.dop_i + 1 - acq.lo[acq.prn_i]) & 4) ? 5 : 1);
		if (acq.dop_i > acq.hi[acq.prn_i]) {
			acq.prn_i++;
			acq_next();
		}
		return STEP_ACQ_CORR;
	}
	/* One Doppler, the cold pass's: 11 steps of 25 Hz across the bin cost 11 times the time,
	 * detected no more on a 5-minute recording, and picked noise for the Doppler. */
	if (m >= acq_threshold(acq_refine_blocks(), 8.0f)) {
		acq.ref_doppler = acq.doppler;
		acq.ref_cell = cell;
		acq.ref_first0 = acq.first0;
		acq_arm(prn);
	}
	if (++acq.cand_i == acq.ncand) {
		acq_pass_done();
	}
	return STEP_ACQ_CORR;
}

/* The cold pass over every PRN index from 0 at its Dopplers (acq.lo..acq.hi, predicted or all). */
static void acq_cold_pass(void)
{
	acq.mode = ACQ_COLD;
	acq.prn_i = 0;
	acq.blk = 0;
	acq.phase = PH_WAIT;
	acq.pass_armed = 0;
	for (int i = 0; i < 32; i++) {
		acq.best[i].metric = 0;
	}
	acq_next();
}

/* A pass is done: idle when it armed one and four channels track (a fix's satellites), else again
 * with twice the blocks (up to the most): the satellites left are weaker than those it found. */
static void acq_pass_done(void)
{
	int active = 0;
	for (int c = 0; c < CHANNELS; c++) {
		active += g.ch[c].active;
	}
	/* A predicted pass searches about a twelfth of the cells: it goes twice as deep. */
	if ((acq.pass_armed && (active >= 4)) ||
	    (acq.blocks >= (acq.predicted ? 2 * ACQ_MAX_BLOCKS : ACQ_MAX_BLOCKS))) {
		acq.mode = ACQ_IDLE;
		return;
	}
	acq.blocks *= 2;
	acq_begin(); /* predicted afresh: a minute on, a satellite's Doppler has moved tens of hertz */
}

/* Message 'C' blocks refine_blocks: the cold search, records on (setting 0 = 1). */
static uint32_t acq_start(uint32_t blocks, uint32_t refine_blocks)
{
	mgr.on = 1;
	if (blocks == 0) {
		return 0; /* the manager without a search */
	}
	if ((refine_blocks == 0) || (acq.mode != ACQ_IDLE) || !send_setting(0, 1)) {
		return 2;
	}
	acq.blocks = (uint8_t) blocks;
	acq.refine_blocks = (uint8_t) refine_blocks;
	acq.code_prn = 0;
	acq.code_op = 0;
	acq.q_n = 0;
	acq.narmed = 0;
	acq_begin();
	return 0;
}

/*
 * A pass. With an almanac, a held position (a survey's mean so far will do: a bin is 250 Hz, 100 km
 * about 80 Hz), an estimate's offset and a GPS time: each PRN predicted first (ACQ_PREDICT); else
 * every Doppler. The time is a usable channel's, which also sets the M4's clock of it (warm_t), or
 * that clock's from the last such or the host's ('W').
 */
static void acq_begin(void)
{
	acq.prn_i = 0;
	int c = 0;
	while ((c < CHANNELS) && !nav_usable(&g.ch[c])) {
		c++;
	}
	if (c < CHANNELS) {
		acq.warm_t = nav_transmit_time(&g.ch[c].nav, g.ch[c].nav.last[0].n0);
		acq.warm = 1;
	}
	/* No position: located from the satellites the manager holds locked, with an almanac entry
	 * each. */
	acq.loc_n = 0;
	for (int i = 0; i < 32; i++) {
		acq.lo[i] = ACQ_NDOP; /* none until a point sees it */
		acq.hi[i] = -1;
		acq.pred_hz[i] = 0;
	}
	for (c = 0; c < CHANNELS; c++) {
		const struct channel* ch = &g.ch[c];
		const struct gps_alm* a = &alm[(ch->prn - 1) & 31];
		if (mgr_locked(c) && a->have && !a->health) {
			acq.loc_k[acq.loc_n] = ch->prn - 1;
			acq.loc_m[acq.loc_n++] = (float) gps_loops_doppler(&ch->loops);
		}
	}
	acq.pred_t = acq.warm_t;
	acq.pred_end = 0;
	acq.predicted = 1;
	acq.loc_ok = 0;
	acq.loc_kept = 0;
	acq.loc_ring = 0;
	acq.loc_i = acq.loc_np = 0;
	acq.loc_cl = 0; /* the south pole */
	acq.loc_sl = -1;
	acq.mode = ACQ_LOCATE;
	if (fr.have_bias && fr.pos_state) {
		acq.loc_n = 0;
		acq.mode = ACQ_PREDICT;
	}
	/* With no almanac entry, every PRN at every Doppler (loc_done) after the grid's few seconds. */
	if (!acq.warm || ((acq.mode == ACQ_LOCATE) && !acq.loc_n)) {
		acq.loc_kept = 0;
		acq.predicted = 0;
		loc_done();
	}
}

static void put_u32(volatile uint8_t* p, uint32_t v);

/*
 * Search reply, 29 bytes: mode (0 idle, 1 cold, 2 refining, 3 predicting from the almanac, 4
 * locating from it), the PRN index and Doppler index in hand, the candidates (count, then up to
 * seven PRN indices), PRNs armed (count, then seven), the blocks a cold PRN and Doppler (byte 19),
 * 1 when the pass
 * searches only the Dopplers a whole almanac predicted (byte 20), then the search's ticks that did
 * an operation (4) and that waited for a snapshot (4): 29 bytes.
 */
static uint32_t acq_reply(void)
{
	volatile uint8_t* r = APP_REPLY_OUT;
	r[0] = acq.mode;
	r[1] = acq.prn_i;
	r[2] = (uint8_t) acq.dop_i;
	r[3] = acq.ncand;
	for (int i = 0; i < CHANNELS; i++) {
		r[4 + i] = acq.cand[i];
		r[12 + i] = acq.armed[i];
	}
	r[11] = acq.narmed;
	r[19] = acq.blocks;
	r[20] = acq.predicted;
	put_u32(&r[21], acq.ticks_work);
	put_u32(&r[25], acq.ticks_wait);
	return 29;
}

/* A tick of the search: counted as work, or as a wait for a snapshot (acq_reply). */
static int acq_tick(void)
{
	if (acq.mode == ACQ_IDLE) {
		return 0;
	}
	const int step = acq_search();
	if (step) {
		acq.ticks_work++;
	} else {
		acq.ticks_wait++;
	}
	return step;
}

/* ---- the satellite manager -------------------------------------------------------------- */

/*
 * Once a second (each frequency window's end), on from the first search (message 'C'). A GPS
 * channel whose symbol-sum lock stays under MGR_LOCK for MGR_LOST seconds, after MGR_SETTLE of
 * pull-in, is kept while its loops' Doppler stays within MGR_HELD Hz of its last locked second's,
 * plus 1 Hz a second unlocked (a satellite's Doppler moves under 1 Hz a second), up to MGR_HOLD
 * seconds; then stopped. Settled loops hold a fading satellite's frequency below the lock gate
 * and often lock again; a re-arm restarted the loops' pull-in, which wandered off a weak signal
 * in a second (none of six re-arms held, on the box). While a channel is free and the search
 * idle, a pass over the PRNs not tracked starts MGR_AGAIN seconds after the last.
 */
#define MGR_LOCK   30 /* hundredths */
#define MGR_SETTLE 15
#define MGR_LOST   10
#define MGR_HELD   20
#define MGR_HOLD   60
#define MGR_AGAIN  60

/* Channel c holds a GPS satellite past its pull-in, locked. */
static int mgr_locked(int c)
{
	const struct channel* ch = &g.ch[c];
	return ch->active && (ch->prn < 120) && (ch->age >= MGR_SETTLE) && (fr.lock[c] >= MGR_LOCK);
}

static void mgr_second(void)
{
	mgr.quiet = 0;
	/* Records a stopped channel sent before the stop took are still read as records for a
	 * second or two (parse_word). */
	g.stopped = g.stopping;
	g.stopping = 0;
	if (!mgr.on) {
		return;
	}
	int free = 0;
	for (int c = 0; c < CHANNELS; c++) {
		struct channel* ch = &g.ch[c];
		if (!ch->active) {
			free = 1;
			continue;
		}
		ch->age += ch->age < MGR_SETTLE; /* saturating: only the settling is counted */
		if ((ch->prn >= 120) || (ch->age < MGR_SETTLE)) {
			continue;
		}
		const float fd = (float) gps_loops_doppler(&ch->loops);
		ch->unlocked = (fr.lock[c] < MGR_LOCK) ? ch->unlocked + 1 : 0;
		if (!ch->unlocked) {
			ch->fd_lock = fd;
		}
		const float wander = fd > ch->fd_lock ? fd - ch->fd_lock : ch->fd_lock - fd;
		if ((ch->unlocked < MGR_LOST) ||
		    ((ch->unlocked < MGR_HOLD) && (wander < (float) (MGR_HELD + ch->unlocked)))) {
			continue;
		}
		if (command(OP_STOP, (uint32_t) c, 0, 0)) {
			ch->active = 0;
			g.stopping |= 1u << c;
		}
	}
	if (!free || (acq.mode != ACQ_IDLE)) {
		mgr.idle = 0;
	} else if (++mgr.idle >= MGR_AGAIN) {
		mgr.idle = 0;
		acq_start(acq.blocks, acq.refine_blocks);
	}
}

static void parse_word(uint16_t w)
{
	if (g.snap_skip) {
		g.snap_skip--;
		if (g.snap_skip < 3) {
			acq.trailer[2 - g.snap_skip] = w;
			if (g.snap_skip == 0) {
				acq_snap_done();
			}
		}
		return;
	}
	if (g.in_record) {
		g.rec[g.rec_n++] = w;
		if (g.rec_n == RECORD_WORDS) {
			g.in_record = 0;
			handle_record();
		}
		return;
	}
	// Records come only from channels the app armed, or stopped in the last second or two (their
	// last records may still be on the way). A snapshot starts at a record boundary; its first
	// word passes for such a channel's record start with a chance of about 1 in 65,536 (then the
	// snapshot is misread).
	const uint32_t c = w >> 8;
	const int looks_like_record = ((w & 0xFF) == MAGIC) && (c < CHANNELS) &&
		(g.ch[c].active || (((g.stopping | g.stopped) >> c) & 1));
	if (g.snap_pending && !looks_like_record) {
		g.snap_pending = 0;
		g.snap_skip = g.snap_words - 1;
		g.snaps++;
		acq_snap_start(g.pos);
		return;
	}
	if (looks_like_record) {
		g.rec[0] = w;
		g.rec_n = 1;
		g.in_record = 1;
		return;
	}
	g.garbage++;
}

static void read_stream(void)
{
	const uint32_t written = APP_M0_COUNT & ~1u;
	if (!g.synced) {
		g.pos = written;
		g.synced = 1;
		return;
	}
	uint32_t avail = written - g.pos;
	if (avail > APP_SAMPLES_SIZE - 2048) {
		// Fell behind the M0: start over at the newest data. The snapshots queued for the search
		// lay in the ring passed over.
		g.overruns++;
		g.pos = written;
		g.in_record = 0;
		g.snap_skip = 0;
		acq.q_n = 0;
		return;
	}
	const volatile uint16_t* words = (const volatile uint16_t*) APP_SAMPLES; /* pos stays even */
	for (uint32_t i = 0; (i < MAX_WORDS) && (avail >= 2); i++, avail -= 2) {
		if (g.snap_skip > 3) {
			/* Inside a snapshot: its samples stay in the ring for the search (acq_snap_start),
			 * so pass them at once, up to the trailer's three words. Word by word, 4,096 a
			 * block cost the search's ticks about 0.1 ms each (measured). */
			uint32_t k = g.snap_skip - 3;
			k = (k < avail / 2) ? k : avail / 2;
			g.snap_skip -= k;
			g.pos += 2 * k;
			avail -= 2 * k;
			if (avail < 2) {
				break;
			}
		}
		parse_word(words[(g.pos % APP_SAMPLES_SIZE) / 2]);
		g.pos += 2;
	}
}

/* ---- messages ------------------------------------------------------------------------- */

static uint32_t get_u32(const volatile uint8_t* p)
{ return p[0] | (p[1] << 8) | (p[2] << 16) | ((uint32_t) p[3] << 24); }

static void put_u32(volatile uint8_t* p, uint32_t v)
{
	p[0] = v;
	p[1] = v >> 8;
	p[2] = v >> 16;
	p[3] = v >> 24;
}

static uint32_t arm(const volatile uint8_t* m)
{
	return arm_channel(
		m[1],
		m[2],
		get_u32(&m[3]),
		get_u32(&m[7]),
		m[11] | (m[12] << 8),
		m[13] | (m[14] << 8),
		m[15],
		get_u32(&m[16]),
		(int32_t) get_u32(&m[20]));
}

static uint32_t arm_channel(
	uint32_t c,
	uint32_t prn,
	uint32_t W,
	uint32_t Wc,
	uint32_t g1,
	uint32_t g2,
	uint32_t win,
	uint32_t start,
	int32_t doppler_mhz)
{
	if (c >= CHANNELS) {
		return 2;
	}
	uint8_t bits[15] = {0};
	put_bits(bits, 0, W, 32);
	put_bits(bits, 32, Wc + (uint32_t) steer.word, 32); /* the steering, as in every word */
	put_bits(bits, 64, g1, 10);
	put_bits(bits, 74, g2, 10);
	put_bits(bits, 84, win, 3);
	put_bits(bits, 87, start, 32);
	struct channel* ch = &g.ch[c];
	const int symbol_ms = (prn >= 120 && prn <= 138) ? 2 : 20;
	struct gps_loops loops;
	if (gps_loops_init(&loops, FS, doppler_mhz * 1e-3, symbol_ms, 300, 1000, 0) !=
	    0) {
		return 2;
	}
	if (!command(OP_ARM, c, bits, 20)) {
		return 1;
	}
	*ch = (struct channel){0};
	ch->loops = loops;
	ch->prn = (uint8_t) prn;
	ch->code = W;
	ch->fd_lock = (float) doppler_mhz * 1e-3f;
	ch->carrier = (int32_t) Wc;
	ch->steer_in_use = steer.word; /* in the arm's word */
	/* The first GPS channel of a run starts the frequency windows' grid afresh. */
	int others = 0;
	for (uint32_t k = 0; k < CHANNELS; k++) {
		others |= (k != c) && g.ch[k].active && (g.ch[k].prn < 120);
	}
	if (!others) {
		fr.have_grid = 0;
		fr.state = 0;
	}
	/* A window's expected frequencies or a fix in progress that used this channel read its
	 * ephemeris, now zeroed (a host's 'A' on a tracking channel): they are dropped, not finished
	 * on garbage (a NaN in the estimate poisoned every satellite's weight). Such an 'A' also hands
	 * the new loops the record or two of the old arm still on the way: arm a stopped channel. */
	for (int i = 0; i < fr.n; i++) {
		if (fr.chs[i] == c) {
			fr.state = 0;
		}
	}
	if (nav.used & (1u << c)) {
		nav.state = 0;
	}
	ch->active = 1;
	return 0;
}

/*
 * Status reply, 128 bytes: per channel 16 bytes (its PRN, 0 when stopped; gaps (saturating), symbol
 * edge (255 until found), lost pending requests (saturating), records (4), Doppler in mHz (4),
 * |IP| + |QP| of the last record (4)), then records (4), garbage words (4), overruns (2), snapshots (2) and refused
 * register writes (4).
 */
static uint32_t status_reply(void)
{
	volatile uint8_t* r = APP_REPLY_OUT;
	for (int c = 0; c < CHANNELS; c++, r += 16) {
		const struct channel* ch = &g.ch[c];
		r[0] = ch->active ? ch->prn : 0;
		r[1] = ch->gaps > 255 ? 255 : ch->gaps;
		r[2] = (uint8_t) ch->loops.edge;
		r[3] = ch->lost > 255 ? 255 : ch->lost;
		put_u32(&r[4], ch->records);
		/* Hz * 1000 = carrier * fs * 1000 / 2^32, fs = 4 MHz */
		put_u32(&r[8], (uint32_t) ((ch->carrier * 4000000000LL) >> 32));
		put_u32(&r[12], (uint32_t) ch->prompt);
	}
	put_u32(&r[0], g.records);
	put_u32(&r[4], g.garbage);
	r[8] = g.overruns;
	r[9] = g.overruns >> 8;
	r[10] = g.snaps;
	r[11] = g.snaps >> 8;
	put_u32(&r[12], g.refused);
	return CHANNELS * 16 + 16;
}

/* Bytes as they lie: a double or float goes to the host little endian, as it reads them. */
static void copy_bytes(volatile void* dst, const volatile void* src, uint32_t n)
{
	for (uint32_t i = 0; i < n; i++) {
		((volatile uint8_t*) dst)[i] = ((const volatile uint8_t*) src)[i];
	}
}

/*
 * Navigation reply, 128 bytes: fixes kept (4), the last kept fix's sample index (4), satellites in
 * it, solve iterations, the channels used (bit mask), satellites usable at the last attempt; GPS
 * time of week at that sample (double), ECEF x, y, z in meters (doubles), G, P, H, V, T DOP
 * (floats); per channel 8 bytes: subframes (2), flags (bit 0 anchored, 1 ephemeris complete, 2
 * healthy, 3 reading a subframe), 0, the latest HOW count (4); then fixes rejected (4) and the last
 * attempt's result: 0 kept, 1 fewer than four usable satellites, 2 not converged in ten
 * iterations, 3 not near the Earth's surface, 4 GDOP above MAX_GDOP, 5 three satellites with the
 * height held, not kept but the held position until the first kept fix.
 */
static uint32_t nav_reply(void)
{
	volatile uint8_t* r = APP_REPLY_OUT;
	put_u32(&r[0], nav.fixes);
	put_u32(&r[4], nav.fix_n);
	r[8] = nav.fix_sats;
	r[9] = nav.fix_iterations;
	r[10] = nav.fix_used;
	r[11] = nav.tried_sats;
	copy_bytes(&r[12], &nav.fix_t, 8);
	copy_bytes(&r[20], nav.fix_x, 24);
	copy_bytes(&r[44], nav.fix_dop, 20);
	put_u32(&r[120], nav.rejected);
	r[124] = nav.result;
	r[125] = r[126] = r[127] = 0;
	for (int c = 0; c < CHANNELS; c++) {
		const struct chnav* v = &g.ch[c].nav;
		volatile uint8_t* e = &r[64 + 8 * c];
		e[0] = (uint8_t) v->subframes;
		e[1] = (uint8_t) (v->subframes >> 8);
		e[2] = (uint8_t) (v->anchored | (gps_nav_eph_complete(&v->eph) << 1) |
				  ((v->eph.health == 0) << 2) | ((v->nwords != 0) << 3));
		e[3] = 0;
		put_u32(&e[4], v->tow);
	}
	return 128;
}

/*
 * Frequency reply, 64 bytes: estimates made (4), the sample index the last window ended at (4; it
 * began FIX_SAMPLES before), the estimate in ppb (float: the oscillator's fractional error, so -L1
 * times it is the carrier offset common to the satellites), the used satellites' rms about it in
 * ppb (float); the channels with a usable window (bit mask), the channels used (bit mask), the last
 * window's result (0 estimated, 1 no position yet, 2 no usable satellite, 3 no window yet, 4 the
 * satellites disagree), the
 * held position (0 none, 1 surveying, 2 surveyed, 3 given); kept fixes in the survey (4); per
 * channel its offset in ppb (7 floats, where usable); per channel its symbol-sum phase lock in
 * hundredths (7 signed bytes, 0 for a channel with no window); then the steering's state
 * (bit 0 a ratio waits, bit 1 the last was refused, bit 2 the loop on) and the steering in effect
 * in ppb (float).
 * The offsets are the TCXO's; CLKOUT's is the estimate plus the steering.
 */
static uint32_t freq_reply(void)
{
	volatile uint8_t* r = APP_REPLY_OUT;
	put_u32(&r[0], fr.estimates);
	put_u32(&r[4], fr.window_end);
	copy_bytes(&r[8], &fr.est, 4);
	copy_bytes(&r[12], &fr.rms, 4);
	r[16] = fr.valid;
	r[17] = fr.used;
	r[18] = fr.result;
	r[19] = fr.pos_state;
	put_u32(&r[20], fr.surveyed);
	copy_bytes(&r[24], fr.ppb, 4 * CHANNELS);
	copy_bytes(&r[52], fr.lock, CHANNELS);
	r[59] = steer.want | (steer.refused << 1) | (servo.on << 2);
	const float d = (float) (steer.delta * 1e9);
	copy_bytes(&r[60], &d, 4);
	return 64;
}

/* ---- steering ---------------------------------------------------------------------------- */

static void steer_tick(uint32_t now)
{
	if (steer.log_fw) {
		uint32_t st[4];
		app_pll_b_status(st);
		if (!st[3]) {
			steer.log[(steer.writes - 1) % 4].start = st[1];
			steer.log[(steer.writes - 1) % 4].end = st[2];
			steer.log_fw = 0;
		}
	}
	if (!steer.want) {
		return;
	}
	const double delta = ((double) ((int32_t) steer.a - PLL_B_NOMINAL) * (double) steer.c +
			      (double) steer.b) /
		(PLL_B_NOMINAL * (double) steer.c);
	const double w = -steer.lo_hz * delta / FS * TWO32;
	const int32_t word = (int32_t) (w < 0 ? w - 0.5 : w + 0.5);
	const int large = (word - steer.word > LARGE_STEP) || (steer.word - word > LARGE_STEP);
	/*
	 * A large step waits for a tick with no request outstanding, so every channel's words follow
	 * within a period or two; a small one goes now, and a channel with a request outstanding
	 * follows when it is applied (under 0.02 ppb in one window at 7.9 Hz).
	 */
	uint32_t latest = 0;
	for (int c = 0; c < CHANNELS; c++) {
		if (large && g.ch[c].active && g.ch[c].outstanding) {
			return;
		}
		if (g.ch[c].active) {
			latest = g.ch[c].last_n;
		}
	}
	const uint32_t r = app_pll_b_set(steer.a, steer.b, steer.c);
	if (r == 1) {
		return;
	}
	steer.want = 0;
	steer.refused = r != 0;
	if (r != 0) {
		return;
	}
	steer.delta = delta;
	steer.word = word;
	for (int c = 0; c < CHANNELS; c++) {
		struct channel* ch = &g.ch[c];
		if (ch->active) {
			ch->send = 1;
			if (large) {
				ch->fw.bad = 1;
				ch->fw_hold = 1;
			}
		}
	}
	const uint32_t k = steer.writes++ % 4;
	steer.log[k].a = steer.a;
	steer.log[k].b = steer.b;
	steer.log[k].c = steer.c;
	steer.log[k].tick = now;
	steer.log[k].n = latest;
	steer.log[k].start = steer.log[k].end = 0;
	steer.log_fw = 1;
}

/*
 * PLL B writes, 116 bytes: the count (4), then the last four, the latest first, 28 bytes each: a,
 * b, c, the cycle counter at the tick that asked, the highest active channel's last period start
 * then, and the cycle counter before and after the firmware's write (0 until known).
 */
static uint32_t steer_reply(void)
{
	volatile uint8_t* r = APP_REPLY_OUT;
	put_u32(&r[0], steer.writes);
	for (uint32_t i = 0; i < 4; i++) {
		const uint32_t k = (steer.writes - 1 - i) % 4;
		copy_bytes(&r[4 + 28 * i], &steer.log[k], 28);
	}
	return 116;
}

/* ---- load -------------------------------------------------------------------------------- */

/*
 * The app cannot read the cycle counter (the system control space faults), but each tick's
 * argument is the counter at the call. Ticks a second, against the count with no channel armed,
 * give the share of the M4 the app and its register writes take, the writes being made in the
 * same main-loop pass; the longest gap between ticks is the worst delay a record can see.
 */
#define CPU_HZ 204000000u

static struct {
	uint8_t started;
	uint32_t window_start, last, ticks, max_gap;
	uint32_t done_ticks, done_cycles, done_max_gap; /* the last full window */
} load APP_STATE;

static void count_tick(uint32_t now)
{
	if (!load.started) {
		load.started = 1;
		load.window_start = now;
	} else if (now - load.last > load.max_gap) {
		load.max_gap = now - load.last;
	}
	load.last = now;
	load.ticks++;
	if (now - load.window_start >= CPU_HZ) {
		/* The manager runs at each frequency window's end; with no channel to end one (no GPS
		 * channel has records), on the M4's clock once three seconds pass without, then each
		 * second, or no search would start again. */
		if (++mgr.quiet >= 3) {
			mgr_second();
			mgr.quiet = 2;
		}
		acq.warm_t += 1.0; /* a second to within a tick */
		load.done_ticks = load.ticks;
		load.done_cycles = now - load.window_start;
		load.done_max_gap = load.max_gap;
		load.window_start = now;
		load.ticks = 0;
		load.max_gap = 0;
	}
}

/*
 * Ticks (4), the window's length in cycles (4), the longest gap between ticks in cycles (4), and the
 * stack's deepest use in bytes (4): the firmware zeroes the 4 KB stack at start, so the lowest
 * word ever written marks it (a word left zero by chance reads a little shallow).
 */
static uint32_t load_reply(void)
{
	volatile uint8_t* r = APP_REPLY_OUT;
	put_u32(&r[0], load.done_ticks);
	put_u32(&r[4], load.done_cycles);
	put_u32(&r[8], load.done_max_gap);
#ifndef APP_HOST
	const volatile uint32_t* sp = (const volatile uint32_t*) 0x10018000;
	uint32_t i = 0;
	while ((i < 1024) && (sp[i] == 0)) {
		i++;
	}
	put_u32(&r[12], 4096 - 4 * i);
#else
	put_u32(&r[12], 0);
#endif
	return 16;
}

/*
 * The warm start (message 'W'): what a run learned, which hackrf_gps keeps in a file and hands back
 * so the first search is a predicted one.
 *   'W'                  reply: the held position (3 doubles), its state, whether the offset below
 *                        is known, and the offset (a double: a tracked satellite's loops' Doppler
 *                        less its expected carrier, in hertz)
 *   'W' prn-1            reply: that PRN's almanac entry (struct gps_alm, its bytes)
 *   'W' prn-1 entry      set it
 *   'W' t(8) offset(8)   GPS time of week now (s) and the offset: predict with no usable channel
 */
static uint32_t warm_message(const volatile uint8_t* m, uint32_t len)
{
	volatile uint8_t* r = APP_REPLY_OUT;
	if (len == 1) {
		copy_bytes(r, fr.rx, 24);
		r[24] = fr.pos_state;
		r[25] = fr.have_bias;
		copy_bytes(&r[26], &fr.bias, 8);
		return 34;
	}
	if (((len == 2) || (len == 2 + sizeof(struct gps_alm))) && (m[1] < 32)) {
		if (len > 2) {
			copy_bytes(&alm[m[1]], &m[2], sizeof(struct gps_alm));
		}
		copy_bytes(r, &alm[m[1]], sizeof(struct gps_alm));
		return sizeof(struct gps_alm);
	}
	r[0] = 2;
	if (len == 17) {
		copy_bytes(&acq.warm_t, &m[1], 8);
		copy_bytes(&fr.bias, &m[9], 8);
		fr.have_bias = 1;
		acq.warm = 1;
		r[0] = 0;
	}
	return 1;
}

static uint32_t handle_message(uint32_t len)
{
	const volatile uint8_t* m = APP_MESSAGE_IN;
	uint32_t status = 2;
	switch (m[0]) {
	case 'M':
		status = (len == 2) ? (send_setting(0, m[1]) ? 0 : 1) : 2;
		break;
	case 'S':
		// Taken while the previous snapshot is still passing (its last words reach the
		// sample buffer only when more data follows), but not while one waits to start.
		if ((len == 2) && !g.snap_pending) {
			status = send_setting(1, m[1]) ? 0 : 1;
			if (status == 0) {
				g.snap_pending = 1;
				g.snap_words = (m[1] + 1u) * 1024 + 3;
			}
		}
		break;
	case 'A':
		status = (len == 24) ? arm(m) : 2;
		break;
	case 'X':
		if ((len == 2) && (m[1] < CHANNELS)) {
			status = command(OP_STOP, m[1], 0, 0) ? 0 : 1;
			if (status == 0) {
				g.ch[m[1]].active = 0;
				g.stopping |= 1u << m[1];
			}
		}
		break;
	case 'Q':
		return status_reply();
	case 'L':
		return load_reply();
	case 'P':
		return nav_reply();
	case 'F':
		if ((len == 5) || (len == 9)) {
			const double off = (int32_t) get_u32(&m[1]) * 1e-6;
			nav_lo_word = off / FS * TWO32;
			steer.lo_hz = 1575.42e6 + off - (len == 9 ? (int32_t) get_u32(&m[5]) : 0);
			status = 0;
		}
		break;
	case 'R':
		return freq_reply();
	case 'J':
		return steer_reply();
	case 'C':
		status = (len == 3) ? acq_start(m[1], m[2]) : 2;
		break;
	case 'K':
		return acq_reply();
	case 'D':
		if (len == 2) {
			servo.on = m[1] != 0;
			servo.have = 0;
			servo.n = 0;
			servo.target = steer.delta * 1e9;
			status = 0;
		}
		break;
	case 'B':
		if (len == 13) {
			status = 1;
			if (!steer.want) {
				steer.a = get_u32(&m[1]);
				steer.b = get_u32(&m[5]);
				steer.c = get_u32(&m[9]);
				steer.want = 1;
				status = 0;
			}
		}
		break;
	case 'U':
		if ((len == 5) || (len == 29)) {
			fr.survey_n = get_u32(&m[len - 4]);
			fr.surveyed = 0;
			fr.sum[0] = fr.sum[1] = fr.sum[2] = 0;
			fr.pos_state = 0;
			status = 0;
		}
		if ((len == 25) || (len == 29)) {
			copy_bytes(fr.rx, &m[1], 24); /* little endian, as the host sends them */
			fr.pos_state = (len == 25) ? 3 : 1; /* a seed: the first kept fix replaces it */
			fr.from3 = 0;
			status = 0;
		}
		break;
	case 'W':
		return warm_message(m, len);
	case 'H':
		if ((len == 5) || (len == 6)) {
			copy_bytes(&fr.height, &m[1], 4);
			fr.rough = (len == 6) && m[5];
			status = 0;
		}
		break;
	}
	APP_REPLY_OUT[0] = status;
	return 1;
}

/* ---- entry -------------------------------------------------------------------------------- */

uint32_t app_main(uint32_t op, uint32_t arg)
{
	switch (op) {
	case APP_OP_INIT:
		/* The firmware zeroes the data region at the start (app.h); only what is not zero. */
		nav.result = 1; /* no attempt yet: no usable satellite */
		fr.survey_n = SURVEY_FIXES;
		fr.result = 3;
		steer.lo_hz = 1575.42e6;
		gps_acq_init();
		return 0;
	case APP_OP_TICK:
		count_tick(arg);
		read_stream();
		steer_tick(arg);
		/* One step of navigation, else of the estimate, else of the search, a tick. */
		if ((nav_tick() == 0) && (freq_tick() == 0)) {
			acq_tick();
		}
		for (int c = 0; c < CHANNELS; c++) {
			struct channel* ch = &g.ch[c];
			if (ch->active && ch->send && !ch->outstanding &&
			    send_pending(c, (uint32_t) ch->code, (uint32_t) ch->carrier + steer.word)) {
				ch->send = 0;
				ch->outstanding = 1;
				ch->waited = 0;
				ch->sent_W = (uint32_t) ch->code;
				ch->sent_Wc = (uint32_t) ch->carrier + steer.word;
				ch->sent_steer = steer.word;
			}
		}
		flush_commands();
		return 0;
	case APP_OP_MESSAGE: {
		const uint32_t n = handle_message(arg);
		flush_commands();
		return n;
	}
	default:
		return 0;
	}
}
