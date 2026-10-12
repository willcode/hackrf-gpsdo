#
# SPDX-License-Identifier: BSD-3-Clause

"""Model of one correlator channel and its tracking loops, laid out as the gateware and the app are.

The channel's arithmetic is the gateware's: 32-bit code and carrier NCOs stepped once a sample, a
64-entry signed 8-bit carrier table, 8-bit samples, and early, prompt and late I and Q sums over one
code period. The loops here are floating point; their integer form, which the app runs, is
gps_loop_int.

Channel, per sample k of a code period (the period's first sample is k = 0):

- code: acc = F + k W, with W = round(chip rate / fs * 2^32). The prompt chip is acc >> 32; early is
  the chip of acc + 2^31, late of acc - 2^31 (one chip apart). The period ends at the first k where
  the prompt chip would reach 1023, and the next period starts with F = acc - 1023 * 2^32 there.
- carrier: theta = T + k Wc (mod 2^32); the table entry is theta >> 26.
  yi = xi C + xq S, yq = xq C - xi S, with C and S the table's cos and sin.
- sums: I and Q of each of E, P and L add yi and yq times the chip's sign.

The loops run on each period's sums (a dump). Their new code and carrier words take effect at the
start of the period after the next, the engine's latency: the loops compute during one period and
the channel picks up the result at its next period boundary.
"""

import numpy as np

from . import gps_model as g

NCO_BITS   = 32
NCO_ONE    = 1 << NCO_BITS
TABLE_BITS = 6
TABLE_AMP  = 127

_angles = 2 * np.pi * (np.arange(1 << TABLE_BITS) + 0.5) / (1 << TABLE_BITS)
COS = np.round(TABLE_AMP * np.cos(_angles)).astype(np.int64)
SIN = np.round(TABLE_AMP * np.sin(_angles)).astype(np.int64)


def code_word(fs, chip_rate):
    return int(round(chip_rate / fs * NCO_ONE))


def carrier_word(fs, hz):
    return int(round(hz / fs * NCO_ONE))


class Channel:
    """One channel's state and its per-period correlation."""

    def __init__(self, prn, fs, start, code_rate, doppler):
        self.chips = 1 - 2 * g.ca_code(prn).astype(np.int64)
        self.fs = fs
        self.n0 = start                  # first sample of the current code period
        self.F = 0                       # code NCO fraction at n0
        self.T = 0                       # carrier NCO phase at n0
        self.W = code_word(fs, code_rate)
        self.Wc = carrier_word(fs, doppler)
        self.pending = None              # (W, Wc) to take at the next period boundary

    def period(self, x):
        """Correlate one code period starting at n0; returns the dump, or None past the data."""
        W, Wc = self.W, self.Wc
        k_end = -(-(1023 * NCO_ONE - self.F) // W)
        if self.n0 + k_end > len(x):
            return None
        k = np.arange(k_end, dtype=np.int64)
        acc = self.F + k * W
        p = self.chips[(acc >> NCO_BITS) % 1023]
        e = self.chips[((acc + (NCO_ONE >> 1)) >> NCO_BITS) % 1023]
        l = self.chips[((acc - (NCO_ONE >> 1)) >> NCO_BITS) % 1023]
        theta = (self.T + k * Wc) & (NCO_ONE - 1)
        idx = theta >> (NCO_BITS - TABLE_BITS)
        c, s = COS[idx], SIN[idx]
        xs = x[self.n0:self.n0 + k_end]
        xi = xs.real.astype(np.int64)
        xq = xs.imag.astype(np.int64)
        yi = xi * c + xq * s
        yq = xq * c - xi * s
        dump = dict(n0=self.n0, F=self.F, W=W, Wc=Wc, samples=k_end,
                    IE=int(yi @ e), QE=int(yq @ e), IP=int(yi @ p), QP=int(yq @ p),
                    IL=int(yi @ l), QL=int(yq @ l))
        self.n0 += k_end
        self.F = self.F + k_end * W - 1023 * NCO_ONE
        self.T = (self.T + k_end * Wc) & (NCO_ONE - 1)
        if self.pending is not None:
            self.W, self.Wc = self.pending
            self.pending = None
        return dump


def arm_values(prn):
    """What the host sends to start a channel at chip 0 (dsp/gps_engine.py, config kind 0).

    g1, g2: the generators after producing chips 0 and 1, so their next output is chip 2 (bit i is
    stage i + 1). win: c(1022), c(0), c(1) as bits 0, 1, 2.
    """
    v = g.g2_start(prn)
    g1, g2 = [1] * 10, [(v >> i) & 1 for i in range(10)]
    for _ in range(2):
        f1 = g1[2] ^ g1[9]
        f2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9]
        g1, g2 = [f1] + g1[:9], [f2] + g2[:9]
    code = g.ca_code(prn)
    win = int(code[1022]) | (int(code[0]) << 1) | (int(code[1]) << 2)
    pack = lambda bits: sum(b << i for i, b in enumerate(bits))
    return pack(g1), pack(g2), win


class Loops:
    """The loops, floating point, in three stages, one dump at a time:

    1. FLL on 1 ms prompt pairs (two-quadrant, blind to data flips), `fll_bw`, for `fll_ms`, to pull in
       the acquisition's Doppler error;
    2. the same FLL at `sync_bw` for `sync_ms`, while sign changes between neighboring prompts are
       counted at each position within a data symbol: the most frequent position is the symbol edge
       (the host's bit sync, which it will send to the gateware);
    3. FLL-assisted Costas PLL on sums over whole symbols (20 ms GPS, 2 ms SBAS). The PLL's bandwidth
       keeps bandwidth times update interval at or below SYMBOL_BT. The assist is a four-quadrant FLL on
       the prompt sums of the symbol's two halves, which no data flip separates, at `assist_bw`; it
       holds the frequency within the PLL's pull-in. It runs only for symbols of at least ASSIST_MIN_MS:
       on SBAS's 2 ms symbols the halves are 1 ms, and the assist's noise broke the message decoding
       (1 of 42 WAAS messages passed CRC). New words come out once a symbol.

    A 1 ms Costas PLL is not used: near 30 dB-Hz it slipped to a false lock 30-70 Hz off on synthetic
    signals. An assist on 1 ms prompt pairs was too noisy there; one on half-symbol sums has 10 ms of
    signal in each term. The DLL (noncoherent early minus late, normalized, carrier aided) runs on the same sums as
    the carrier loop.
    """

    SYMBOL_BT = 0.16
    ASSIST_MIN_MS = 10

    def __init__(self, fs, doppler, symbol_ms=20, fll_ms=300, sync_ms=1000, pll_bw=15.0,
                 fll_bw=8.0, sync_bw=3.0, assist_bw=1.0, dll_bw=1.0):
        self.fs = fs
        self.doppler = doppler           # carrier NCO frequency, Hz
        self.integ = doppler             # PLL integrator, Hz
        self.code_rate_corr = 0.0        # DLL correction, chips/s
        self.symbol_ms = symbol_ms
        self.fll_ms = fll_ms
        self.sync_ms = sync_ms
        self.pll_bw = pll_bw
        self.fll_bw = fll_bw
        self.sync_bw = sync_bw
        self.assist_bw = assist_bw
        self.dll_bw = dll_bw
        self.prev = None
        self.count = 0                   # dumps seen
        self.flips = np.zeros(symbol_ms)
        self.edge = None                 # dump index mod symbol_ms at which a symbol starts
        self.acc = None
        self.half = None
        self.words = None

    def _pll(self, ip, qp, T, bw):
        w0 = bw / 0.53
        e = np.arctan(qp / ip) if ip != 0 else 0.0          # Costas, radians
        self.integ += w0 ** 2 * T * e / (2 * np.pi)
        self.doppler = self.integ + 1.414 * w0 * e / (2 * np.pi)

    def _dll(self, d):
        ae = np.hypot(d["IE"], d["QE"])
        al = np.hypot(d["IL"], d["QL"])
        tau = (ae - al) / (ae + al) / 2 if ae + al > 0 else 0.0   # chips; E-L spacing 1 chip
        self.code_rate_corr = 4 * self.dll_bw * tau

    def _words(self):
        code_rate = g.CHIP_RATE * (1 + self.doppler / g.L1_HZ) + self.code_rate_corr
        self.words = (code_word(self.fs, code_rate), carrier_word(self.fs, self.doppler))
        return self.words

    def update(self, d):
        idx = self.count
        self.count += 1
        ip, qp = d["IP"], d["QP"]
        prev, self.prev = self.prev, (ip, qp)
        cross = dot = 0
        if prev is not None and self.edge is None:
            cross = prev[0] * qp - prev[1] * ip
            dot = prev[0] * ip + prev[1] * qp
        if self.edge is None:
            bw = self.fll_bw if idx < self.fll_ms else self.sync_bw
            if dot != 0:
                self.doppler += 4 * bw * 1e-3 * np.arctan(cross / dot) / (2 * np.pi * 1e-3)
                self.integ = self.doppler
            if idx >= self.fll_ms and dot < 0:
                self.flips[idx % self.symbol_ms] += 1
            if idx + 1 >= self.fll_ms + self.sync_ms:
                self.edge = int(np.argmax(self.flips))
            self._dll(d)
            return self._words()
        pos = (idx - self.edge) % self.symbol_ms
        if pos == 0:
            self.acc = {k: 0 for k in ("IE", "QE", "IP", "QP", "IL", "QL")}
            self.half = [0j, 0j]
        if self.acc is not None:
            for k in self.acc:
                self.acc[k] += d[k]
            self.half[pos >= self.symbol_ms // 2] += complex(ip, qp)
            if pos == self.symbol_ms - 1:
                T = self.symbol_ms * 1e-3
                h = self.half[1] * np.conj(self.half[0])
                if h != 0 and self.symbol_ms >= self.ASSIST_MIN_MS:
                    df = np.angle(h) / (2 * np.pi * T / 2)
                    self.integ += 4 * self.assist_bw * T * df
                self._pll(self.acc["IP"], self.acc["QP"], T, min(self.pll_bw, self.SYMBOL_BT / T))
                self._dll(self.acc)
                self.acc = None
                return self._words()
        return self.words


def track(x, fs, prn, start, doppler, periods=None, integer=False, latency=0, **loop_kw):
    """Track one satellite from a code epoch at sample `start` with an initial Doppler.

    integer: the bit-exact integer loops (gps_loop_int.IntLoops) instead of Loops.
    latency: periods by which each update arrives later than the correlator's deadline (a dump's
    words otherwise take effect at the period after next); updates keep their order, none is lost.
    Returns a list of dumps, each with the loops' Doppler after it.
    """
    from .gps_loop_int import IntLoops
    code_rate = g.CHIP_RATE * (1 + doppler / g.L1_HZ)
    ch = Channel(prn, fs, start, code_rate, doppler)
    loops = (IntLoops if integer else Loops)(fs, doppler, symbol_ms=g.codes_per_symbol(prn), **loop_kw)
    dumps, late = [], []
    while periods is None or len(dumps) < periods:
        d = ch.period(x)
        if d is None:
            break
        late.append(loops.update(d))
        if len(late) > latency:
            ch.pending = late.pop(0)
        d["doppler"] = loops.doppler
        d["edge"] = -1 if loops.edge is None else loops.edge
        dumps.append(d)
    return dumps


def cn0_dbhz(dumps, T=1e-3):
    """C/N0 from the prompt sums by the second and fourth moments (insensitive to carrier phase)."""
    p = np.array([complex(d["IP"], d["QP"]) for d in dumps])
    m2 = np.mean(np.abs(p) ** 2)
    m4 = np.mean(np.abs(p) ** 4)
    pc = np.sqrt(max(2 * m2 ** 2 - m4, 0.0))
    pn = m2 - pc
    if pc <= 0 or pn <= 0:
        return float("nan")
    return float(10 * np.log10(pc / pn / T))


def cn0_windows(dumps, window=1000, T=1e-3):
    """cn0_dbhz over successive windows of dumps: the median and the list. A long span breaks the
    estimator's assumption of steady power when the signal fades."""
    vals = [cn0_dbhz(dumps[a:a + window], T) for a in range(0, len(dumps) - window + 1, window)]
    good = [v for v in vals if np.isfinite(v)]
    return (float(np.median(good)) if good else float("nan")), vals


def phase_lock(dumps, group=1, edge=0):
    """Phase-lock indicator (sum I^2 - Q^2) / (sum I^2 + Q^2) of prompt sums over `group` dumps from dump
    `edge` on (whole symbols when group is the symbol length); 1 is locked. Short sums understate lock
    on weak signals, where noise fills Q."""
    p = np.array([complex(d["IP"], d["QP"]) for d in dumps[edge:]])
    n = len(p) // group
    p = p[:n * group].reshape(n, group).sum(axis=1)
    i, q = p.real, p.imag
    return float((np.sum(i * i) - np.sum(q * q)) / (np.sum(i * i) + np.sum(q * q)))
