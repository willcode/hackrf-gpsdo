#
# SPDX-License-Identifier: BSD-3-Clause

"""Bit-exact model of the integer loops: the floating-point loops of gps_track_model.Loops in the
integer arithmetic the M4 app runs (gps/app/gps_loops.c, checked against this line for line).

The loops see one channel's sums and write its two NCO words. Every operation here is one the C
performs, so it must reproduce these words exactly.

- Angles: a 16-iteration CORDIC in vectoring mode on the sums shifted right by PRE_SHIFT (a 20 ms sum
  of products of 8-bit samples and the 8-bit carrier table is below 2^31; the shift keeps the CORDIC's
  1.65 gain inside 32 bits). The angle unit is a full circle / 2^16; the vector's length comes out
  scaled by the CORDIC gain, the same for early and late, so their ratio is unaffected.
- Costas error: the angle folded into (-1/4, 1/4] of a circle (atan(Q/I)).
- FLL error: the difference of two angles, wrapped into [-1/2, 1/2) and, for the two-quadrant FLL on
  1 ms pairs, folded into (-1/4, 1/4]. A sign change between neighbors (for the symbol edge) is a
  difference outside [-1/4, 1/4].
- DLL error: tau = ((|E| - |L|) << 14) / (|E| + |L|), a division truncated toward zero: chips in Q15,
  with the early-late spacing of 1 chip already halved in.
- Frequency state: the carrier word in Q16 (48 bits), so gains below one word per angle unit exist.
  Every gain is an integer the host computes from a bandwidth (gains()), at most 16 bits wide.
- Code word: W_NOM + ((carrier word * AID) >> 24) + ((K_dll * tau) >> 15), AID = round(2^24 / 1540)
  being the carrier aiding (L1 is 1540 chip rates).

The symbol-edge histogram is the host's part, kept here so the stages run as in Loops.
"""

import math

from . import gps_model as g

ANG_BITS = 16
FULL = 1 << ANG_BITS
HALF = FULL >> 1
QUARTER = FULL >> 2
CORDIC_ITERS = 16
ATAN = [round(math.atan(2.0 ** -i) / (2 * math.pi) * FULL) for i in range(CORDIC_ITERS)]
PRE_SHIFT = 2
Q = 16                                   # fraction bits of the frequency state
AID = round((1 << 24) / 1540)


def wrap(a):
    """An angle into [-1/2, 1/2) of a circle."""
    return ((a + HALF) & (FULL - 1)) - HALF


def fold(a):
    """A wrapped angle into (-1/4, 1/4]: the two-quadrant angle, blind to a sign flip."""
    if a > QUARTER:
        return a - HALF
    if a <= -QUARTER:
        return a + HALF
    return a


def cordic(x, y):
    """(angle, length times the CORDIC gain) of the vector (x, y), vectoring mode."""
    x >>= PRE_SHIFT
    y >>= PRE_SHIFT
    a = 0
    if x < 0:
        x, y, a = -x, -y, HALF
    for i in range(CORDIC_ITERS):
        if y > 0:
            x, y, a = x + (y >> i), y - (x >> i), a + ATAN[i]
        else:
            x, y, a = x - (y >> i), y + (x >> i), a - ATAN[i]
    return wrap(a), x


def div_q14(n, d):
    """(n << 14) / d truncated toward zero; 0 when d is 0."""
    if d == 0:
        return 0
    q = (abs(n) << 14) // abs(d)
    return q if (n >= 0) == (d > 0) else -q


def gains(fs, symbol_ms, fll_bw=8.0, sync_bw=3.0, assist_bw=1.0, dll_bw=1.0, pll_bw=15.0,
          symbol_bt=0.16):
    """The integer loops' gains for one channel, from the bandwidths Loops takes."""
    w = (1 << 32) / fs                   # carrier word per hertz
    T = symbol_ms * 1e-3
    w0 = min(pll_bw, symbol_bt / T) / 0.53
    return dict(
        fll=round(4 * fll_bw * w), sync=round(4 * sync_bw * w),
        pll_i=round(w0 * w0 * T * w), pll_p=round(1.414 * w0 * w),
        assist=round(8 * assist_bw * w), dll=round(4 * dll_bw * w))


class IntLoops:
    """The integer loops for one channel; same stages and interface as gps_track_model.Loops."""

    ASSIST_MIN_MS = 10

    def __init__(self, fs, doppler, symbol_ms=20, fll_ms=300, sync_ms=1000, **bw):
        self.k = gains(fs, symbol_ms, **bw)
        for name, v in self.k.items():
            assert 0 <= v < (1 << 16), f"gain {name} = {v} does not fit 16 bits"
        self.w_nom = round(g.CHIP_RATE / fs * (1 << 32))
        self.fs = fs
        self.freq = round(doppler / fs * (1 << 32)) << Q      # carrier word, Q16
        self.integ = self.freq
        self.corr = 0                                         # DLL code-word correction
        self.symbol_ms = symbol_ms
        self.fll_ms = fll_ms
        self.sync_ms = sync_ms
        self.prev_angle = None
        self.count = 0
        self.flips = [0] * symbol_ms
        self.edge = None
        self.acc = None
        self.half = None
        self.words = None

    @property
    def doppler(self):
        return (self.freq >> Q) * self.fs / (1 << 32)

    def _dll(self, d):
        _, me = cordic(d["IE"], d["QE"])
        _, ml = cordic(d["IL"], d["QL"])
        tau = div_q14(me - ml, me + ml)
        self.corr = (self.k["dll"] * tau) >> 15

    def _words(self):
        carrier = self.freq >> Q
        code = self.w_nom + ((carrier * AID) >> 24) + self.corr
        self.words = (code, carrier)
        return self.words

    def update(self, d):
        idx = self.count
        self.count += 1
        if self.edge is None:
            angle, _ = cordic(d["IP"], d["QP"])
            prev, self.prev_angle = self.prev_angle, angle
            if prev is not None:
                diff = wrap(angle - prev)
                gain = self.k["fll"] if idx < self.fll_ms else self.k["sync"]
                self.freq += gain * fold(diff)
                self.integ = self.freq
                if idx >= self.fll_ms and not -QUARTER <= diff <= QUARTER:
                    self.flips[idx % self.symbol_ms] += 1
            if idx + 1 >= self.fll_ms + self.sync_ms:
                self.edge = max(range(self.symbol_ms), key=lambda i: (self.flips[i], -i))
            self._dll(d)
            return self._words()
        pos = (idx - self.edge) % self.symbol_ms
        if pos == 0:
            self.acc = {k: 0 for k in ("IE", "QE", "IP", "QP", "IL", "QL")}
            self.half = [[0, 0], [0, 0]]
        if self.acc is not None:
            for k in self.acc:
                self.acc[k] += d[k]
            h = self.half[pos >= self.symbol_ms // 2]
            h[0] += d["IP"]
            h[1] += d["QP"]
            if pos == self.symbol_ms - 1:
                if self.symbol_ms >= self.ASSIST_MIN_MS:
                    a0, _ = cordic(*self.half[0])
                    a1, _ = cordic(*self.half[1])
                    self.integ += self.k["assist"] * wrap(a1 - a0)
                e, _ = cordic(self.acc["IP"], self.acc["QP"])
                e = fold(e)
                self.integ += self.k["pll_i"] * e
                self.freq = self.integ + self.k["pll_p"] * e
                self._dll(self.acc)
                self.acc = None
                return self._words()
        return self.words
