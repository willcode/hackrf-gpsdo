#!/usr/bin/env python3
#
# SPDX-License-Identifier: BSD-3-Clause

"""The oscillator's frequency offset from each satellite's carrier, against a held position: the
reference for the GPS app's estimator (gps/app/gps_freq.c).

Per window of records and satellite: the carrier cycles the engine applied (each record's carrier
word times the samples to the next record's start, summed), plus the tuner's offset, give the
received carrier frequency. The ephemeris and the held position give the frequency expected with
a perfect oscillator: the range rate with the Earth's rotation, and the satellite clock's rate, by
differences over transmit times one second apart. Measured less expected is -L1 times the
oscillator's fractional error, which every satellite shares; the window's estimate is the
satellites' mean weighted by each one's recent scatter, those more than REJECT_PPB from the median
left out (Combiner).

With the oscillator common to the tuner and the sample clock (no steering), the fractional error
is the same in both and -L1 times it is exact to parts in 1e7 of itself.
"""

import argparse
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import fix  # noqa: E402
import nav  # noqa: E402

L1 = 1575.42e6
ONE = float(1 << 32)
H = 0.5           # half the transmit-time interval of the expected frequency's difference, s
GAP = 12500       # samples between records beyond which a window is not used: two missing
# A satellite whose offset is this far from the window's median is left out: the satellites agree
# to about 0.11 ppb per second when locked, so 1 ppb is many times their scatter and far under what
# a lost lock gives.
REJECT_PPB = 1.0
MIN_LOCK = 0.5    # symbol-sum phase lock below which a window is not used (locked: 0.90-0.98)


def expected(eph, t, rx):
    """Carrier frequency (Hz) from satellite eph received at ECEF rx with a perfect oscillator, at
    transmit time t: L1 times (the satellite clock's rate) times (transmit time's rate in receive
    time) less L1."""
    off, clk = [], []
    for s in (-1.0, 1.0):
        pos, c = nav.sat_position(eph, t + s * H)
        rot = nav.OMEGA_E * np.linalg.norm(pos - rx) / nav.C
        sat = np.array([pos[0] * np.cos(rot) + pos[1] * np.sin(rot),
                        -pos[0] * np.sin(rot) + pos[1] * np.cos(rot), pos[2]])
        off.append(s * H + np.linalg.norm(sat - rx) / nav.C)   # receive time less t
        clk.append(c)
    return L1 * ((2 * H + clk[1] - clk[0]) / (off[1] - off[0]) - 1.0)


def phase_lock(ip, qp, first, edge):
    """(sum I^2 - sum Q^2) / (sum I^2 + sum Q^2) over the whole 20 ms symbol sums of records
    [first, first + len(ip)), edge being the record index of a symbol's start mod 20."""
    a = (edge - first) % nav.BIT_MS
    n = (len(ip) - a) // nav.BIT_MS
    if n < 1:
        return 0.0
    i = ip[a:a + n * nav.BIT_MS].reshape(n, nav.BIT_MS).sum(axis=1)
    q = qp[a:a + n * nav.BIT_MS].reshape(n, nav.BIT_MS).sum(axis=1)
    return float((np.sum(i * i) - np.sum(q * q)) / (np.sum(i * i) + np.sum(q * q)))


def satellite_offsets(z, s, boundaries, rx, lo_offset):
    """Per window [boundaries[j], boundaries[j + 1]): the satellite's offset in ppb (nan where the
    window is not usable), its measured frequency (Hz) and its phase lock. A window runs from the
    first record starting at or after one boundary to the first at or after the next."""
    prn = s["prn"]
    n0, Wc = z[f"{prn}_n0"], z[f"{prn}_Wc"]
    ip, qp = z[f"{prn}_IP"], z[f"{prn}_QP"]
    fs = float(z["fs"])
    edge = s["subframes"][0]["dump"] % nav.BIT_MS
    lo_word = round(lo_offset * 1e6) * 1e-6 / fs * ONE
    idx = np.searchsorted(n0, boundaries, side="left")
    out = np.full((len(boundaries) - 1, 3), np.nan)
    for j in range(len(boundaries) - 1):
        a, b = int(idx[j]), int(idx[j + 1])
        if b >= len(n0) or b - a < 500 or np.any(np.diff(n0[a:b + 1]) > GAP):
            continue
        dn = np.diff(n0[a:b + 1]).astype(np.float64)
        span = float(n0[b] - n0[a])
        cycles = float(np.sum(dn * Wc[a:b])) + span * lo_word
        f = cycles / span * fs / ONE
        mid = (n0[a] + n0[b]) // 2
        t = fix.transmit_time(z, prn, s["subframes"][0], mid)
        lock = phase_lock(ip[a:b], qp[a:b], a, edge)
        ppb = -(f - expected(s["eph"], t, rx)) / L1 * 1e9
        out[j] = (ppb if lock >= MIN_LOCK else np.nan, f, lock)
    return out


PRIOR_VAR = 0.1 ** 2  # a satellite's variance before any history, ppb^2
FLOOR = 0.01          # no satellite's standard deviation counts as less, ppb
VAR_K = 16            # the variance follows the squared residual with weight 1/VAR_K


class Combiner:
    """The common offset window by window, as gps_freq.c does it: the weighted mean of the
    satellites within REJECT_PPB of the median, each weighted by the inverse of its offset's recent
    variance about the estimate (PRIOR_VAR before any history, never under FLOOR^2); every usable
    satellite's variance then follows its squared residual, limited to REJECT_PPB. A NaN or infinite
    offset takes no part. None is used, and no variance changes, when none is within REJECT_PPB of
    the median."""

    def __init__(self):
        self.var = {}

    def __call__(self, offsets):
        """offsets: {key: ppb} of the usable satellites. Returns (estimate, keys used)."""
        offsets = {k: v for k, v in offsets.items() if np.isfinite(v)}
        if not offsets:
            return np.nan, []
        med = np.median(list(offsets.values()))
        used = [k for k, v in offsets.items() if abs(v - med) <= REJECT_PPB]
        if not used:  # two apart by more than twice REJECT_PPB: nothing says which is wrong
            return np.nan, []
        w = {k: 1 / max(self.var.get(k, 0) or PRIOR_VAR, FLOOR ** 2) for k in used}
        est = sum(w[k] * offsets[k] for k in used) / sum(w.values())
        for k, v in offsets.items():
            r = min(max(v - est, -REJECT_PPB), REJECT_PPB)
            old = self.var.get(k, 0) or PRIOR_VAR
            self.var[k] = old + (r * r - old) / VAR_K
        return est, used


def predict(t, e, at):
    """The least-squares line through estimates e at times t, at time at; the latest when fewer
    than three (gps_freq_predict, operation for operation)."""
    n = len(t)
    if n < 3:
        return e[-1]
    tm, em = sum(t) / n, sum(e) / n
    stt = sum((x - tm) * (x - tm) for x in t)
    ste = sum((x - tm) * (y - em) for x, y in zip(t, e))
    return em + ste / stt * (at - tm)


class Servo:
    """The GPS app's steering loop (gps_app.c, servo_estimate): each estimate of the TCXO's error
    (ppb) at a window's end joins the last N; PLL B's steering is the negative of the line through
    them at the middle of the coming second, within MAX_STEP of the last asked for, as the nearest
    ratio a + b/c (c up to 2^20 - 1), not asked for when it changes the steering less than
    MIN_CHANGE. A target in the gap where no ratio lies (within GAP of nominal) becomes whichever of
    -GAP, 0 and +GAP is nearest the line's own. After a pause longer than N windows the line starts
    over. Returns the ratio asked for, or None."""

    N, HORIZON, MAX_STEP, MIN_CHANGE, MAX, MAX_C = 10, 1.0, 5.0, 0.002, 4000.0, 1048575
    GAP = 29.80232

    def __init__(self, steering_ppb=0.0, fs=4e6):
        self.t, self.e, self.fs = [], [], fs
        self.target = steering_ppb
        self.t_last = self.end_last = None

    def __call__(self, end, est):
        d = None
        if self.end_last is not None:
            d = (end - self.end_last) & 0xFFFFFFFF
            d = d - (1 << 32) if d >= 1 << 31 else d
            if d <= 0 or d > self.N * self.fs:  # a pause longer than the history: start over
                d, self.t, self.e = None, [], []
        t = 0.0 if d is None else self.t_last + d / self.fs
        self.t_last, self.end_last = t, end
        self.t = (self.t + [t - 0.5])[-self.N:]
        self.e = (self.e + [est])[-self.N:]
        raw = -predict(self.t, self.e, t - 0.5 + self.HORIZON)
        target = min(max(raw, self.target - self.MAX_STEP), self.target + self.MAX_STEP)
        target = min(max(target, -self.MAX), self.MAX)
        if abs(target) < self.GAP:
            target = self.GAP if raw > self.GAP / 2 else -self.GAP if raw < -self.GAP / 2 else 0.0
        a = 32 if target >= 0 else 31
        x = 32 * target * 1e-9 if target >= 0 else 1.0 + 32 * target * 1e-9
        f = Fraction(x).limit_denominator(self.MAX_C)
        b, c = f.numerator, f.denominator
        if b == c:
            a, b, c = a + 1, 0, 1
        got = ((a - 32) * c + b) / (32 * c) * 1e9
        if abs(got - self.target) < self.MIN_CHANGE:
            return None
        self.target = got
        return a, b, c


def combine(ppb):
    """Per window (rows), the common offset over satellites (columns, nan unusable), by Combiner;
    also how many were used."""
    est, used = np.full(len(ppb), np.nan), np.zeros(len(ppb), dtype=int)
    c = Combiner()
    for j, row in enumerate(ppb):
        e, u = c({k: v for k, v in enumerate(row) if not np.isnan(v)})
        if u:
            est[j], used[j] = e, len(u)
    return est, used


def allan(y, taus):
    """Overlapping Allan deviation of fractional frequencies y sampled each second."""
    out = []
    for m in taus:
        if 2 * m >= len(y):
            break
        c = np.cumsum(np.r_[0.0, y])
        avg = (c[m:] - c[:-m]) / m
        d = avg[m:] - avg[:-m]
        out.append((m, np.sqrt(0.5 * np.mean(d * d))))
    return out


def survey(z, sats, fs, step):
    """The mean of position fixes at each step of samples, from the satellites' code alone."""
    n_end = min(int(z[f"{s['prn']}_n0"][-2]) for s in sats)
    n_start = max(int(z[f"{s['prn']}_n0"][2500]) for s in sats)
    pos = []
    for n_r in range(n_start, n_end, step):
        t_tx = [fix.transmit_time(z, s["prn"], s["subframes"][0], n_r) for s in sats]
        p, _, _, _ = nav.solve(t_tx, [s["eph"] for s in sats], max(t_tx) + 0.075)
        pos.append(p)
    return np.mean(pos, axis=0), len(pos)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dumps", help="track.py's .npz, or with --channels a hackrf_gps stream")
    ap.add_argument("--channels", help="ch:prn,... : read the file as a hackrf_gps stream")
    ap.add_argument("--lo-offset", type=float, default=0.0, metavar="HZ",
                    help="the tuner's offset from 1575.42 MHz, as hackrf_gps printed it")
    ap.add_argument("--hold", help="ECEF x,y,z (m) of the held position; else the mean of fixes")
    ap.add_argument("--window-s", type=float, default=1.0)
    ap.add_argument("--skip-s", type=float, default=0.0, help="leave out the first seconds")
    ap.add_argument("--csv", help="write each window's estimate and per-satellite offsets here")
    args = ap.parse_args()
    z = fix.load_stream(args.dumps, args.channels) if args.channels else np.load(args.dumps)
    fs = float(z["fs"])
    prns = sorted({int(k.split("_")[0]) for k in z if k.endswith("_IP")})
    sats = [s for s in (fix.satellite(z, p) for p in prns if p < 120) if s["eph"] is not None]
    print("satellites with ephemeris:", " ".join(str(s["prn"]) for s in sats))
    if not sats:
        return
    if args.hold:
        rx = np.array([float(v) for v in args.hold.split(",")])
    else:
        if len(sats) < 4:
            print("a survey needs four satellites; give --hold")
            return
        rx, n = survey(z, sats, fs, int(fs))
        print(f"held position: the mean of {n} fixes")
    start = max(int(z[f"{s['prn']}_n0"][0]) for s in sats) + int(args.skip_s * fs)
    end = min(int(z[f"{s['prn']}_n0"][-1]) for s in sats)
    boundaries = np.arange(start, end, int(args.window_s * fs))
    per = [satellite_offsets(z, s, boundaries, rx, args.lo_offset) for s in sats]
    ppb = np.stack([p[:, 0] for p in per], axis=1)
    est, used = combine(ppb)
    ok = ~np.isnan(est)
    t = (boundaries[:-1] - boundaries[0]) / fs
    print(f"{ok.sum()} of {len(est)} windows of {args.window_s:g} s with an estimate")
    for k, s in enumerate(sats):
        v = ppb[:, k] - est
        g = ~np.isnan(v)
        locks = per[k][:, 2]
        print(f"  PRN {s['prn']}: {g.sum()} windows used, residual {np.std(v[g]) if g.any() else 0:.3f} "
              f"ppb rms, phase lock {np.nanmedian(locks):.2f}")
    if ok.sum() < 2:
        return
    print(f"offset {est[ok][0]:.2f} ppb at {t[ok][0]:.0f} s, {est[ok][-1]:.2f} ppb at {t[ok][-1]:.0f} s; "
          f"satellites used per window {used[ok].min()}-{used[ok].max()}")
    spread = np.nanstd(ppb - est[:, None], axis=1)
    print(f"satellite scatter about the estimate: {np.nanmedian(spread[ok & (used > 1)]):.3f} ppb rms "
          "(median over windows)")
    if ok.all() and args.window_s == 1.0:
        print("Allan deviation: " + ", ".join(f"{m} s {a * 1e-9:.2e}" for m, a in
                                              allan(est, (1, 2, 5, 10, 20, 50, 100))))
    if args.csv:
        np.savetxt(args.csv, np.column_stack([t, est, used, ppb]), fmt="%.6f",
                   header="t_s estimate_ppb used " + " ".join(f"prn{s['prn']}" for s in sats))


if __name__ == "__main__":
    main()
