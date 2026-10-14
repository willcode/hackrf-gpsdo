#!/usr/bin/env python3
#
# SPDX-License-Identifier: BSD-3-Clause

"""Position, receiver clock bias and clock rate from the dumps gps/track.py saved, or from a stream
hackrf_gps saved (--channels names each channel's PRN, as hackrf_gps printed them).

Per satellite: bit timing, every subframe whose ten words pass parity, and the ephemeris from
subframes 1-3. A subframe's start fixes the GPS time of one code period; each dump's code NCO state then
gives the transmit time at any sample. Solutions at one-second steps give the clock bias against the
sample count, whose slope is the sample clock's (and so the TCXO's) frequency offset.
"""

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import nav  # noqa: E402
import stream  # noqa: E402

CHIP_RATE = 1.023e6
ONE = float(1 << 32)


def load_stream(path, channels):
    """A hackrf_gps stream as track.py's dumps, for the channels in "ch:prn,...". A record carries the
    start of the period it opens and the sums of the one it closes, so its sums go with the previous
    record's start."""
    _, recs = stream.split(np.fromfile(path, dtype="<u2"))
    z = {"fs": 4e6}
    for item in channels.split(","):
        ch, prn = (int(v) for v in item.split(":"))
        rs = stream.fill_gaps(recs[ch])  # missing records back as placeholders, indices on time
        n = np.array([r["n_start"] for r in rs], dtype=np.int64)
        n = n[0] + np.concatenate(([0], np.cumsum(np.diff(n) % (1 << 32))))
        z[f"{prn}_n0"] = n[:-1]
        for key, name in (("F", "F_start"), ("W", "W"), ("Wc", "Wc")):
            z[f"{prn}_{key}"] = np.array([r[name] for r in rs[:-1]], dtype=np.float64)
        for key in ("IP", "QP"):
            z[f"{prn}_{key}"] = np.array([r[key] for r in rs[1:]], dtype=np.float64)
    return z


def satellite(z, prn):
    p = z[f"{prn}_IP"] + 1j * z[f"{prn}_QP"]
    edge, _ = nav.bit_sync(p[1000:])
    edge = (edge + 1000) % nav.BIT_MS
    best = None
    for coherent in (True, False):
        soft = nav.bits(p.astype(complex), edge, coherent=coherent)
        hard = [int(v < 0) for v in soft]
        if not coherent:
            hard = list(np.cumsum([0] + hard) % 2)
        subs = []
        for inv in (0, 1):
            b = [v ^ inv for v in hard]
            for i in range(2, len(b) - 300 + 1):
                if tuple(b[i:i + 8]) != nav.PREAMBLE:
                    continue
                ws = nav.words(b, i)
                if ws is None:
                    continue
                sf = nav.decode_subframe(ws)
                if 1 <= sf["id"] <= 5:
                    sf["dump"] = edge + i * nav.BIT_MS
                    sf["bit"] = i
                    subs.append(sf)
        subs = nav.consistent(sorted(subs, key=lambda sf: sf["bit"]))
        if best is None or len(subs) > len(best):
            best = subs
    eph = {}
    for sf in best:
        if sf["id"] in (1, 2, 3):
            eph.setdefault(sf["id"], sf)
    ok = len(eph) == 3 and eph[2]["iode"] == eph[3]["iode"] == eph[1]["iodc"] & 0xFF
    merged = {**eph.get(1, {}), **eph.get(2, {}), **eph.get(3, {})} if ok else None
    if merged is not None and not nav.plausible(merged):
        merged = None   # parity passes for anything sent: an orbit no satellite has
    return dict(prn=prn, subframes=best, eph=merged)


def smoothing(z, prn, smooth_max, lo_offset=0.0):
    """Carrier smoothing's correction to each record's code phase, in chips, for offline work on
    unsteered runs (the GPS app smooths nothing: PLL B's steering is in the recorded carrier words,
    which smoothing would read as code-carrier divergence): div sums the code's
    advance less the carrier's (the nominal code word's samples plus the carrier cycles / 1540),
    smooth follows it with weight 1/M, M rising to smooth_max, and a step over 0.5 chip starts it
    over. lo_offset is the tuner's offset from 1575.42 MHz in hertz (hackrf_gps prints it), which is
    in the carrier words and not in the code. Needs the stream's carrier words (load_stream)."""
    n0, F, Wc = z[f"{prn}_n0"], z[f"{prn}_F"], z[f"{prn}_Wc"]
    w_nom = float(round(CHIP_RATE / float(z["fs"]) * ONE))
    lo_word = round(lo_offset * 1e6) * 1e-6 / float(z["fs"]) * ONE
    delta = np.zeros(len(n0))
    div = smooth = 0.0
    count = 0
    for k in range(1, len(n0)):
        dn = float(n0[k] - n0[k - 1])
        code = 1023.0 + (F[k] - F[k - 1]) / ONE
        carrier = (dn * w_nom + dn * (Wc[k - 1] + lo_word) / 1540.0) / ONE
        div += code - carrier
        step = div - smooth
        if count == 0 or step > 0.5 or step < -0.5:
            smooth, count = div, 1
        else:
            if count < smooth_max:
                count += 1
            smooth += step / count
        delta[k] = smooth - div
    return delta


def transmit_time(z, prn, anchor, n_r, delta=None):
    """GPS time of week at which the code received at sample n_r left the satellite; with delta
    (smoothing()), from the carrier-smoothed code phase."""
    n0, F, W = z[f"{prn}_n0"], z[f"{prn}_F"], z[f"{prn}_W"]
    j = int(np.searchsorted(n0, n_r, side="right") - 1)
    t_sf = anchor["tow"] * 6.0 - 6.0
    d = 0.0 if delta is None else delta[j]
    t = t_sf + (j - anchor["dump"]) * 1e-3 + (F[j] / ONE + d) / CHIP_RATE
    return t + (n_r - n0[j]) * (W[j] / ONE) / CHIP_RATE


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("dumps", help="track.py's .npz, or with --channels a hackrf_gps stream")
    ap.add_argument("--channels", help="ch:prn,... : read the file as a hackrf_gps stream")
    ap.add_argument("--smooth", type=int, default=0, metavar="M",
                    help="carrier smoothing over at most M records, for unsteered runs "
                         "(the GPS app smooths nothing; a stream only)")
    ap.add_argument("--lo-offset", type=float, default=0.0, metavar="HZ",
                    help="the tuner's offset from 1575.42 MHz, as hackrf_gps printed it")
    ap.add_argument("--step-s", type=float, default=1.0)
    ap.add_argument("--hold", help="ECEF x,y,z (m) of a surveyed position: solve the clock only, from "
                                   "any number of satellites")
    args = ap.parse_args()
    z = load_stream(args.dumps, args.channels) if args.channels else np.load(args.dumps)
    fs = float(z["fs"])
    prns = sorted({int(k.split("_")[0]) for k in z if k.endswith("_IP")})
    sats = []
    for prn in prns:
        s = satellite(z, prn)
        ids = [int(sf["id"]) for sf in s["subframes"]]
        tows = sorted({int(sf["tow"]) for sf in s["subframes"]})
        print(f"PRN {prn}: subframes {ids}, TOW counts {tows}, "
              f"ephemeris {'complete' if s['eph'] else 'incomplete'}")
        if s["eph"]:
            e = s["eph"]
            print(f"   week {e['week']} (mod 1024), health {e['health']}, IODE {e['iode']}, "
                  f"sqrt(A) {e['sqrta']:.3f}, e {e['e']:.5f}, i0 {np.degrees(e['i0']):.2f} deg")
            sats.append(s)
    hold = np.array([float(v) for v in args.hold.split(",")]) if args.hold else None
    deltas = ({s["prn"]: smoothing(z, s["prn"], args.smooth, args.lo_offset) for s in sats}
              if args.smooth else {})
    if len(sats) < (1 if hold is not None else 4):
        print(f"{len(sats)} satellites with ephemeris; a fix needs 4, a held position 1")
        return
    n_end = min(int(z[f"{s['prn']}_n0"][-2]) for s in sats)
    n_start = max(int(z[f"{s['prn']}_n0"][2500]) for s in sats)
    times = np.arange(n_start, n_end, int(args.step_s * fs))
    t0 = None
    rows, kept = [], []
    for n_r in times if hold is None else []:
        t_tx = [transmit_time(z, s["prn"], s["subframes"][0], n_r, deltas.get(s["prn"])) for s in sats]
        if t0 is None:
            t0 = max(t_tx) + 0.075 - n_r / fs
        t_rx = t0 + n_r / fs
        pos, bias, res, els = nav.solve(t_tx, [s["eph"] for s in sats], t_rx)
        if np.all(np.isfinite(pos)):
            rows.append((n_r / fs, pos, bias, res, els))
            kept.append(n_r)
    if hold is not None:
        t = times / fs
        held = []
        for n_r in times:
            t_tx = [transmit_time(z, s["prn"], s["subframes"][0], n_r, deltas.get(s["prn"])) for s in sats]
            if t0 is None:
                t0 = max(t_tx) + 0.075 - n_r / fs
            b, _ = nav.solve_clock(t_tx, [s["eph"] for s in sats], t0 + n_r / fs, hold)
            held.append(b)
        report_held(t, np.array(held))
        return
    if not rows:
        print(f"no solution from {len(sats)} satellites")
        return
    pos = np.array([r[1] for r in rows])
    bias = np.array([r[2] for r in rows])
    t = np.array([r[0] for r in rows])
    mean = pos.mean(axis=0)
    lat, lon, h = nav.ecef_to_geodetic(mean)
    spread = np.sqrt(((pos - mean) ** 2).sum(axis=1)).max()
    print(f"\n{len(rows)} solutions from {len(sats)} satellites over {t[-1] - t[0]:.0f} s")
    print(f"position {lat:.5f}, {lon:.5f}, height {h:.0f} m; largest distance from the mean {spread:.1f} m")
    print(f"  ECEF for --hold: {mean[0]:.1f},{mean[1]:.1f},{mean[2]:.1f}")
    if rows[0][4] is not None:   # none when the solve stayed off the Earth
        print("elevations (deg):", " ".join(f"{s['prn']}:{e:.0f}" for s, e in zip(sats, rows[0][4])))
    print("residuals at the first solution (m):", " ".join(f"{r:.1f}" for r in rows[0][3]))
    t_tx = [transmit_time(z, s["prn"], s["subframes"][0], int(kept[0]), deltas.get(s["prn"]))
            for s in sats]
    at = [nav.sat_position(s["eph"], t)[0] for s, t in zip(sats, t_tx)]
    print("DOP at the first solution: G {:.1f}, P {:.1f}, H {:.1f}, V {:.1f}, T {:.1f}".format(
        *nav.dop(rows[0][1], at)))
    report_clock("full solution", t, bias)
    held = []
    for n_r in kept:
        t_tx = [transmit_time(z, s["prn"], s["subframes"][0], n_r, deltas.get(s["prn"])) for s in sats]
        b, _ = nav.solve_clock(t_tx, [s["eph"] for s in sats], t0 + n_r / fs, mean)
        held.append(b)
    report_held(t, np.array(held))


def report_held(t, held):
    report_clock("position held", t, held)
    span = t[-1] - t[0]
    if span >= 60:
        print("oscillator offset (ppb) in 20 s windows, position held:")
        out = []
        for a in np.arange(t[0], t[-1] - 20 + 1e-9, 20):
            m = (t >= a) & (t < a + 20)
            if m.sum() >= 5:
                out.append(f"{a - t[0]:.0f}s:{np.polyfit(t[m], held[m], 1)[0] * 1e9:.1f}")
        print("  " + "  ".join(out))


def report_clock(label, t, bias):
    slope, icpt = np.polyfit(t, bias, 1)
    fit = bias - (slope * t + icpt)
    print(f"clock, {label}: rate {slope * 1e9:.2f} ns/s = {slope * 1e6:.4f} ppm (the sample clock "
          f"runs {'fast' if slope > 0 else 'slow'}); rms about the line {fit.std() * 1e9:.1f} ns")


if __name__ == "__main__":
    main()
