#!/usr/bin/env python3
#
# SPDX-License-Identifier: BSD-3-Clause

"""Read a stream that hackrf_gps saved (records and snapshots of the GPS image) and report each
channel's tracking: records, gaps, C/N0, phase lock on symbol sums, and the navigation subframes
the prompt sums carry.

Records (21 words) and snapshots (fpga-build/dsp/gps_engine.py) are read in order.
"""

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fpga-build"))
sys.path.insert(0, HERE)

from dsp import gps_engine as ge, gps_track_model as tm  # noqa: E402
import nav  # noqa: E402

KEYS = ("IE", "QE", "IP", "QP", "IL", "QL")


# The snapshot sizes the receiver asks for (setting 1 = v, (v + 1) * 1,024 samples): the app's own
# search takes 4,096 samples, the host's (hackrf_gps -H) 32 ms. A stream may hold both.
SNAP_VALUES = (3, 127)


def split(words, snap_values=SNAP_VALUES):
    """(snapshots as (first index, samples), records by channel as lists of dicts), read in order:
    a record where its magic word stands, a snapshot of one of snap_values where its trailer stands
    one snapshot later. Where both readings fit (a trailer value inside the records, or a snapshot
    whose first sample looks like a record's magic word), the one that something valid follows
    wins, the record when both do."""
    sizes = [((v + 1) * 1024, ge.SNAP_MAGIC | (v << 8)) for v in snap_values]

    def record_at(j):
        w = int(words[j])
        return (w & 0xFF) == ge.MAGIC and (w >> 8) < 7 and j + ge.RECORD_WORDS <= len(words)

    def snapshot_at(j):
        for n, trailer in sizes:
            if j + n + 3 <= len(words) and int(words[j + n]) == trailer:
                return n
        return None

    def starts(j):
        return j + 3 > len(words) or record_at(j) or snapshot_at(j) is not None

    snaps, recs, i = [], {}, 0
    while i + 3 <= len(words):
        n = snapshot_at(i)
        if record_at(i) and (n is None or starts(i + ge.RECORD_WORDS) or not starts(i + n + 3)):
            r = ge.decode_record(words[i:i + ge.RECORD_WORDS])
            recs.setdefault(r["ch"], []).append(r)
            i += ge.RECORD_WORDS
        elif n is not None:
            snaps.append(ge.decode_snapshot(words[i:i + n + 3]))
            i += n + 3
        else:
            i += 1
    return snaps, recs


def fill_gaps(rs):
    """A channel's records with each missing one (a snapshot pauses the stream for a record or
    two) put back as a placeholder: its period start one period on, the words of the record before,
    the sums zero, so record indices stay on time as the app keeps them (gps_loops_skip). A jump of
    more than 6 periods (a misread record) is left as it is. Of two records with the same start, the
    later is kept: before the engine's fix, a record that waited out a snapshot could start on its
    channel's next boundary and come out with that start and a few samples' sums, the true record
    after it."""
    out = []
    for r in rs:
        if out and r["n_start"] == out[-1]["n_start"]:
            out[-1] = r
            continue
        if out:
            d = (r["n_start"] - out[-1]["n_start"]) % (1 << 32)
            missed = (d + 2000) // 4000 - 1
            if 0 < missed <= 5:
                step = d // (missed + 1)
                for k in range(missed):
                    p = dict(out[-1])
                    p["n_start"] = (p["n_start"] + step) & 0xFFFFFFFF
                    for key in KEYS:
                        p[key] = 0
                    out.append(p)
        out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("file")
    ap.add_argument("--snap-value", type=int, action="append",
                    help="a snapshot size to recognize, (v + 1) * 1,024 samples; repeatable "
                         f"(default {', '.join(map(str, SNAP_VALUES))})")
    args = ap.parse_args()
    words = np.fromfile(args.file, dtype="<u2")
    snaps, recs = split(words, args.snap_value or SNAP_VALUES)
    print(f"{len(words)} words: {len(snaps)} snapshots"
          + "".join(f" (first sample {s[0]})" for s in snaps))
    for ch, rs in sorted(recs.items()):
        filled = fill_gaps(rs)
        gaps = len(filled) - len(rs)  # records missing, put back as placeholders
        rs = filled
        dumps = [{k: r[k] for k in KEYS} for r in rs]
        late = dumps[2000:] if len(dumps) > 3000 else dumps
        cn0, _ = tm.cn0_windows(late)
        p = np.array([d["IP"] + 1j * d["QP"] for d in dumps])
        # The symbol edge: where the prompt's sign changes cluster, after the loops settle.
        flips = np.zeros(20)
        ang = np.angle(p[1300:])
        d = np.abs(np.angle(np.exp(1j * np.diff(ang))))
        for k in np.flatnonzero(d > np.pi / 2):
            flips[(1300 + k + 1) % 20] += 1
        edge = int(np.argmax(flips))
        lock = tm.phase_lock(late, 20, (edge - (len(dumps) - len(late))) % 20)
        print(f"channel {ch}: {len(rs)} records ({len(rs) / 1000:.1f} s), {gaps} gaps, "
              f"C/N0 {cn0:.1f} dB-Hz, symbol-sum lock {lock:.2f}, edge {edge}")
        for coherent in (True, False):
            soft = nav.bits(p, edge, coherent=coherent)
            hard = [int(v < 0) for v in soft]
            if not coherent:
                hard = list(np.cumsum([0] + hard) % 2)
            subs = nav.consistent(nav.subframes(hard))
            tows = sorted({int(x["tow"]) for x in subs})
            print(f"  {'coherent' if coherent else 'differential'}: {len(subs)} subframe headers "
                  f"pass parity in {len(hard)} bits" + (f", TOW {tows[0]}-{tows[-1]}" if tows else ""))


if __name__ == "__main__":
    main()
