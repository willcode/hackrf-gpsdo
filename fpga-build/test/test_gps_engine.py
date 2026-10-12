#
# SPDX-License-Identifier: BSD-3-Clause

"""Simulation of the time-shared correlator (dsp/gps_engine.py) against gps_track_model.Channel.

Channels are armed on different PRNs with different code and carrier words and staggered starts.
Some receive new words after a dump, as a loop would send them. A record holds the finished period's
sums and the next period's start; so record j must equal the model's dump j in its sums, and the
model's dump j + 1 in its start sample, code fraction and words.
"""

import unittest

import numpy as np
from amaranth.sim import Simulator

from dsp import gps_model as g, gps_track_model as tm
from dsp.gps_engine import GPSEngine, MAGIC, RECORD_WORDS, SNAP_MAGIC, decode_record, decode_snapshot

MASK32 = (1 << 32) - 1


def run(chans, x, gap=8, updates=None, out_ready=lambda cyc: True, snap_at=None, snap_value=0,
        before=(), stalls=None):
    """chans: list of dicts prn, start, W, Wc. updates: {(ch, record number): (W, Wc)} sent after that
    record arrives. snap_at: a cycle at which to request a snapshot of (snap_value + 1) * 1,024
    samples, which the records are read around. before: requests (kind, ch, W, Wc) sent ahead of the
    arms. stalls: {sample index: extra idle cycles before that sample}."""
    dut = GPSEngine()
    sim = Simulator(dut)
    sim.add_clock(1 / 32e6)
    words, recs, seen, raw = [], [], {}, []
    updates = dict(updates or {})
    flags = {}

    async def bench(ctx):
        queue = [dict(kind=k, ch=c, W=W, Wc=Wc, g1=0, g2=0, win=0, start=0) for k, c, W, Wc in before]
        for c, chan in enumerate(chans):
            g1, g2, win = tm.arm_values(chan["prn"])
            queue.append(dict(kind=0, ch=c, W=chan["W"], Wc=chan["Wc"] & MASK32, g1=g1, g2=g2,
                              win=win, start=chan["start"] - 1))   # the sample before the first
        n, next_at, cyc, taken = 0, 0, 0, False
        while n < len(x):
            # A request is taken at the clock edge where cfg_ready is high, so valid stays up through
            # that edge and drops after it.
            if taken:
                ctx.set(dut.cfg_valid, 0)
                taken = False
            elif queue and not ctx.get(dut.cfg_valid):
                req = queue.pop(0)
                ctx.set(dut.cfg_valid, 1)
                for k in ("kind", "ch", "W", "Wc", "g1", "g2", "win", "start"):
                    ctx.set(getattr(dut, f"cfg_{k}"), req[k])
            if cyc == next_at:
                ctx.set(dut.in_valid, 1)
                ctx.set(dut.in_i, int(x[n].real))
                ctx.set(dut.in_q, int(x[n].imag))
                n += 1
                next_at = cyc + gap + (stalls or {}).get(n, 0)
            else:
                ctx.set(dut.in_valid, 0)
            rdy = out_ready(cyc)
            ctx.set(dut.out_ready, rdy)
            ctx.set(dut.snap_valid, cyc == snap_at)
            ctx.set(dut.snap_value, snap_value)
            if rdy and ctx.get(dut.out_valid):
                words.append(ctx.get(dut.out_data))
                if snap_at is not None:
                    raw.append(words.pop())        # read after the run, around the snapshot
                if len(words) == RECORD_WORDS:
                    rec = decode_record(words)
                    words.clear()
                    recs.append(rec)
                    k = seen.get(rec["ch"], 0)
                    seen[rec["ch"]] = k + 1
                    if (rec["ch"], k) in updates:
                        W, Wc = updates.pop((rec["ch"], k))
                        queue.append(dict(kind=1, ch=rec["ch"], W=W, Wc=Wc & MASK32, g1=0, g2=0,
                                          win=0, start=0))
            if ctx.get(dut.cfg_valid) and ctx.get(dut.cfg_ready):
                taken = True
            flags["held"] = bool(queue) or bool(ctx.get(dut.cfg_valid))
            flags["overrun"] = ctx.get(dut.overrun)
            flags["dump_lost"] = ctx.get(dut.dump_lost)
            cyc += 1
            await ctx.tick()

    sim.add_testbench(bench)
    sim.run()
    i = 0
    while i < len(raw):
        k = (snap_value + 1) * 1024
        if i + k + 3 <= len(raw) and raw[i + k] == SNAP_MAGIC | (snap_value << 8):
            flags["snapshot"] = decode_snapshot(raw[i:i + k + 3])
            i += k + 3
        else:
            recs.append(decode_record(raw[i:i + RECORD_WORDS]))
            i += RECORD_WORDS
    return recs, flags


def model_dumps(chan, x, fs, updates_for_ch):
    ch = tm.Channel(chan["prn"], fs, chan["start"], g.CHIP_RATE, 0.0)
    ch.W, ch.Wc = chan["W"], chan["Wc"]
    dumps, j = [], 0
    while True:
        d = ch.period(x)
        if d is None:
            break
        if j in updates_for_ch:
            ch.pending = updates_for_ch[j]
        dumps.append(d)
        j += 1
    return dumps


class TestEngine(unittest.TestCase):
    def check(self, chans, x, fs, updates, min_records, **kw):
        recs, flags = run(chans, x, updates=updates, **kw)
        self.assertFalse(flags["held"], "a request never taken")
        self.assertFalse(flags["overrun"], "input overrun")
        self.assertFalse(flags["dump_lost"], "dump lost")
        for c, chan in enumerate(chans):
            want = model_dumps(chan, x, fs, {j: v for (cc, j), v in updates.items() if cc == c})
            got = [r for r in recs if r["ch"] == c]
            self.assertGreaterEqual(len(got), min_records, f"channel {c}: too few records")
            for j, rec in enumerate(got):
                if j + 1 >= len(want):
                    break
                for k in ("IE", "QE", "IP", "QP", "IL", "QL"):
                    self.assertEqual(rec[k], want[j][k], f"channel {c} record {j} {k}")
                nxt = want[j + 1]
                self.assertEqual(rec["n_start"], nxt["n0"] & MASK32, f"channel {c} record {j} start")
                self.assertEqual(rec["F_start"], nxt["F"], f"channel {c} record {j} fraction")
                self.assertEqual(rec["W"], nxt["W"], f"channel {c} record {j} W")
                self.assertEqual(rec["Wc"], nxt["Wc"], f"channel {c} record {j} Wc")
        return recs

    def test_three_channels_with_new_words(self):
        fs = 4e6
        x = g.synthesize(fs, 9000, [dict(prn=5, code_phase=100.3, doppler=4000.0, cn0=70)],
                         seed=11, noise_rms=20)
        # About 2.2 samples a chip, so a period is about 2,270 samples.
        W0 = int(0.45 * (1 << 32)) + 12345
        chans = [
            dict(prn=5, start=37, W=W0, Wc=tm.carrier_word(fs, 4000.0)),
            dict(prn=131, start=500, W=W0 + 999_999, Wc=tm.carrier_word(fs, -2500.0)),
            dict(prn=17, start=1203, W=W0 - 777_777, Wc=tm.carrier_word(fs, 123.4)),
        ]
        updates = {(0, 0): (W0 + 5000, tm.carrier_word(fs, 4100.0)),
                   (2, 0): (W0 - 700_000, tm.carrier_word(fs, -50.0))}
        self.check(chans, x, fs, updates, min_records=2)

    def test_seven_channels_shared_boundaries_and_slow_output(self):
        # Starts one sample apart put every boundary in the same few samples, so dumps queue behind
        # the readout; the output also stalls a third of the time.
        fs = 4e6
        x = g.synthesize(fs, 7400, [dict(prn=3, code_phase=5.0, doppler=-1500.0, cn0=65)],
                         seed=12, noise_rms=20)
        W0 = int(0.45 * (1 << 32))
        chans = [dict(prn=p, start=100 + c, W=W0 + 1000 * c, Wc=tm.carrier_word(fs, 300.0 * c - 900))
                 for c, p in enumerate((3, 7, 12, 19, 24, 131, 135))]
        updates = {(c, 0): (W0 + 1000 * c - 4000, tm.carrier_word(fs, 300.0 * c - 880)) for c in range(7)}
        self.check(chans, x, fs, updates, min_records=2, out_ready=lambda cyc: cyc % 3 != 0)

    def test_requests_not_applied_and_a_start_after_an_input_gap(self):
        # Pending requests to a channel never armed are taken and discarded, so the second does not
        # wait for the first to apply and hold the port; and channel 0 starts on time just after an
        # input gap, while the slot rests at 0 between samples.
        fs = 4e6
        x = g.synthesize(fs, 5000, [dict(prn=9, code_phase=40.0, doppler=-700.0, cn0=70)],
                         seed=14, noise_rms=20)
        W0 = int(0.45 * (1 << 32)) + 4321
        chans = [dict(prn=9, start=400, W=W0, Wc=tm.carrier_word(fs, -700.0)),
                 dict(prn=22, start=430, W=W0 + 20_000, Wc=tm.carrier_word(fs, 900.0))]
        before = [(1, 3, W0, 0), (1, 3, W0 + 1, 0), (1, 1, W0 + 2, 0)]
        updates = {(0, 0): (W0 + 3000, tm.carrier_word(fs, -650.0))}
        self.check(chans, x, fs, updates, min_records=2, before=before, stalls={399: 5})

    def test_record_waiting_out_a_snapshot(self):
        # A period of about 4,000 samples and a snapshot of 4,099 words: the record of the boundary
        # inside the snapshot waits for its end, which falls just before the next boundary or on
        # it. Each record sent must be whole: its sums the period before its start, as the model's.
        fs = 4e6
        x = g.synthesize(fs, 12300, [dict(prn=5, code_phase=10.0, doppler=500.0, cn0=70)],
                         seed=13, noise_rms=20)
        chan = dict(prn=5, start=60, W=int(1023 / 4000.3 * (1 << 32)), Wc=tm.carrier_word(fs, 500.0))
        want = model_dumps(chan, x, fs, {})
        sums = {want[j + 1]["n0"]: want[j] for j in range(len(want) - 1)}
        b2 = want[2]["n0"]
        for delta in (-10, 0):                         # the snapshot's end against the boundary
            recs, flags = run([chan], x[:b2 + 300], snap_at=8 * (b2 - 4099 + delta) + 3, snap_value=3)
            self.assertIn("snapshot", flags)
            starts = [r["n_start"] for r in recs]
            self.assertEqual(starts, sorted(set(starts)), f"{delta}: a start repeated")
            self.assertIn(b2, starts, f"{delta}")
            for r in recs:
                for k in ("IE", "QE", "IP", "QP", "IL", "QL"):
                    self.assertEqual(r[k], sums[r["n_start"]][k], f"{delta}: record at {r['n_start']} {k}")


class TestSnapshot(unittest.TestCase):
    def test_snapshot_samples_and_index(self):
        # A snapshot requested at an arbitrary cycle: its words are the input samples from the
        # index its trailer gives, and the stream stalls part of the time.
        rng = np.random.default_rng(5)
        x = rng.integers(-128, 128, 3600) + 1j * rng.integers(-128, 128, 3600)
        dut = GPSEngine()
        sim = Simulator(dut)
        sim.add_clock(1 / 32e6)
        words = []

        async def bench(ctx):
            n, cyc = 0, 0
            while n < len(x):
                # A second request, for a longer one, in the middle of the first snapshot.
                ctx.set(dut.snap_valid, cyc in (8 * 333 + 3, 8 * 800 + 5))
                ctx.set(dut.snap_value, int(cyc > 8 * 333 + 3))
                if cyc % 8 == 0:
                    ctx.set(dut.in_valid, 1)
                    ctx.set(dut.in_i, int(x[n].real))
                    ctx.set(dut.in_q, int(x[n].imag))
                    n += 1
                else:
                    ctx.set(dut.in_valid, 0)
                rdy = cyc % 5 != 0
                ctx.set(dut.out_ready, rdy)
                if rdy and ctx.get(dut.out_valid):
                    words.append(ctx.get(dut.out_data))
                cyc += 1
                await ctx.tick()
            self.assertFalse(ctx.get(dut.snap_lost))

        sim.add_testbench(bench)
        sim.run()
        # Two snapshots, one after the other, each with its own trailer and length.
        self.assertEqual(len(words), (1024 + 3) + (2048 + 3))
        first, samples = decode_snapshot(words[:1027])
        self.assertGreater(first, 333)
        self.assertLess(first, 340)
        np.testing.assert_array_equal(np.array(samples), x[first:first + 1024])
        second, samples2 = decode_snapshot(words[1027:])
        self.assertGreaterEqual(second, first + 1024 + 3)
        np.testing.assert_array_equal(np.array(samples2), x[second:second + 2048])


if __name__ == "__main__":
    unittest.main()
