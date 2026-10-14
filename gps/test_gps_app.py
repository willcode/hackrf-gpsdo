#!/usr/bin/env python3
#
# SPDX-License-Identifier: BSD-3-Clause

"""The GPS app (apps/gps/gps_app.c) on the host, against the models.

The app reads records from a simulated sample buffer and writes register 6 through a simulated
service. Checked: an arm message becomes the engine's arm command with every field; each record
becomes a pending command with the words the integer loops (gps_loop_int.IntLoops) give for it,
only when they change, and only once the records show the channel's previous request in use; a
request the engine never applies is counted lost and sent again; a snapshot in the stream is passed
over; the status reply counts it all. And with device runs in data/: drop2.bin (four satellites for
90 s) replayed through the app gives fixes once a second equal to nav.py's from the same records;
nav4.bin (four satellites in a near-degenerate geometry) gives none kept, every one rejected.
"""

import os
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "fpga-build"))

from dsp import gps_cmd, gps_engine as ge, gps_loop_int as li, gps_model as g_model, \
    gps_track_model as tm  # noqa: E402

sys.path.insert(0, HERE)
import fix, freq, nav  # noqa: E402
from test_gps_nav import PAGES, PRN2, to_bits  # noqa: E402

FS = 4e6
M32 = 0xFFFFFFFF

# The device run and its arms, as hackrf_gps printed them: channel, PRN, Doppler (Hz).
RUN = os.path.join(HERE, "data", "drop2.bin")
REC45 = os.path.join(HERE, "data", "l1-45s.iq")
RUN_ARMS = [(0, 27, -25), (1, 16, -2075), (2, 8, 2675), (3, 4, -2475), (4, 120, 750),
            (5, 128, 1750), (6, 18, -975)]
# The tuner's offset from 1575.42 MHz in these runs (radio register freq_rf, applied bank).
LO_OFFSET = -9.46044921875
# A run whose four satellites (27, 8, 16, 2) stood in a near-degenerate geometry (DOP about 1,400):
# the solve leaves the Earth.
BAD_RUN = os.path.join(HERE, "data", "nav4.bin")
BAD_RUN_ARMS = [(0, 27, -1900), (1, 8, 900), (2, 16, -2600), (3, 2, 3550), (4, 10, -3750),
                (5, 122, 750)]


def decode_commands(data):
    """Register-6 bytes as a list of (op, ch, payload bits as an int); checks the toggle."""
    out, cur, toggle = [], None, 0
    for b in data:
        assert (b >> 7) == toggle ^ 1, "toggle bit did not alternate"
        toggle = b >> 7
        if b & 0x40:
            op, ch = (b >> 4) & 3, b & 15
            cur = [op, ch, 0, 0, gps_cmd.LENGTHS[op]]
            if cur[4] == 0:
                out.append((op, ch, 0))
                cur = None
        else:
            assert cur is not None, "payload byte outside a command"
            cur[2] |= (b & 0x3F) << (6 * cur[3])
            cur[3] += 1
            if cur[3] == cur[4]:
                out.append((cur[0], cur[1], cur[2]))
                cur = None
    return out


def field(bits, lo, hi):
    return (bits >> lo) & ((1 << (hi - lo)) - 1)


def arm_message(ch, prn, doppler_mhz, start):
    g1, g2, win = tm.arm_values(prn)
    W = tm.code_word(FS, 1.023e6 * (1 + doppler_mhz * 1e-3 / 1575.42e6))
    Wc = tm.carrier_word(FS, doppler_mhz * 1e-3) & M32
    return (b"A" + bytes([ch, prn]) + W.to_bytes(4, "little") + Wc.to_bytes(4, "little")
            + g1.to_bytes(2, "little") + g2.to_bytes(2, "little") + bytes([win])
            + start.to_bytes(4, "little") + doppler_mhz.to_bytes(4, "little", signed=True))


def after_snapshots(words, snap_value=127):
    """The word index just past the snapshots that open a hackrf_gps stream, back to back (no channel
    runs before the arms); the records of the armed channels follow."""
    n = (snap_value + 1) * 1024
    trailer = ge.SNAP_MAGIC | (snap_value << 8)
    end = 0
    while end + n < len(words) and words[end + n] == trailer:
        end += n + 3
    return end


def record_words(ch, n0, F, W, Wc, sums):
    w = [0xD5 | (ch << 8)]
    for v in (n0, F, W, Wc) + tuple(sums):
        v &= M32
        w += [v & 0xFFFF, v >> 16]
    return w


class TestGPSApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.exe = os.path.join(cls.tmp, "harness")
        a = os.path.join(ROOT, "apps")
        subprocess.run(["cc", "-std=c11", "-O1", "-Wall", "-Wextra", "-Werror", "-DAPP_HOST",
                        "-I", os.path.join(a, "sdk"), "-I", os.path.join(HERE, "app"),
                        "-o", cls.exe, os.path.join(a, "gps", "gps_app.c"),
                        os.path.join(HERE, "app", "gps_freq.c"),
                        os.path.join(HERE, "app", "gps_acq.c"),
                        os.path.join(a, "gps", "host_harness.c"),
                        os.path.join(HERE, "app", "gps_loops.c"),
                        os.path.join(HERE, "app", "gps_nav.c"), "-lm"], check=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def run_app(self, lines):
        out = subprocess.run([self.exe], input="\n".join(lines) + "\n", capture_output=True,
                             text=True, check=True).stdout.splitlines()
        writes = [int(l[2:], 16) for l in out if l.startswith("6 ")]
        replies = [bytes.fromhex(l[2:]) for l in out if l.startswith("r ")]
        self.pll_b = [tuple(int(v) for v in l[2:].split()) for l in out if l.startswith("b ")]
        return writes, replies

    def test_arm_track_snapshot_status(self):
        prn, ch, doppler_mhz = 7, 3, 1234567
        g1, g2, win = tm.arm_values(prn)
        W = tm.code_word(FS, 1.023e6 * (1 + 1.234567 / 1540))
        Wc = tm.carrier_word(FS, 1234.567) & M32
        start = 123456789
        arm = (b"A" + bytes([ch, prn]) + W.to_bytes(4, "little") + Wc.to_bytes(4, "little")
               + g1.to_bytes(2, "little") + g2.to_bytes(2, "little") + bytes([win])
               + start.to_bytes(4, "little") + doppler_mhz.to_bytes(4, "little", signed=True))
        lines = ["m " + arm.hex(), "m " + b"M\x01".hex()]

        # Records with rotating prompt sums and data bits, as the integer loops' test uses.
        rng = np.random.default_rng(7)
        n, sums_seq = 1600, []
        bits = rng.choice([-1, 1], size=n // 20 + 1).repeat(20)[:n]
        for k in range(n):
            p = 2e5 * bits[k] * np.exp(1j * (2 * np.pi * 2.5 * k * 1e-3 + 0.3))
            row = []
            for scale in (0.6, 1.0, 0.55):
                z = scale * p + 6e4 * (rng.standard_normal() + 1j * rng.standard_normal())
                row += [int(round(z.real)), int(round(z.imag))]
            sums_seq.append(row)
        # The app's rule, run here to give each record the words in use: a request sent after record
        # k shows in record k + 1, except one the engine is made to lose (sent after record LOSE),
        # which never shows; after three records without it the app counts it lost and sends its
        # latest words again.
        LOSE = 700
        loops = li.IntLoops(FS, doppler_mhz * 1e-3, symbol_ms=20, fll_ms=300, sync_ms=1000)
        in_use, cur = (W, Wc), (W, tm.carrier_word(FS, 1234.567))
        send, outstanding, waited, sent, lost = False, None, 0, [], 0
        record_in_use = []
        for k, sums in enumerate(sums_seq):
            record_in_use.append(in_use)
            if outstanding is not None:
                if in_use == outstanding:
                    outstanding = None
                else:
                    waited += 1
                    if waited >= 3:
                        outstanding, lost, send = None, lost + 1, True
            words = loops.update(dict(zip(("IE", "QE", "IP", "QP", "IL", "QL"), sums)))
            if words != cur:
                cur, send = words, True
            if send and outstanding is None:
                w = (cur[0] & M32, cur[1] & M32)
                sent.append(w)
                send, outstanding, waited = False, w, 0
                if k != LOSE:
                    in_use = w
        self.assertEqual(lost, 1)

        for k, sums in enumerate(sums_seq):
            words = record_words(ch, start + 1 + 4000 * k, 0, *record_in_use[k], sums)
            if k == 500:
                lines.append("m " + b"S\x00".hex())      # a snapshot, then its 1,027 words
                snap = list(rng.integers(0, 65536, 1024)) + [0x00D7, 0x1234, 0x0001]
                snap[0] = 0x00D5                            # first word passes for a record start
                lines.append("w " + " ".join(f"{v:x}" for v in snap))
            lines.append("w " + " ".join(f"{v:x}" for v in words))
        lines.append("m " + b"Q".hex())
        lines.append("m " + b"L".hex())
        writes, replies = self.run_app(lines)
        cmds = decode_commands(writes)

        # The arm, field by field (dsp/gps_cmd.py's payload layout).
        op, c, b = cmds[0]
        self.assertEqual((op, c), (0, ch))
        self.assertEqual([field(b, 0, 32), field(b, 32, 64), field(b, 64, 74), field(b, 74, 84),
                          field(b, 84, 87), field(b, 87, 119)], [W, Wc, g1, g2, win, start])
        # Settings: 0 = 1 (records), then 1 = 0 (the snapshot, asked for at record 500).
        self.assertEqual([x for x in cmds if x[0] == 3], [(3, 0, 1 << 4), (3, 0, 1 | (0 << 4))])

        # Pending commands: the loops' words, by the rule above.
        got = [(field(b, 0, 32), field(b, 32, 64)) for op, c, b in cmds if op == 1 and c == ch]
        self.assertEqual(got, sent)
        self.assertEqual(len(cmds), 3 + len(sent))

        # Status: the channel's records, the snapshot, no garbage or refused writes.
        # The load reply: ticks, window, longest gap and stack depth; the harness's ticks carry no cycle
        # count, so no window has completed.
        self.assertEqual(replies[-1], bytes(16))
        r = replies[-2]
        self.assertEqual(len(r), 128)
        self.assertEqual(r[16 * ch], prn)
        self.assertEqual(int.from_bytes(r[16 * ch + 4:16 * ch + 8], "little"), n)
        self.assertEqual(r[16 * ch + 2], loops.edge)
        self.assertEqual(r[16 * ch + 3], 1)                             # the lost request
        tail = r[112:]
        self.assertEqual(int.from_bytes(tail[0:4], "little"), n)        # records
        self.assertEqual(int.from_bytes(tail[4:8], "little"), 0)        # garbage words
        self.assertEqual(int.from_bytes(tail[10:12], "little"), 1)      # snapshots
        self.assertEqual(int.from_bytes(tail[12:16], "little"), 0)      # refused writes

    def test_steering_words(self):
        """A PLL B ratio asked for: written at the first tick with no request outstanding, and
        every carrier word sent from that tick on carries the step, the code word unchanged."""
        prn, ch, doppler_mhz = 7, 2, -800000
        W = tm.code_word(FS, 1.023e6 * (1 + doppler_mhz * 1e-3 / 1575.42e6))
        Wc0 = tm.carrier_word(FS, doppler_mhz * 1e-3) & M32
        start = 5000
        q = 8000000
        lines = ["m " + (b"F" + round(LO_OFFSET * 1e6).to_bytes(4, "little", signed=True)
                         + q.to_bytes(4, "little", signed=True)).hex(),
                 "m " + arm_message(ch, prn, doppler_mhz, start).hex()]
        rng = np.random.default_rng(3)
        n, STEP = 900, 400
        bits = rng.choice([-1, 1], size=n // 20 + 1).repeat(20)[:n]
        sums_seq = []
        for k in range(n):
            p = 2e5 * bits[k] * np.exp(1j * (2 * np.pi * 1.5 * k * 1e-3 + 0.2))
            row = []
            for scale in (0.6, 1.0, 0.55):
                z = scale * p + 5e4 * (rng.standard_normal() + 1j * rng.standard_normal())
                row += [int(round(z.real)), int(round(z.imag))]
            sums_seq.append(row)
        lo = 1575.42e6 + round(LO_OFFSET * 1e6) * 1e-6 - q
        delta = 1 / (32 * 1048575)
        w = -lo * delta / FS * 2 ** 32
        step_word = int(w - 0.5) if w < 0 else int(w + 0.5)
        # The app's rule (as in test_arm_track_snapshot_status), with the write: asked for before
        # record STEP's tick, made at the first tick after it with no request outstanding.
        loops = li.IntLoops(FS, doppler_mhz * 1e-3, symbol_ms=20, fll_ms=300, sync_ms=1000)
        in_use, cur = (W, Wc0), (W, tm.carrier_word(FS, doppler_mhz * 1e-3))
        send, outstanding, sent, record_in_use = False, None, [], []
        steer, want, written_at = 0, False, None
        for k, sums in enumerate(sums_seq):
            record_in_use.append(in_use)
            if outstanding is not None and in_use == outstanding:
                outstanding = None
            words = loops.update(dict(zip(("IE", "QE", "IP", "QP", "IL", "QL"), sums)))
            if words != cur:
                cur, send = words, True
            want |= k == STEP
            if want and outstanding is None:
                steer, want, send, written_at = step_word, False, True, k
                sent_before = len(sent)
            if send and outstanding is None:
                w2 = (cur[0] & M32, (cur[1] + steer) & M32)
                sent.append(w2)
                send, outstanding, in_use = False, w2, w2
        self.assertIsNotNone(written_at)
        for k, sums in enumerate(sums_seq):
            if k == STEP:
                lines.append("m " + (b"B" + (32).to_bytes(4, "little") + (1).to_bytes(4, "little")
                                     + (1048575).to_bytes(4, "little")).hex())
            words = record_words(ch, start + 1 + 4000 * k, 0, *record_in_use[k], sums)
            lines.append("w " + " ".join(f"{v:x}" for v in words))
        lines += ["m " + b"J".hex(), "m " + b"R".hex()]
        writes, replies = self.run_app(lines)
        self.assertEqual(self.pll_b, [(32, 1, 1048575)])
        cmds = decode_commands(writes)
        got = [(field(b, 0, 32), field(b, 32, 64)) for op, c, b in cmds if op == 1 and c == ch]
        self.assertEqual(got, sent)
        self.assertGreater(len(sent) - sent_before, 10)                 # words sent with the step
        j, r = replies[-2], replies[-1]
        self.assertEqual(int.from_bytes(j[0:4], "little"), 1)
        e = [int.from_bytes(j[4 + 4 * i:8 + 4 * i], "little") for i in range(7)]
        self.assertEqual(e[0:3], [32, 1, 1048575])
        self.assertEqual(e[5:7], [1000, 1500])                           # the harness's counts
        self.assertEqual(r[59], 0)
        self.assertAlmostEqual(np.frombuffer(r[60:64], "<f4")[0], delta * 1e9, places=4)

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_steering_estimate(self):
        """drop2.bin with a ratio asked for at 45 s, and from 48.5 s every record's carrier word
        carrying the steering word, as the engine's would once the app's steered words reach it.
        The app takes each record's steering word out once a record shows its request's words, so
        every estimate from the cut on equals the unsteered run's (the first record after the cut
        already shows the app's steered words); the windows across the large step (141 Hz) are
        left out."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        q = -8000000
        f9 = (b"F" + round(LO_OFFSET * 1e6).to_bytes(4, "little", signed=True)
              + q.to_bytes(4, "little", signed=True))
        hold = b"U" + rx.astype("<f8").tobytes()
        base = self.replay(RUN, RUN_ARMS, b"R", [f9, hold])
        b = (b"B" + (32).to_bytes(4, "little") + (3).to_bytes(4, "little")
             + (1048575).to_bytes(4, "little"))
        delta = 3 / (32 * 1048575)
        lo = 1575.42e6 + round(LO_OFFSET * 1e6) * 1e-6 - q
        w = -lo * delta / FS * 2 ** 32
        word = int(w - 0.5) if w < 0 else int(w + 0.5)
        words = np.fromfile(RUN, dtype="<u2")
        start = after_snapshots(words)
        n_first = int(words[start + 1]) | (int(words[start + 2]) << 16)
        cut = n_first + int(48.5 * FS)
        i, changed = start, 0
        while i + 21 <= len(words):          # records only, back to back after the snapshots
            n0 = int(words[i + 1]) | (int(words[i + 2]) << 16)
            if ((n0 - cut) & M32) < (1 << 31):
                wc = ((int(words[i + 7]) | (int(words[i + 8]) << 16)) + word) & M32
                words[i + 7], words[i + 8] = wc & 0xFFFF, wc >> 16
                changed += 1
            i += 21
        self.assertGreater(changed, 100000)
        path = os.path.join(self.tmp, "steered.bin")
        words.tofile(path)
        steered = self.replay(path, RUN_ARMS, b"R", [f9, hold], at={45: b})
        self.assertEqual(self.pll_b, [(32, 3, 1048575)])
        self.assertAlmostEqual(np.frombuffer(steered[-1][60:64], "<f4")[0], delta * 1e9, places=3)
        est = lambda rs: {int.from_bytes(r[4:8], "little"): np.frombuffer(r[8:12], "<f4")[0]
                          for r in rs if r[18] == 0}
        e0, e1 = est(base), est(steered)
        after = [k for k in e1 if ((k - cut) & M32) < (1 << 31)]
        self.assertGreater(len(after), 30)
        # Each satellite's offset (the weights then differ a little: the cut's window changed the
        # satellites' variances).
        sat = lambda rs: {int.from_bytes(r[4:8], "little"): (r[16], np.frombuffer(r[24:52], "<f4"))
                          for r in rs if r[18] == 0}
        s0, s1 = sat(base), sat(steered)
        for k in after:
            self.assertEqual(s1[k][0], s0[k][0])
            for c in range(4):
                if s1[k][0] & (1 << c):
                    self.assertAlmostEqual(s1[k][1][c], s0[k][1][c], delta=1e-4)
            self.assertAlmostEqual(e1[k], e0[k], delta=2e-3)
        self.assertEqual(len(set(e0) - set(e1)), 2)                     # the two windows

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_steering_loop(self):
        """drop2.bin with the steering loop on: every ratio the app asks for is the one freq.Servo
        asks for from the same estimates. The recorded carriers never carry the app's steering,
        so the estimates stay the unsteered run's, and the loop slews toward the TCXO's error at
        5 ppb a second."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        f9 = (b"F" + round(LO_OFFSET * 1e6).to_bytes(4, "little", signed=True)
              + (-8000000).to_bytes(4, "little", signed=True))
        replies = self.replay(RUN, RUN_ARMS, b"R", [f9, b"U" + rx.astype("<f8").tobytes()],
                              at={5: b"D\x01"})
        estimates = {int.from_bytes(r[0:4], "little"): r for r in replies if r[18] == 0}
        counts = sorted(estimates)
        self.assertEqual(counts, list(range(1, len(counts) + 1)))
        servo, asks = freq.Servo(), []
        for k in counts:
            r = estimates[k]
            ask = servo(int.from_bytes(r[4:8], "little"), float(np.frombuffer(r[8:12], "<f4")[0]))
            if ask:
                asks.append(ask)
        # The estimates went through float in the reply; the app's own are double. Compare the
        # ratios' steering, which a float's rounding moves by far less than a step.
        self.assertEqual(len(self.pll_b), len(asks))
        self.assertGreater(len(asks), 40)
        ppb = lambda a, b, c: ((a - 32) * c + b) / (32 * c) * 1e9
        for got, want in zip(self.pll_b, asks):
            self.assertAlmostEqual(ppb(*got), ppb(*want), delta=1e-3)
        # The first ask crosses the gap at nominal (29.8 ppb); every later one is within 5 ppb.
        self.assertEqual(self.pll_b[0], (32, 1, 1048575))
        steps = np.diff([ppb(*a) for a in self.pll_b])
        self.assertLessEqual(np.max(np.abs(steps)), 5.001)            # the ratio rounds
        self.assertEqual(replies[-1][59] & 4, 4)                          # the loop on

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_steering_crosses_nominal(self):
        """The tuner's offset ('F') moved so the estimates read about -40 ppb, then about +120 from
        50 s: the steering follows from +40 across the gap at nominal to about -120, as freq.Servo
        does from the same estimates (clamped by the step limit to the gap's edge, it once stayed
        at +29.8 ppb for good)."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        # The run's TCXO reads -126.4 ppb; a tuner offset dHz more reads dHz / L1 lower.
        f = lambda est: (b"F" + round((LO_OFFSET - (est + 126.4) * 1e-9 * 1575.42e6) * 1e6)
                         .to_bytes(4, "little", signed=True) + (-8000000).to_bytes(4, "little", signed=True))
        replies = self.replay(RUN, RUN_ARMS, b"R", [f(-40.0), b"U" + rx.astype("<f8").tobytes()],
                              at={5: b"D\x01", 50: f(120.0)})
        estimates = {int.from_bytes(r[0:4], "little"): r for r in replies if r[18] == 0}
        servo, asks = freq.Servo(), []
        for k in sorted(estimates):
            r = estimates[k]
            ask = servo(int.from_bytes(r[4:8], "little"), float(np.frombuffer(r[8:12], "<f4")[0]))
            if ask:
                asks.append(ask)
        ppb = lambda a, b, c: ((a - 32) * c + b) / (32 * c) * 1e9
        self.assertEqual(len(self.pll_b), len(asks))
        for got, want in zip(self.pll_b, asks):
            self.assertAlmostEqual(ppb(*got), ppb(*want), delta=1e-3)
        steering = [ppb(*a) for a in self.pll_b]
        self.assertGreater(max(steering), 35)
        self.assertLess(steering[-1], -100)

    @unittest.skipUnless(os.path.exists(REC45), "no recording in data/")
    def test_cold_search(self):
        """The search on the box (message 'C' 16 64) against the 45 s recording, the harness
        answering the app's snapshot requests from it: the satellites armed are real (the model's
        search at the arm's start finds chip 0 there, within a few samples, at the arm's Doppler
        within a bin), and include the strongest; the 16-block pass arms two, and the search goes
        on at 32 blocks until four channels track."""
        offset = 400000                                     # the engine's index 0, in the file
        lines = [f"e {REC45} {offset}", "m " + b"C\x10\x40".hex()]
        for _ in range(50):                             # the search's state as it goes
            lines += ["n 10000", "m " + b"K".hex()]
        writes, replies = self.run_app(lines)
        # A candidate's check is one Doppler, the cold pass's (the index 0 throughout).
        refining = [r[2] for r in replies[1:] if r[0] == 2]
        self.assertTrue(refining)
        self.assertEqual(set(refining), {0})
        k = replies[-1]
        self.assertEqual(k[0], 0)                                         # done
        armed = list(k[12:12 + k[11]])
        self.assertIn(15, armed)                                          # the strongest
        self.assertEqual(len(armed), 4)
        self.assertEqual([r[11] for r in replies[1:] if r[19] == 32][0], 2)
        self.assertEqual(k[19], 32)
        arms = [(c, b) for op, c, b in decode_commands(writes) if op == 0]
        self.assertEqual(len(arms), len(armed))
        for (c, bits), prn in zip(arms, armed):
            W, Wc = field(bits, 0, 32), field(bits, 32, 64)
            g1, g2, win = field(bits, 64, 74), field(bits, 74, 84), field(bits, 84, 87)
            start1 = field(bits, 87, 119) + 1
            self.assertEqual((g1, g2, win), tuple(int(v) for v in tm.arm_values(prn)))
            fd = (Wc - (1 << 32) if Wc >= 1 << 31 else Wc) / (1 << 32) * FS
            self.assertAlmostEqual(W / (1 << 32) * FS, 1.023e6 * (1 + fd / 1575.42e6), delta=1e-3)
            x = g_model.read_hackrf_iq(REC45, count=10 * 4000, offset=offset + start1)
            f = g_model.acquire(x, FS, prn, doppler_max=abs(fd) + 500, doppler_step=125, blocks=10)
            # the model's grid runs from -(fd + 500); the bin nearest fd
            self.assertLessEqual(abs(f["doppler"] - fd), 250, f"PRN {prn}")
            e = min(f["epoch"], 4000 - f["epoch"])
            self.assertLessEqual(e, 4, f"PRN {prn}: chip 0 {f['epoch']} samples from the start")

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_missing_record(self):
        """drop2.bin with one of channel 0's records cut out at about 50 s (as a snapshot pauses
        the stream): the app counts the missed period, so the fixes after the gap equal the uncut
        run's (one bit and its subframe lost, the anchor kept)."""
        words = np.fromfile(RUN, dtype="<u2")
        start = after_snapshots(words)
        n_first = int(words[start + 1]) | (int(words[start + 2]) << 16)
        cut = None
        i = start
        while i + 21 <= len(words):
            n0 = int(words[i + 1]) | (int(words[i + 2]) << 16)
            if (words[i] >> 8) == 0 and ((n0 - n_first) & M32) > 50 * FS:
                cut = i
                break
            i += 21
        self.assertIsNotNone(cut)
        path = os.path.join(self.tmp, "cut.bin")
        np.concatenate((words[:cut], words[cut + 21:])).tofile(path)
        base, cutrun = self.replay(RUN, RUN_ARMS), self.replay(path, RUN_ARMS)
        fix = lambda r: (int.from_bytes(r[0:4], "little"), np.frombuffer(r[20:44], "<f8"),
                         np.frombuffer(r[12:20], "<f8")[0])
        n0, x0, t0 = fix(base[-1])
        n1, x1, t1 = fix(cutrun[-1])
        self.assertGreater(n1, 40)
        self.assertLess(np.linalg.norm(x1 - x0), 1e-3)
        self.assertAlmostEqual(t1, t0, delta=1e-12)
        self.assertEqual(int.from_bytes(cutrun[-1][120:124], "little"), 0)  # none rejected

    # PRN 2's almanac entry and the time it serves (2026-10-07; test_gps_nav.py's pages), every
    # other PRN the same orbit called unhealthy: the search has PRN 2 alone to predict.
    ALM_KEYS = ("toa", "e", "i0", "omegadot", "sqrta", "omega0", "omega", "m0")

    def warm_setup(self, healthy, have, bias=0.0):
        """The held position under PRN 2, the almanac (message 'W', entries with have as given)
        and GPS time with the carrier offset bias; also that position and time."""
        alm = nav.decode_almanac([to_bits(w, 24) for w in PAGES[1]])
        t = nav.field([to_bits(w, 24) for w in PRN2[0]], (2, 1, 17)) * 6.0
        pos = nav.sat_position(alm, t)[0]
        rx = pos / np.linalg.norm(pos) * 6371e3
        lines = ["m " + (b"U" + struct.pack("<3d", *rx)).hex()]
        for prn in range(1, 33):
            e = struct.pack("<8fBBxx", *(alm[k] for k in self.ALM_KEYS),
                            0 if (prn == 2 and healthy) else 1, have)
            lines.append("m " + (b"W" + bytes([prn - 1]) + e).hex())
        lines.append("m " + (b"W" + struct.pack("<dd", t, bias)).hex())
        return lines, alm, t, rx

    def test_search_goes_cold_after_predicted_passes(self):
        """An almanac that calls every PRN unhealthy (a spoofed page passes parity like any other):
        the predicted passes search nothing, up to their deepest (128 blocks), and the search then
        goes on over every PRN at every Doppler (64 blocks) instead of stopping, so no almanac
        entry keeps a satellite from the search."""
        lines, _, _, _ = self.warm_setup(healthy=False, have=2)
        lines += ["m " + b"C\x08\x10".hex(), "n 2000", "m " + b"K".hex()]
        k = self.run_app(lines)[1][-1]
        # Cold at PRN 1's first Doppler (no snapshot comes), not predicted, at most 64 blocks.
        self.assertEqual((k[0], k[1], k[2], k[19], k[20]), (1, 0, 0, 64, 0))

    def test_search_marks_the_entry_that_found_its_satellite(self):
        """PRN 2 synthesized at the Doppler its almanac entry predicts (with the carrier offset the
        warm start gives): the predicted pass searches its three bins alone, finds and arms it, and
        its entry reads back marked as having led to the satellite (have 3), which hackrf_gps keeps
        it for, where one carried and never confirmed reads back as sent (2). One satellite is too
        few: the predicted passes go on to their deepest, then the search goes cold."""
        _, alm, t, rx = self.warm_setup(healthy=True, have=2)
        doppler = 1250.0                                 # a bin's center
        bias = doppler - freq.expected(alm, t, rx)
        lines, _, _, _ = self.warm_setup(healthy=True, have=2, bias=bias)
        x = g_model.synthesize(FS, 2_000_000, [dict(prn=2, code_phase=100.0, doppler=doppler,
                                                    cn0=48.0)], seed=3)
        iq = np.empty(2 * len(x), dtype=np.int8)
        iq[0::2], iq[1::2] = x.real, x.imag
        path = os.path.join(self.tmp, "prn2.iq")
        iq.tofile(path)
        lines = [f"e {path} 0"] + lines + ["m " + b"C\x08\x10".hex()]
        for _ in range(12):
            lines += ["n 500", "m " + b"K".hex()]
        lines += ["m " + (b"W" + bytes([1])).hex(), "m " + (b"W" + bytes([2])).hex()]
        replies = self.run_app(lines)[1]
        k = replies[-3]
        self.assertEqual(list(k[12:12 + k[11]]), [2])                     # PRN 2 armed
        self.assertEqual((k[0], k[19], k[20]), (1, 64, 0))                # then cold
        self.assertEqual(replies[-2][33], 3)                               # PRN 2's entry
        self.assertEqual(replies[-1][33], 2)                               # PRN 3's, unchanged

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_search_from_the_almanac(self):
        """drop2.bin with its position held and the search started at 70 s, when the app has
        almanac pages (subframes 4 and 5), estimates and a usable channel: it predicts first, then
        stops at the first PRN worth a search (the replay sends no snapshot), on the first of the
        three Doppler bins around its prediction; nav.py's almanac, freq.expected and the offset
        the tracked satellites show (their loops' Doppler less their expected carrier) give the
        same PRN and bin, within one. The almanac is not yet whole, so the pass is not deepened."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        at = 70
        # The almanac pages the app has read by then (a page is read with its tenth word).
        alm = {}
        for c, prn, _ in RUN_ARMS[:4]:
            p = z[f"{prn}_IP"] + 1j * z[f"{prn}_QP"]
            edge, _ = nav.bit_sync(p[1000:])
            edge = (edge + 1000) % nav.BIT_MS
            hard = [int(v < 0) for v in nav.bits(p.astype(complex), edge)]
            for sf in nav.subframes(hard):
                if sf["id"] in (4, 5) and (edge + 20 * (sf["bit"] + 300)) / 1000 < at - 2:
                    b = [v ^ int(sf["inverted"]) for v in hard]
                    ws = nav.words(b, sf["bit"])
                    a = ws and nav.decode_almanac(ws)
                    if a:
                        alm[a["prn"]] = a
        self.assertTrue(alm)
        # GPS time then, and the offset: PRN 27's loops (its carrier word) less its expected carrier.
        s27 = sats[0]
        n_r = int(z["27_n0"][0] + at * FS)
        t = fix.transmit_time(z, 27, s27["subframes"][0], n_r)
        j = int(np.searchsorted(z["27_n0"], n_r))
        wc = z["27_Wc"][j]
        doppler = (wc - 2 ** 32 if wc >= 2 ** 31 else wc) / 2 ** 32 * FS
        bias = doppler - freq.expected(s27["eph"], t, rx)
        tracked = {prn for _, prn, _ in RUN_ARMS}
        want = None
        for prn in range(1, 33):
            if prn in tracked:
                continue
            if prn not in alm:
                want = (prn, [0])
                break
            a = alm[prn]
            pos, _ = nav.sat_position(a, t)
            if a["health"] or np.linalg.norm(pos - rx) > 25.3e6:
                continue
            k = int(np.floor((freq.expected(a, t, rx) + bias + 4500) / 250 + 0.5))
            want = (prn, [k - 2, k - 1, k])
            break
        hold = b"U" + rx.astype("<f8").tobytes()
        k = self.replay(RUN, RUN_ARMS, b"K", [hold], at={at: b"C\x10\x20"})[-1]
        self.assertEqual(k[0], 1)                       # cold, waiting for a snapshot
        self.assertEqual(k[1] + 1, want[0])
        self.assertIn(k[2], want[1])
        # Pages for some PRNs only: the pass keeps the cold depth (not flagged predicted).
        self.assertLess(len(alm), 32)
        self.assertEqual(k[20], 0)

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_search_located(self):
        """No position, an almanac and GPS time: the search started at 70 s locates the receiver
        from the tracked satellites' Dopplers (their almanac entries here from their ephemerides,
        every PRN without one unhealthy) and searches only what the kept points see. PRN 1, below
        the horizon, is left out; PRN 2, about 3 degrees up, is the first searched, from a bin at
        most 6 below the one its almanac gives at the true position, and the pass is flagged
        predicted (the almanac is whole)."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = {prn: fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]}
        rx, _ = freq.survey(z, list(sats.values()), FS, int(FS))
        alm = {}
        for c, prn, _ in RUN_ARMS[:4]:
            p = z[f"{prn}_IP"] + 1j * z[f"{prn}_QP"]
            edge, _ = nav.bit_sync(p[1000:])
            edge = (edge + 1000) % nav.BIT_MS
            hard = [int(v < 0) for v in nav.bits(p.astype(complex), edge)]
            for sf in nav.subframes(hard):
                if sf["id"] in (4, 5):
                    ws = nav.words([v ^ int(sf["inverted"]) for v in hard], sf["bit"])
                    a = ws and nav.decode_almanac(ws)
                    if a:
                        alm[a["prn"]] = a
        alm.update({prn: dict(sat["eph"], toa=sat["eph"]["toe"], health=0)
                    for prn, sat in sats.items()})
        keys = ("toa", "e", "i0", "omegadot", "sqrta", "omega0", "omega", "m0")
        entry = {prn: struct.pack("<8fBBxx", *(alm[prn][k] for k in keys), alm[prn]["health"], 1)
                 if prn in alm else struct.pack("<8fBBxx", *[0.0] * 8, 1, 1)
                 for prn in range(1, 33)}
        t2 = fix.transmit_time(z, 27, sats[27]["subframes"][0], int(z["27_n0"][0] + 2 * FS))
        setup = [b"W" + bytes([prn - 1]) + e for prn, e in entry.items()]
        setup.append(b"W" + struct.pack("<dd", t2, 0.0))
        k = self.replay(RUN, RUN_ARMS, b"K", setup, at={2: b"C\x00\x00", 70: b"C\x10\x20"})[-1]
        # The truth: PRN 2's bin from its almanac at the true position, with the offset PRN 27's
        # loops show (its carrier word less its expected carrier).
        at = 70
        n_r = int(z["27_n0"][0] + at * FS)
        t = fix.transmit_time(z, 27, sats[27]["subframes"][0], n_r)
        wc = z["27_Wc"][int(np.searchsorted(z["27_n0"], n_r))]
        bias = (wc - 2 ** 32 if wc >= 2 ** 31 else wc) / 2 ** 32 * FS - \
            freq.expected(sats[27]["eph"], t, rx)
        true_bin = int(np.floor((freq.expected(alm[2], t, rx) + bias + 4500) / 250 + 0.5))
        self.assertEqual((k[0], k[1] + 1), (1, 2))     # cold, at PRN 2: PRN 1 left out
        self.assertLessEqual(k[2], true_bin)
        self.assertLessEqual(true_bin - k[2], 6)
        self.assertEqual(k[20], 1)

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_warm_start(self):
        """The warm start (message 'W'): an almanac entry written reads back unchanged, the status
        reply holds the position and the offset, and with GPS time given the first search predicts
        before any channel is usable (2 s in, no ephemeris yet): it skips PRN 1, below the horizon,
        and stops at the PRN and bins nav.py predicts; without the time it starts cold at PRN 1;
        with PLL B steered the predicted bins move with the steering's shift of the carriers."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        alm = {}
        for c, prn, _ in RUN_ARMS[:4]:
            p = z[f"{prn}_IP"] + 1j * z[f"{prn}_QP"]
            edge, _ = nav.bit_sync(p[1000:])
            edge = (edge + 1000) % nav.BIT_MS
            hard = [int(v < 0) for v in nav.bits(p.astype(complex), edge)]
            for sf in nav.subframes(hard):
                if sf["id"] in (4, 5):
                    ws = nav.words([v ^ int(sf["inverted"]) for v in hard], sf["bit"])
                    a = ws and nav.decode_almanac(ws)
                    if a:
                        alm[a["prn"]] = a
        self.assertIn(1, alm)
        at = 2
        s27 = sats[0]
        n_r = int(z["27_n0"][0] + at * FS)
        t = fix.transmit_time(z, 27, s27["subframes"][0], n_r)
        wc = z["27_Wc"][int(np.searchsorted(z["27_n0"], n_r))]
        bias = (wc - 2 ** 32 if wc >= 2 ** 31 else wc) / 2 ** 32 * FS - \
            freq.expected(s27["eph"], t, rx)
        entry = {prn: struct.pack("<8fBBxx", *(a[k] for k in ("toa", "e", "i0", "omegadot", "sqrta",
                                                             "omega0", "omega", "m0")), a["health"], 1)
                 for prn, a in alm.items()}
        setup = [b"U" + rx.astype("<f8").tobytes()]
        setup += [b"W" + bytes([prn - 1]) + e for prn, e in entry.items()]
        timed = setup + [b"W" + struct.pack("<dd", t, bias)]
        one = sorted(entry)[0]
        # Read back in the first second, before the app decodes that page itself.
        self.assertEqual(self.replay(RUN, RUN_ARMS, b"W" + bytes([one - 1]), timed)[0], entry[one])
        w = self.replay(RUN, RUN_ARMS, b"W", timed)[-1]
        self.assertEqual(w[:24], rx.astype("<f8").tobytes())
        self.assertEqual((w[24], w[25]), (3, 2))       # given; the offset now this run's own
        tracked = {prn for _, prn, _ in RUN_ARMS}
        want = None
        for prn in range(1, 33):
            if prn in tracked:
                continue
            if prn not in alm:
                want = (prn, [0])
                break
            a = alm[prn]
            pos, _ = nav.sat_position(a, t)
            if a["health"] or np.linalg.norm(pos - rx) > 25.3e6:
                continue
            k = int(np.floor((freq.expected(a, t, rx) + bias + 4500) / 250 + 0.5))
            want = (prn, [k - 2, k - 1, k])
            break
        self.assertGreater(want[0], 1)
        k = self.replay(RUN, RUN_ARMS, b"K", timed, at={at: b"C\x10\x20"})[-1]
        self.assertEqual((k[0], k[1] + 1), (1, want[0]))
        self.assertIn(k[2], want[1])
        k = self.replay(RUN, RUN_ARMS, b"K", setup, at={at: b"C\x10\x20"})[-1]
        self.assertEqual((k[0], k[1] + 1, k[2]), (1, 1, 0))         # cold: PRN 1, every Doppler
        # Steered (PLL B 125 ppb up a second before): the samples' carriers move by -L1 times
        # that, about -197 Hz, and the prediction with them, though the loops' offset does not.
        prn = want[0]
        a = alm[prn]
        shift = -1575.42e6 * 4 / (32 * 1e6)
        ks = int(np.floor((freq.expected(a, t, rx) + bias + shift + 4500) / 250 + 0.5))
        self.assertNotEqual(ks, want[1][2])
        b = b"B" + struct.pack("<III", 32, 4, 1000000)
        k = self.replay(RUN, RUN_ARMS, b"K", timed, at={at - 1: b, at: b"C\x10\x20"})[-1]
        self.assertEqual((k[0], k[1] + 1), (1, prn))
        self.assertEqual(k[2], ks - 1)                                # the first of its three

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_many_missing_records(self):
        """drop2.bin with 1% of channel 0's records cut out from 10 s, single ones at random (as
        snapshots cut them, about ten a second): no bit is dropped, so the channel reads every
        subframe the uncut run reads, or all but one."""
        words = np.fromfile(RUN, dtype="<u2")
        start = after_snapshots(words)
        n_first = int(words[start + 1]) | (int(words[start + 2]) << 16)
        rec = np.arange(start, len(words) - 20, 21)
        n0 = (words[rec + 1].astype(np.int64) | (words[rec + 2].astype(np.int64) << 16)) - n_first
        ours = rec[((words[rec] >> 8) == 0) & ((n0 & M32) > 10 * FS)]
        cut = np.random.default_rng(9).choice(ours, len(ours) // 100, replace=False)
        keep = np.ones(len(words), bool)
        for i in cut:
            keep[i:i + 21] = False
        path = os.path.join(self.tmp, "cuts.bin")
        words[keep].tofile(path)
        subframes = lambda r: int.from_bytes(r[64:66], "little")
        base, cutrun = self.replay(RUN, RUN_ARMS), self.replay(path, RUN_ARMS)
        self.assertGreater(subframes(base[-1]), 10)
        self.assertGreaterEqual(subframes(cutrun[-1]), subframes(base[-1]) - 1)

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_manager_drops_a_lost_channel(self):
        """drop2.bin with channel 0's sums replaced by noise from 40 s (PRN 27 blocked), the manager
        on without a search (message 'C' 0 0: a replay brings no snapshots): unlocked, its loops
        wander off on the noise, so it is stopped (not re-armed: its record count runs on to the
        stop) once 10 s unlocked and further than 20 Hz plus 1 Hz a second from its last locked
        Doppler; the other channels go on."""
        words = np.fromfile(RUN, dtype="<u2")
        start = after_snapshots(words)
        n_first = int(words[start + 1]) | (int(words[start + 2]) << 16)
        rng = np.random.default_rng(5)
        i, changed = start, 0
        while i + 21 <= len(words):
            n0 = int(words[i + 1]) | (int(words[i + 2]) << 16)
            if (words[i] >> 8) == 0 and ((n0 - n_first) & M32) > 40 * FS:
                for k in range(6):
                    v = int(rng.normal(0, 3e4)) & M32
                    words[i + 9 + 2 * k], words[i + 10 + 2 * k] = v & 0xFFFF, v >> 16
                changed += 1
            i += 21
        self.assertGreater(changed, 40000)
        path = os.path.join(self.tmp, "blocked.bin")
        words.tofile(path)
        replies = self.replay(path, RUN_ARMS, b"Q", [b"C\x00\x00"])   # the manager alone
        q = replies[-1]
        active = [q[16 * c] for c in range(4)]
        self.assertEqual(active, [0] + [prn for _, prn, _ in RUN_ARMS[1:4]])
        records = [int.from_bytes(r[4:8], "little") for r in replies if r[0]]
        self.assertTrue(all(b > a for a, b in zip(records, records[1:])))   # never restarted
        stop = next(i for i, r in enumerate(replies) if i > 10 and not r[0])
        self.assertTrue(48 <= stop <= 52)                             # 40 s, then 10 unlocked

    def replay(self, path, arms, query=b"P", setup=(), at=None):
        """The stream through the app, with its reply to query (navigation by default) after each
        second; setup, messages before the arms."""
        words = np.fromfile(path, dtype="<u2")
        start = after_snapshots(words)
        lines = ["m " + (b"F" + round(LO_OFFSET * 1e6).to_bytes(4, "little", signed=True)).hex()]
        lines += ["m " + m.hex() for m in setup]
        lines += ["m " + arm_message(c, prn, d * 1000, 0).hex() for c, prn, d in arms]
        step = len(arms) * 21 * 2 * 1000             # a second of records, in bytes
        remaining = 2 * (len(words) - start)
        lines.append(f"f {path} {2 * start} {step}")
        remaining -= step
        second = 1
        while remaining > 0:
            lines += ["m " + query.hex()]
            if at and second in at:
                lines += ["m " + at[second].hex()]
            lines += [f"f - 0 {step}"]
            remaining -= step
            second += 1
        lines.append("m " + query.hex())
        replies = self.run_app(lines)[1][1 + len(setup) + len(arms):]
        if at:  # the extra messages' replies out: each follows its second's query reply
            skip = {s + i for i, s in enumerate(sorted(at))}
            replies = [r for i, r in enumerate(replies) if i not in skip]
        return replies

    @unittest.skipUnless(os.path.exists(BAD_RUN), "no device run in data/")
    def test_replay_bad_geometry(self):
        last = self.replay(BAD_RUN, BAD_RUN_ARMS)[-1]
        self.assertEqual(int.from_bytes(last[0:4], "little"), 0)            # no fix kept
        self.assertGreater(int.from_bytes(last[120:124], "little"), 20)     # rejected
        self.assertIn(last[124], (2, 3, 4))

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_navigation(self):
        replies = self.replay(RUN, RUN_ARMS)
        fixes = {}
        for r in replies:
            count = int.from_bytes(r[0:4], "little")
            if count:
                fixes[count] = r
        self.assertGreater(len(fixes), 50)

        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        first = int(z["27_n0"][0])
        for r in fixes.values():
            self.assertEqual((r[8], r[10]), (4, 0x0F))
            n_r = first + ((int.from_bytes(r[4:8], "little") - first) & M32)
            t_app = np.frombuffer(r[12:20], "<f8")[0]
            x_app = np.frombuffer(r[20:44], "<f8")
            dop_app = np.frombuffer(r[44:64], "<f4")
            t_tx = [fix.transmit_time(z, s["prn"], s["subframes"][0], n_r) for s in sats]
            t_rx = max(t_tx) + 0.075
            pos, bias, _, _ = nav.solve(t_tx, [s["eph"] for s in sats], t_rx)
            # The same records, the same period for each transmit time, the same arithmetic.
            self.assertLess(np.linalg.norm(x_app - pos), 1e-3)
            self.assertAlmostEqual(t_app, t_rx - bias, delta=1e-12)
            at = [nav.sat_position(s["eph"], t)[0] for s, t in zip(sats, t_tx)]
            want = nav.dop(pos, at)                    # G, P and T; the app sends H and V as 0
            np.testing.assert_allclose(dop_app[[0, 1, 4]], [want[i] for i in (0, 1, 4)], rtol=1e-3)
            np.testing.assert_array_equal(dop_app[2:4], [0, 0])
        last = replies[-1]
        self.assertEqual(len(last), 128)
        self.assertEqual((int.from_bytes(last[120:124], "little"), last[124]), (0, 0))
        for c, prn, _ in RUN_ARMS:
            e = last[64 + 8 * c:72 + 8 * c]
            subframes, flags = int.from_bytes(e[0:2], "little"), e[2]
            if prn in (27, 16, 8, 4):
                self.assertGreaterEqual(subframes, 12)
                self.assertEqual(flags & 7, 7)
            else:
                self.assertEqual(subframes, 0)

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_frequency(self):
        """The app's frequency estimate each second against freq.py's on the same windows, with
        the same held position given to both."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        replies = self.replay(RUN, RUN_ARMS, b"R", [b"U" + rx.astype("<f8").tobytes()])
        estimates = {}
        for r in replies:
            self.assertEqual(len(r), 64)
            if r[18] == 0:
                estimates[int.from_bytes(r[0:4], "little")] = r
        # One reply a second of records, one estimate a second: each seen, so freq.Combiner can
        # follow the app's weights.
        counts = sorted(estimates)
        self.assertEqual(counts, list(range(counts[0], counts[0] + len(counts))))
        self.assertEqual(counts[0], 1)
        combiner = freq.Combiner()
        # Subframes 1-3 take up to 30 s; the run is 90 s.
        self.assertGreater(len(estimates), 40)
        first = int(z["27_n0"][0])
        for r in estimates.values():
            self.assertEqual(r[19], 3)                                   # position given
            end = first + ((int.from_bytes(r[4:8], "little") - first) & M32)
            bounds = np.array([end - int(FS), end])
            ppb_app = np.frombuffer(r[24:52], "<f4")
            lock_app = np.frombuffer(r[52:59], "i1") / 100.0
            valid, used = r[16], r[17]
            self.assertEqual(valid & ~0x0F, 0)                           # GPS with ephemeris only
            for k, (c, prn, _) in enumerate(RUN_ARMS[:4]):
                ppb, _, lock = freq.satellite_offsets(z, sats[k], bounds, rx, LO_OFFSET)[0]
                self.assertAlmostEqual(lock_app[c], lock, delta=0.05)
                if valid & (1 << c):
                    # The same cycles, the expected frequency within microseconds of the same
                    # transmit time.
                    self.assertAlmostEqual(ppb_app[c], ppb, delta=1e-4)
        for count in counts:
            r = estimates[count]
            valid, used = r[16], r[17]
            ppb_app = np.frombuffer(r[24:52], "<f4").astype(np.float64)
            est, keys = combiner({c: ppb_app[c] for c in range(7) if valid & (1 << c)})
            self.assertEqual(sum(1 << c for c in keys), used)
            self.assertAlmostEqual(np.frombuffer(r[8:12], "<f4")[0], est, delta=1e-4)

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_frequency_survey(self):
        """With no position given, the estimate waits for the first kept fix, then uses the mean of
        the fixes kept so far (60 hold it; the replay keeps 58, every one counted)."""
        replies = self.replay(RUN, RUN_ARMS, b"R")
        states = [(r[18], r[19], int.from_bytes(r[20:24], "little")) for r in replies]
        self.assertIn((1, 0, 0), states)                                 # no position yet
        last = replies[-1]
        self.assertEqual((last[18], last[19]), (0, 1))                   # estimated, surveying
        self.assertEqual(int.from_bytes(last[20:24], "little"), 58)
        self.assertEqual(last[17], 0x0F)
        self.assertLess(abs(np.frombuffer(last[8:12], "<f4")[0] + 126.4), 1.0)

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_frequency_seeded_survey(self):
        """A rough position (a grid square's center, here about 3 km off) with the survey: with two
        satellites, which never give a fix, the estimate runs from it (the survey alone gives none);
        with four, the first kept fix replaces it and the survey ends where it would without it."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        seed = rx + np.array([2000.0, -2000.0, 1000.0])
        msg = b"U" + seed.astype("<f8").tobytes() + (60).to_bytes(4, "little")
        two = self.replay(RUN, RUN_ARMS[:2], b"R", [msg])
        self.assertFalse([r for r in self.replay(RUN, RUN_ARMS[:2], b"R") if r[18] == 0])
        est = [r for r in two if r[18] == 0]
        self.assertGreater(len(est), 30)
        for r in est:
            self.assertEqual((r[19], int.from_bytes(r[20:24], "little")), (1, 0))
        self.assertLess(abs(np.frombuffer(est[-1][8:12], "<f4")[0] + 126.4), 3.0)
        seeded = self.replay(RUN, RUN_ARMS, b"R", [msg])
        plain = self.replay(RUN, RUN_ARMS, b"R")
        self.assertEqual(int.from_bytes(seeded[-1][20:24], "little"), 58)
        # The same but for the weights' memory of the first windows (their variances).
        ppb = lambda r: np.frombuffer(r[8:12], "<f4")[0]
        self.assertLess(abs(ppb(seeded[-1]) - ppb(plain[-1])), 0.01)

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_frequency_disagreeing(self):
        """Two satellites from a seed 300 km off disagree by far more than twice REJECT_PPB: nothing
        says which is wrong, so no estimate (result 4, none counted) and no nan, where the seed 3 km
        off gives one every second."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        far = rx + np.array([200e3, -200e3, 100e3])
        replies = self.replay(RUN, RUN_ARMS[:2], b"R", [b"U" + far.astype("<f8").tobytes()])
        both = [r for r in replies if r[16] == 0x03]
        self.assertGreater(len(both), 30)
        for r in both:
            self.assertEqual(r[18], 4)
            self.assertEqual(int.from_bytes(r[0:4], "little"), 0)
            self.assertFalse(np.isnan(np.frombuffer(r[8:12], "<f4")[0]))
        e, used = freq.Combiner()({0: 0.0, 1: 2.5})
        self.assertTrue(np.isnan(e))
        self.assertEqual(used, [])

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_three_satellites_replace_a_grid_square(self):
        """A grid square's center (here 3 km off) held, three satellites: marked rough ('H' with 1),
        their fix with the height held replaces it and the estimate meets the surveyed position's;
        not marked (a warm start's surveyed position), the seed stays and so does its error."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        h = nav.ecef_to_geodetic(rx)[2]
        seed = b"U" + (rx + np.array([2000.0, -2000.0, 1000.0])).astype("<f8").tobytes() + (60).to_bytes(4, "little")
        given = {int.from_bytes(r[4:8], "little"): np.frombuffer(r[8:12], "<f4")[0]
                 for r in self.replay(RUN, RUN_ARMS[:3], b"R", [b"U" + rx.astype("<f8").tobytes()])
                 if r[18] == 0}
        off = {}
        for rough in (1, 0):
            reps = self.replay(RUN, RUN_ARMS[:3], b"R", [seed, b"H" + struct.pack("<f", h) + bytes([rough])])
            d = [abs(np.frombuffer(r[8:12], "<f4")[0] - given[int.from_bytes(r[4:8], "little")])
                 for r in reps[-20:] if r[18] == 0 and int.from_bytes(r[4:8], "little") in given]
            self.assertGreater(len(d), 10)
            off[rough] = np.median(d)
        self.assertLess(off[1], 0.05)
        self.assertGreater(off[0], 3 * off[1])

    @unittest.skipUnless(os.path.exists(RUN), "no device run in data/")
    def test_replay_three_satellites_height_held(self):
        """Three satellites and no position: with the height held ('H'), a fix from the three is
        the held position (not kept: no fix counted, none surveyed), and the estimate runs from it,
        against the same three satellites with the surveyed position given, window by window. Here
        the first such fixes find the mirror solution (3,600 km off), the satellites disagree, and
        the next fix, from the mirror's reflection, finds the right one. The right height gives the
        given position's estimate within 0.05 ppb; 1 km off moves it, by under 1.5 ppb here."""
        z = fix.load_stream(RUN, ",".join(f"{c}:{prn}" for c, prn, _ in RUN_ARMS[:4]))
        sats = [fix.satellite(z, prn) for _, prn, _ in RUN_ARMS[:4]]
        rx, _ = freq.survey(z, sats, FS, int(FS))
        h = nav.ecef_to_geodetic(rx)[2]
        held = lambda height: b"H" + struct.pack("<f", height)
        given = {int.from_bytes(r[4:8], "little"): np.frombuffer(r[8:12], "<f4")[0]
                 for r in self.replay(RUN, RUN_ARMS[:3], b"R", [b"U" + rx.astype("<f8").tobytes()])
                 if r[18] == 0}
        fixes = self.replay(RUN, RUN_ARMS[:3], b"P", [held(h)])
        self.assertIn(5, [r[124] for r in fixes])                        # three, the height held
        self.assertEqual(int.from_bytes(fixes[-1][0:4], "little"), 0)    # none kept
        diffs = {}
        for name, height in (("right", h), ("1 km off", h + 1000.0)):
            reps = self.replay(RUN, RUN_ARMS[:3], b"R", [held(height)])
            est = [r for r in reps if r[18] == 0]
            self.assertGreater(len(est), 30, name)
            first = next(i for i, r in enumerate(reps) if r[18] == 0)
            self.assertIn(4, [r[18] for r in reps[:first]], name)       # the mirror, disagreeing
            for r in est:
                self.assertEqual((r[19], int.from_bytes(r[20:24], "little")), (1, 0))
            d = [abs(np.frombuffer(r[8:12], "<f4")[0] - given[int.from_bytes(r[4:8], "little")])
                 for r in est if int.from_bytes(r[4:8], "little") in given]
            self.assertGreater(len(d), 30, name)
            diffs[name] = np.median(d)
        self.assertLess(diffs["right"], 0.05)
        self.assertGreater(diffs["1 km off"], diffs["right"])
        self.assertLess(diffs["1 km off"], 1.5)

if __name__ == "__main__":
    unittest.main()
