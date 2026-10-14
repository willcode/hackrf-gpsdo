#!/usr/bin/env python3
#
# SPDX-License-Identifier: BSD-3-Clause

"""The C navigation (app/gps_nav.c) against nav.py: its own square root, sine, cosine and arctangent
against numpy; word parity; the subframe fields; the almanac pages; and the fix and DOP from
transmit times.

The subframes are four satellites' broadcast words of 2026-10-07 (GPS week 2439), as the device
tracked them. The receiver is placed where all four stand high in the sky, with a clock bias, and its
transmit times are made from nav.py's own measurement model, so both solvers have an exact answer.
gps_nav.c calls no math library function (it has its own square root, sine and cosine): the M4 app
links none. The host build links -lm as every test build here does.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import nav  # noqa: E402

# (PRN, subframe ID, the ten 24-bit data words)
WORDS = [
    (27, 1, [0x8B061C, 0x6E9BA5, 0x61D000, 0x1D75FE, 0x8F6582, 0x415090, 0xC82105, 0x10545F, 0x000000, 0x000BC6]),
    (27, 2, [0x8B061C, 0x6E99AA, 0x10005B, 0x3B9922, 0xC1F220, 0x00B907, 0x6DCBB3, 0x0738A1, 0x0D5D2C, 0x545F32]),
    (27, 3, [0x8B061C, 0x6E9A2E, 0xFFD6B3, 0x03DA53, 0x004A26, 0xB038A7, 0x268F24, 0xDBDC2F, 0xFF9E40, 0x10F683]),
    (16, 1, [0x8B061C, 0x6E9BA5, 0x61D000, 0x1D75FE, 0x8F6582, 0x415090, 0xC821E9, 0x465460, 0x000039, 0x2D745C]),
    (16, 2, [0x8B061C, 0x6E99AA, 0x460858, 0x33A746, 0xE2FC14, 0x072607, 0xD99214, 0x147AA1, 0x0E0B6D, 0x54607E]),
    (16, 3, [0x8B061C, 0x6E9A2E, 0xFF198C, 0xB26246, 0xFFA227, 0x07A869, 0x183026, 0xB062E3, 0xFFA679, 0x4615CD]),
    (8, 1, [0x8B061C, 0x6E9BA5, 0x61D000, 0x1D75FE, 0x8F6582, 0x415090, 0xC8210A, 0x195460, 0x00FF9E, 0x2A9186]),
    (8, 2, [0x8B061C, 0x6E99AA, 0x19011D, 0x3D611A, 0x7F5EC6, 0x017706, 0x1E4B6B, 0x0642A1, 0x0E3DBB, 0x54607E]),
    (8, 3, [0x8B061C, 0x6E9A2E, 0x0006B1, 0xA20996, 0x003B26, 0x5186AA, 0x267816, 0x978819, 0xFF9DF5, 0x19F71E]),
    (4, 1, [0x8B061C, 0x6E9BA5, 0x61D002, 0x1D75FE, 0x8F6582, 0x415090, 0xC821F6, 0x435460, 0x000000, 0x088AD5]),
    (4, 2, [0x8B061C, 0x6E99AA, 0x43FD17, 0x322C13, 0x459E99, 0xFD8902, 0x0BBC9A, 0x064CA1, 0x0D364C, 0x54607D]),
    (4, 3, [0x8B061C, 0x6E9A2E, 0x001735, 0x67703E, 0x001D27, 0xA1EB05, 0x28DD89, 0x582748, 0xFFA441, 0x43FED4]),
]

# Subframe 4 and 5 pages of 2026-10-07: the almanac of SV 1, 2, 3, 25,
# 26 and 21-23, and pages that hold none (SV IDs 57, 60, 61: other data).
PAGES = [
    [0x8B061C, 0x6E9B36, 0x4110CC, 0x7B08F3, 0xFD4100, 0xA10D89, 0xE15DE8, 0x0688DA, 0x9CE8FF, 0x14FFC4],
    [0x8B061C, 0x6E9DB7, 0x428C5C, 0x7B0BF3, 0xFD3D00, 0xA10D2C, 0xDAC87B, 0xE38F22, 0xD12747, 0x190041],
    [0x8B061C, 0x6EA034, 0x43383A, 0x7B22FC, 0xFD7900, 0xA10D7F, 0x09E519, 0x34F884, 0x422BF1, 0x1DFF46],
    [0x8B061C, 0x6E9D31, 0x595FED, 0x7B0337, 0xFD4F00, 0xA10DDF, 0x880BA2, 0x2813F1, 0x41D862, 0x31FFFF],
    [0x8B061C, 0x6E9FB3, 0x5A5EF3, 0x7BF6B2, 0xFD4300, 0xA10C93, 0x84D784, 0x1E3E59, 0x1490C7, 0xCC0008],
    [0x8B061C, 0x6F88B5, 0x550755, 0x7B0F06, 0xFD6000, 0xA10CC9, 0x09E11B, 0xF5D09F, 0x3049F4, 0x2CFFBD],
    [0x8B061C, 0x6F8B34, 0x565E4A, 0x7B09C5, 0xFD5700, 0xA10D61, 0x8CC8A2, 0xD8F1AD, 0xE04291, 0x0D0045],
    [0x8B061C, 0x6F8DB5, 0x5737CB, 0x7B1F52, 0xFD7600, 0xA10C7E, 0x08714D, 0x91B4BF, 0x45AC09, 0x5B002B],
    [0x8B061C, 0x6E9AB3, 0x79E00F, 0x58C42B, 0x0255D5, 0xE6FBBD, 0x80C004, 0x353004, 0x040868, 0xDEE4D0],
    [0x8B061C, 0x6F8AB1, 0x7CC032, 0x49B1CA, 0xD59B18, 0x6A1F2B, 0x71B88B, 0x4B7B05, 0xB36D47, 0x6A8E61],
    [0x8B061C, 0x6F8D33, 0x7D3B7C, 0x87A338, 0x7A90CB, 0xA3EA85, 0xF08B64, 0x0ABA05, 0xD254E8, 0xEA1282],
]
# PRN 2's subframes 1-3 (nav4), for its almanac's position against its ephemeris.
PRN2 = [
    [0x8B061C, 0x6F8925, 0x61D000, 0x1D75FE, 0x8F6582, 0x415090, 0xC821DA, 0x495460, 0x00003C, 0x18E154],
    [0x8B061C, 0x6F89A8, 0x49F85F, 0x311A24, 0xFDE54E, 0xF88608, 0xC4C562, 0x12F7A1, 0x0D9F3F, 0x54607E],
    [0x8B061C, 0x6F8A2C, 0x00A7DA, 0xD617E5, 0xFF9E27, 0x265CAC, 0x1959E3, 0x88A12F, 0xFFA695, 0x49E6A9],
]
ALM_FIELDS = ("toa", "e", "i0", "omegadot", "sqrta", "omega0", "omega", "m0")   # the C keeps no clock

FIELDS = ("tgd", "toc", "af0", "af1", "af2", "crs", "dn", "m0", "cuc", "e", "cus", "sqrta", "toe",
          "cic", "omega0", "cis", "i0", "crc", "omega", "omegadot", "idot")


def parity_bits(src, d29, d30):
    """The six parity bits of 24 source bits (D1 first), IS-GPS-200 20.3.5."""
    def x(*idx):
        v = 0
        for i in idx:
            v ^= src[i - 1]
        return v
    return [
        d29 ^ x(1, 2, 3, 5, 6, 10, 11, 12, 13, 14, 17, 18, 20, 23),
        d30 ^ x(2, 3, 4, 6, 7, 11, 12, 13, 14, 15, 18, 19, 21, 24),
        d29 ^ x(1, 3, 4, 5, 7, 8, 12, 13, 14, 15, 16, 19, 20, 22),
        d30 ^ x(2, 4, 5, 6, 8, 9, 13, 14, 15, 16, 17, 20, 21, 23),
        d30 ^ x(1, 3, 5, 6, 7, 9, 10, 14, 15, 16, 17, 18, 21, 22, 24),
        d29 ^ x(3, 5, 6, 8, 9, 10, 11, 13, 15, 19, 22, 23, 24),
    ]


def to_bits(v, n):
    return [(v >> (n - 1 - i)) & 1 for i in range(n)]


class TestNavC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.exe = os.path.join(cls.tmp, "nav_check")
        src = os.path.join(HERE, "app")
        subprocess.run(["cc", "-std=gnu11", "-O2", "-Wall", "-Wextra", "-Werror", "-o", cls.exe,
                        os.path.join(src, "gps_nav.c"), os.path.join(src, "gps_freq.c"),
                        os.path.join(src, "nav_check.c"), "-lm"],
                       check=True)
        cls.ephs = {}
        for prn, _, ws in WORDS:
            sf = nav.decode_subframe([to_bits(w, 24) for w in ws])
            cls.ephs.setdefault(prn, {}).update(sf)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def run_c(self, lines):
        out = subprocess.run([self.exe], input="\n".join(lines) + "\n", capture_output=True,
                             text=True, check=True)
        return [line.split() for line in out.stdout.splitlines()]

    def subframe_lines(self):
        return [f"sub {prn} " + " ".join(f"{w:06x}" for w in ws) for prn, _, ws in WORDS]

    def test_math(self):
        rng = np.random.default_rng(1)
        x = np.r_[rng.uniform(-100, 100, 300), rng.uniform(-1, 1, 100), 0.0, np.pi / 2, -np.pi]
        y = np.r_[rng.uniform(-100, 100, 300), rng.uniform(-1e-3, 1e-3, 100), 1.0, 0.0, 0.0]
        out = np.array(self.run_c([f"math {float(a)!r} {float(b)!r}" for a, b in zip(x, y)]), dtype=float)
        np.testing.assert_allclose(out[:, 0], np.sin(x), rtol=0, atol=4e-16)
        np.testing.assert_allclose(out[:, 1], np.cos(x), rtol=0, atol=4e-16)
        np.testing.assert_allclose(out[:, 2], np.arctan2(y, x), rtol=0, atol=1e-15)
        np.testing.assert_allclose(out[:, 3], np.sqrt(np.abs(x)), rtol=2.3e-16, atol=0)  # one ulp

    def test_parity(self):
        rng = np.random.default_rng(2)
        lines, want = [], []
        for _ in range(400):
            src = to_bits(int(rng.integers(1 << 24)), 24)
            d29, d30 = (int(v) for v in rng.integers(2, size=2))
            word = [b ^ d30 for b in src] + parity_bits(src, d29, d30)
            if rng.random() < 0.5:
                word[int(rng.integers(30))] ^= 1
            got = nav.parity_ok(word, d29, d30)
            want.append(["0", "000000"] if got is None
                        else ["1", f"{int(''.join(map(str, got)), 2):06x}"])
            lines.append(f"parity {int(''.join(map(str, word)), 2):08x} {d29} {d30}")
        self.assertEqual(self.run_c(lines), want)
        self.assertIn(["0", "000000"], want)

    def test_subframes(self):
        out = self.run_c(self.subframe_lines() + [f"eph {prn}" for prn in sorted(self.ephs)])
        self.assertEqual(len(out), len(WORDS) + len(self.ephs))
        for (prn, sid, ws), line in zip(WORDS, out):
            self.assertEqual(int(line[0]), sid)
            self.assertEqual(int(line[1]), nav.field([to_bits(w, 24) for w in ws], (2, 1, 17)))
        for prn, line in zip(sorted(self.ephs), out[len(WORDS):]):
            e = self.ephs[prn]
            self.assertEqual([int(v) for v in line[:6]],
                             [1, e["week"], e["health"], e["iodc"], e["iode"], e["iode"]])
            self.assertEqual(len(line), 6 + len(FIELDS))
            for name, v in zip(FIELDS, line[6:]):
                self.assertEqual(float(v), e[name], f"PRN {prn} {name}")

    def test_almanac(self):
        """Each page's entry as nav.py decodes it, in single precision; the pages of other data
        refused; and PRN 2's almanac within 3 km of its ephemeris over six hours."""
        lines = ["alm " + " ".join(f"{w:06x}" for w in ws) for ws in PAGES]
        out = self.run_c(lines)
        prns = []
        for ws, line in zip(PAGES, out):
            want = nav.decode_almanac([to_bits(w, 24) for w in ws])
            if want is None:
                self.assertEqual(line, ["0"])
                continue
            prns.append(want["prn"])
            self.assertEqual(int(line[0]), want["prn"])
            self.assertEqual(int(line[1]), want["health"])
            # Scaled in single precision on the M4: within a float's rounding of nav.py's doubles.
            for name, v in zip(ALM_FIELDS, line[2:]):
                self.assertLessEqual(abs(float(v) - want[name]), 2.4e-7 * abs(want[name]) + 1e-30,
                                     f"SV {want['prn']} {name}")
        self.assertEqual(prns, [1, 2, 3, 25, 26, 21, 22, 23])
        eph = {}
        for ws in PRN2:
            eph.update(nav.decode_subframe([to_bits(w, 24) for w in ws]))
        alm = nav.decode_almanac([to_bits(w, 24) for w in PAGES[1]])
        t0 = eph["tow"] * 6.0
        times = [t0 + dt for dt in (0, 3600, 6 * 3600)]
        got = self.run_c(lines + [f"almpos 2 {t!r}" for t in times])[len(lines):]
        for t, line in zip(times, got):
            x = np.array([float(v) for v in line[:3]])
            self.assertLess(np.linalg.norm(x - nav.sat_position(alm, t)[0]), 10.0)   # the floats
            self.assertLess(np.linalg.norm(x - nav.sat_position(eph, t)[0]), 3000.0)

    def synthetic(self, t_true, bias):
        """A receiver where the satellites stand high, and its transmit times at t_true."""
        prns = sorted(self.ephs)
        ephs = [self.ephs[p] for p in prns]
        mean = np.mean([nav.sat_position(e, t_true)[0] for e in ephs], axis=0)
        rx = mean / np.linalg.norm(mean) * 6371e3
        t_tx = []
        for e in ephs:
            t = t_true - 0.075
            for _ in range(20):
                _, clk = nav.sat_position(e, t)
                pos, clk = nav.sat_position(e, t - clk)
                rot = nav.OMEGA_E * np.linalg.norm(pos - rx) / nav.C
                sat = np.array([pos[0] * np.cos(rot) + pos[1] * np.sin(rot),
                                -pos[0] * np.sin(rot) + pos[1] * np.cos(rot), pos[2]])
                rng = np.linalg.norm(sat - rx)
                tropo = 2.47 / (np.sin(nav.elevation(rx, sat)) + 0.0121)
                t = t_true - (rng + tropo) / nav.C + clk
            t_tx.append(t)
        return prns, ephs, rx, t_tx, t_true + bias

    def test_solve_and_dop(self):
        tow = min(nav.field([to_bits(w, 24) for w in ws], (2, 1, 17)) for _, _, ws in WORDS)
        prns, ephs, rx, t_tx, t_rx = self.synthetic(tow * 6.0 + 10.0, 1.234e-3)
        pos, bias, res, els = nav.solve(t_tx, ephs, t_rx)
        self.assertTrue(np.all(els > 30), els)
        line = " ".join(f"{p} {float(t)!r}" for p, t in zip(prns, t_tx))
        out = self.run_c(self.subframe_lines() + [f"solve {float(t_rx)!r} {len(prns)} {line}"])[-1]
        x = np.array([float(v) for v in out[1:5]])
        self.assertLess(np.linalg.norm(x[:3] - pos), 1e-3)
        self.assertAlmostEqual(x[3] / nav.C, bias, delta=1e-12)
        # Against the truth, the times themselves limit: a double near 3.4e5 s resolves 6e-11 s, 1.7 cm
        # of range, and four high satellites (GDOP near 37) make that a decimeter.
        self.assertLess(np.linalg.norm(x[:3] - rx), 0.5)
        self.assertAlmostEqual(bias, 1.234e-3, delta=1e-9)
        np.testing.assert_allclose([float(v) for v in out[5:5 + len(prns)]], res, atol=1e-3)
        sats = [nav.sat_position(e, t)[0] for e, t in zip(ephs, t_tx)]
        # G, P and T; H and V are left out of the C (0) for room in the app. The C takes them from
        # the solve's last normal matrix: its satellites at the transmit time less their clocks and
        # turned with the Earth, against these at the transmit time (a few parts in 1e6 apart).
        got, want = [float(v) for v in out[5 + len(prns):]], nav.dop(pos, sats)
        np.testing.assert_allclose([got[i] for i in (0, 1, 4)], [want[i] for i in (0, 1, 4)], rtol=1e-5)
        self.assertEqual(got[2:4], [0.0, 0.0])

    def test_frequency_survives_a_nan_ephemeris(self):
        """A fifth satellite whose subframe 2 says sqrt(A) = 0 (any parity-valid subframe can be sent
        on the air) gives a NaN expected frequency. Its NaN offset once made its variance NaN; once
        its offset was finite again, its weight made the estimate NaN, every satellite's variance
        followed, and no estimate reached the steering until the channels were armed again. Such an
        ephemeris is now refused, and a non-finite offset is left out of the median and the
        variances whatever its cause."""
        spoof = [list(ws) for prn, _, ws in WORDS if prn == 27]
        spoof[1][7] &= 0xFFFF00   # sqrt(A): word 8's last 8 bits, then word 9
        spoof[1][8] = 0
        lines = self.subframe_lines() + ["sub 40 " + " ".join(f"{w:06x}" for w in ws) for ws in spoof]
        out = self.run_c(lines + ["eph 40", "eph 27"])
        self.assertEqual(float(out[-2][6 + FIELDS.index("sqrta")]), 0.0)
        self.assertEqual([out[-2][0], out[-1][0]], ["0", "1"])   # refused; the real one complete
        tow = min(nav.field([to_bits(w, 24) for w in ws], (2, 1, 17)) for _, _, ws in WORDS)
        _, _, rx, _, _ = self.synthetic(tow * 6.0 + 10.0, 0.0)
        windows = []
        for k in range(12):
            # Channel 0 is PRN 27's carrier against the zero sqrt(A) for three windows, then
            # against PRN 27's own ephemeris (the spoofer stopped, the next issue read).
            sats = [(27, 40 if k < 3 else 27), (16, 16), (8, 8), (4, 4)]
            windows.append(f"freq 0.15 {rx[0]!r} {rx[1]!r} {rx[2]!r} {tow * 6.0 + 10.0 + k!r} 4 "
                           + " ".join(f"{m} {e}" for m, e in sats))
        got = np.array(self.run_c(lines + windows)[len(lines):], dtype=float)
        self.assertTrue(np.all(np.isfinite(got[:, 0])), got[:, 0])          # an estimate each
        np.testing.assert_allclose(got[:, 0], 0.15, atol=1e-6)
        self.assertEqual([int(u) for u in got[:, 1]], [14] * 3 + [15] * 9)  # channel 0 back
        self.assertTrue(np.all(np.isfinite(got[:, 2:])), got[:, 2:])        # every variance

if __name__ == "__main__":
    unittest.main()
