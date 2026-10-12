#!/usr/bin/env python3
#
# SPDX-License-Identifier: BSD-3-Clause

"""The C loops (app/gps_loops.c) against the bit-exact Python model (gps_loop_int.IntLoops): the same
dumps in, the same code and carrier words out, every period.

The dumps are replayed open loop, so the inputs are identical by construction: synthetic phasors with
data bits, a residual frequency and noise, for GPS and SBAS symbol lengths; sums at the correlator's
full 28-bit width; all-zero sums; and the recorded dumps in data/ when they are present.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fpga-build"))

from dsp import gps_loop_int as li  # noqa: E402

KEYS = ("IE", "QE", "IP", "QP", "IL", "QL")
FS = 4e6
RECORDED = os.path.join(HERE, "data", "l1-45s.int.npz")


def synthetic(n, symbol_ms, amp, freq_hz, seed):
    """n dumps: prompt a rotating phasor with data bits, early and late 0.6 of it, plus noise."""
    rng = np.random.default_rng(seed)
    t = np.arange(n) * 1e-3
    bits = rng.choice([-1, 1], size=n // symbol_ms + 1).repeat(symbol_ms)[:n]
    p = amp * bits * np.exp(1j * (2 * np.pi * freq_hz * t + 0.4))
    out = []
    for scale in (0.6, 1.0, 0.55):
        z = scale * p + amp * 0.3 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
        out += [np.round(z.real).astype(np.int64), np.round(z.imag).astype(np.int64)]
    return np.stack(out, axis=1)


class TestLoopsC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.exe = os.path.join(cls.tmp, "loops_replay")
        src = os.path.join(HERE, "app")
        subprocess.run(["cc", "-std=c11", "-O2", "-Wall", "-Wextra", "-Werror", "-o", cls.exe,
                        os.path.join(src, "gps_loops.c"), os.path.join(src, "loops_replay.c"), "-lm"],
                       check=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def compare(self, sums, doppler, symbol_ms, fll_ms=300, sync_ms=1000):
        text = f"{FS} {doppler} {symbol_ms} {fll_ms} {sync_ms}\n"
        text += "\n".join(" ".join(str(int(v)) for v in row) for row in sums) + "\n"
        out = subprocess.run([self.exe], input=text, capture_output=True, text=True, check=True)
        lines = out.stdout.splitlines()
        k = li.gains(FS, symbol_ms)
        self.assertEqual([int(v) for v in lines[0].split()],
                         [k[n] for n in ("fll", "sync", "pll_i", "pll_p", "assist", "dll")])
        loops = li.IntLoops(FS, doppler, symbol_ms=symbol_ms, fll_ms=fll_ms, sync_ms=sync_ms)
        for j, (row, line) in enumerate(zip(sums, lines[1:])):
            code, carrier = loops.update({name: int(v) for name, v in zip(KEYS, row)})
            edge = -1 if loops.edge is None else loops.edge
            self.assertEqual([int(v) for v in line.split()], [code, carrier, edge], f"dump {j}")
        self.assertEqual(len(lines) - 1, len(sums))

    def test_synthetic_gps(self):
        self.compare(synthetic(4000, 20, 2e5, 3.7, 1), 1234.5, 20)

    def test_synthetic_sbas(self):
        self.compare(synthetic(3000, 2, 1e5, -2.2, 2), -2500.0, 2)

    def test_full_width_sums(self):
        # A period's sums fill 28 bits; the symbol sums then pass 2^31.
        self.compare(synthetic(2000, 20, 2 ** 26, 40.0, 3), 4000.0, 20, fll_ms=100, sync_ms=400)

    def test_zero_sums(self):
        self.compare(np.zeros((1500, 6), dtype=np.int64), -900.0, 20, fll_ms=100, sync_ms=400)

    @unittest.skipUnless(os.path.exists(RECORDED), "recorded dumps not present")
    def test_recorded(self):
        data = np.load(RECORDED)
        for prn in (15, 131):
            if f"{prn}_IP" not in data:
                continue
            sums = np.stack([data[f"{prn}_{k}"] for k in KEYS], axis=1)[:6000]
            doppler = float(data[f"{prn}_doppler"][0])
            with self.subTest(prn=prn):
                self.compare(sums, doppler, 2 if prn >= 120 else 20)


if __name__ == "__main__":
    unittest.main()
