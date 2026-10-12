#
# SPDX-License-Identifier: BSD-3-Clause

"""The integer loops' model: CORDIC accuracy, gains that fit the unit, and tracking that matches the
floating-point loops it stands for."""

import math
import unittest

import numpy as np

from dsp import gps_model as g, gps_track_model as tm, gps_loop_int as li


class TestPieces(unittest.TestCase):
    def test_cordic_angle_and_ratio(self):
        # Inputs lose PRE_SHIFT bits, so small vectors lose angle precision. The recordings' 1 ms prompt
        # sums have medians of 2e5 to 4e5 (signal and noise), where the error is 3 units or less.
        for r, bound in ((3e3, 64), (3e4, 8), (1e6, 4), (3e8, 4)):
            for a in np.linspace(-math.pi, math.pi, 499, endpoint=False):
                x, y = int(r * math.cos(a)), int(r * math.sin(a))
                ang, mag = li.cordic(x, y)
                self.assertLessEqual(abs(li.wrap(ang - round(a / (2 * math.pi) * li.FULL))), bound)
        _, m1 = li.cordic(400000, 300000)
        _, m2 = li.cordic(0, -250000)
        self.assertAlmostEqual(m1 / m2, 2.0, places=3)

    def test_fold_and_division(self):
        self.assertEqual(li.fold(li.QUARTER + 5), 5 - li.QUARTER)
        self.assertEqual(li.fold(-li.QUARTER), li.QUARTER)
        self.assertEqual(li.div_q14(-3, 4), -12288)
        self.assertEqual(li.div_q14(1, 3), 5461)
        self.assertEqual(li.div_q14(5, 0), 0)

    def test_gains_fit_16_bits(self):
        for sym in (2, 20):
            for v in li.gains(4e6, sym).values():
                self.assertTrue(0 <= v < 1 << 16)


class TestTracking(unittest.TestCase):
    def test_matches_float_loops_at_30_dbhz(self):
        fs = 4e6
        sat = dict(prn=7, code_phase=311.4, doppler=1234.0, cn0=30, phase=0.7)
        x = g.synthesize(fs, int(fs * 4.0), [sat], seed=3, noise_rms=8)
        a = g.acquire(x[:int(fs * 0.2)], fs, 7, blocks=200, doppler_step=250)
        out = {}
        for integer in (False, True):
            d = tm.track(x, fs, 7, a["epoch"], a["doppler"], integer=integer)
            edge = d[-1]["edge"]
            out[integer] = (tm.phase_lock(d[2300:], 20, (edge - 2300) % 20), d[-1]["doppler"], edge)
        self.assertEqual(out[True][2], out[False][2])
        self.assertGreater(out[True][0], 0.75)
        self.assertLess(abs(out[True][0] - out[False][0]), 0.05)
        self.assertLess(abs(out[True][1] - sat["doppler"]), 3.0)


if __name__ == "__main__":
    unittest.main()
