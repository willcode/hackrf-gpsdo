#
# SPDX-License-Identifier: BSD-3-Clause

"""Checks of the host-side decoders: the SBAS Viterbi decoder and CRC-24Q, and the TOW sequence filter;
and freq.py's steering loop (the app's reference). Run from the repository's root as
python -m unittest gps.test_gps_host, or from this directory as python -m unittest test_gps_host."""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import freq  # noqa: E402
import nav  # noqa: E402
import sbas  # noqa: E402


class TestSbas(unittest.TestCase):
    def test_crc24q_check_value(self):
        bits = [int(c) for byte in b"123456789" for c in format(byte, "08b")]
        self.assertEqual(sbas.crc24q(bits), 0xCDE703)

    def test_viterbi_corrects_noise(self):
        rng = np.random.default_rng(1)
        u = rng.integers(0, 2, 600)
        s, sym = 0, []
        for b in u:
            o = sbas._OUT[s, b]
            sym += [1 - 2 * int(o[0]), 1 - 2 * int(o[1])]
            s = (int(b) << 5) | (s >> 1)
        noisy = np.array(sym, float) + rng.normal(0, 0.7, len(sym))
        self.assertEqual(int((sbas.viterbi(noisy)[:590] != u[:590]).sum()), 0)


class TestNav(unittest.TestCase):
    def test_consistent_drops_a_chance_header(self):
        subs = [dict(bit=10 + 300 * k, tow=500 + k) for k in range(5)]
        subs.insert(2, dict(bit=377, tow=99999))
        kept = nav.consistent(subs)
        self.assertEqual([s["tow"] for s in kept], [500, 501, 502, 503, 504])


class TestServo(unittest.TestCase):
    @staticmethod
    def run_servo(estimates, start=0, step=4000000):
        s, end, out = freq.Servo(), start, []
        for e in estimates:
            if e is None:  # a second with no estimate
                end += step
                continue
            end += step
            s(end, e)
            out.append(s.target)
        return out

    def test_crosses_nominal(self):
        """A TCXO from -40 to +120 ppb: the steering goes from +40 across the gap at nominal (no
        ratio within 29.8 ppb) to -120; from the gap's edge the step limit once held it there."""
        out = self.run_servo([-40.0] * 20 + [120.0] * 60)
        self.assertAlmostEqual(out[19], 40.0, delta=0.01)
        self.assertAlmostEqual(out[-1], -120.0, delta=0.01)
        out = self.run_servo([-20.0] * 20 + [-5.0] * 20)
        self.assertEqual(out[-1], 0.0)                    # within half the gap: nominal

    def test_starts_over_after_a_pause(self):
        """After a pause longer than the history the line starts over from the new estimates
        alone: through the old ones and the new it overshot to 67.8 ppb (and past 537 s the
        pause's sample count wraps)."""
        out = self.run_servo([0.0] * 10 + [None] * 30 + [-60.0] * 20)
        self.assertLessEqual(max(out), 60.01)
        self.assertAlmostEqual(out[-1], 60.0, delta=0.01)


if __name__ == "__main__":
    unittest.main()
