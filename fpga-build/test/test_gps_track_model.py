#
# SPDX-License-Identifier: BSD-3-Clause

"""The channel model and its loops on a synthetic satellite: from an acquisition's coarse code epoch and
Doppler bin, the code must converge to the true phase and the carrier must lock."""

import unittest

import numpy as np

from dsp import gps_model as g, gps_track_model as tm


class TestTracking(unittest.TestCase):
    FS = 4e6

    def track(self, cn0, seconds, seed):
        sat = dict(prn=7, code_phase=311.4, doppler=1234.0, cn0=cn0, phase=0.7)
        x = g.synthesize(self.FS, int(self.FS * seconds), [sat], seed=seed, noise_rms=8)
        a = g.acquire(x[:int(self.FS * 0.2)], self.FS, 7, blocks=200, doppler_step=250)
        return sat, tm.track(x, self.FS, 7, a["epoch"], a["doppler"])

    def code_errors(self, sat, dumps):
        rate = g.CHIP_RATE * (1 + sat["doppler"] / g.L1_HZ)
        truth = [(sat["code_phase"] + d["n0"] / self.FS * rate) % g.CODE_LENGTH for d in dumps]
        return np.array([(d["F"] / tm.NCO_ONE - t + 511.5) % g.CODE_LENGTH - 511.5
                         for d, t in zip(dumps, truth)])

    def test_strong_signal(self):
        sat, dumps = self.track(45, 2.0, 3)
        late = dumps[1500:]
        self.assertLess(np.max(np.abs(self.code_errors(sat, late))), 0.05)
        self.assertLess(abs(dumps[-1]["doppler"] - sat["doppler"]), 2.0)
        edge = dumps[-1]["edge"]
        self.assertGreater(tm.phase_lock(late, 20, (edge - 1500) % 20), 0.95)
        self.assertLess(abs(tm.cn0_dbhz(dumps[300:]) - sat["cn0"]), 1.5)

    def test_weak_signal_holds_phase_lock_on_symbol_sums(self):
        # 30 dB-Hz, the weakest satellites of the first recordings. Seed 3 put the earlier loops
        # (1 ms PLL, or a 1 ms-pair FLL assist) off lock for seconds.
        sat, dumps = self.track(30, 4.0, 3)
        late = dumps[2300:]
        edge = dumps[-1]["edge"]
        self.assertGreater(tm.phase_lock(late, 20, (edge - 2300) % 20), 0.75)
        self.assertLess(abs(dumps[-1]["doppler"] - sat["doppler"]), 3.0)
        self.assertLess(np.std(self.code_errors(sat, late)), 0.05)


if __name__ == "__main__":
    unittest.main()
