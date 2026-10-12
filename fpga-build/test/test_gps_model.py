#
# SPDX-License-Identifier: BSD-3-Clause

"""Checks of the GPS reference model: the C/A codes against the spec three ways, the Gold-code
correlation property, and acquisition of synthetic satellites in 8-bit noise."""

import unittest

import numpy as np

from dsp import gps_model as g


class TestCodes(unittest.TestCase):
    def test_gps_codes_three_ways(self):
        for prn in g.GPS_PRNS:
            code = g.ca_code(prn)
            self.assertTrue((code == g.ca_code_by_taps(prn)).all(), f"PRN {prn}: G2 stage pair")
            self.assertTrue((code == g.ca_code_by_delay(prn)).all(), f"PRN {prn}: G2 delay")
            first = int("".join(str(c) for c in code[:10]), 2)
            self.assertEqual(first, g.FIRST_CHIPS_OCTAL[prn], f"PRN {prn}: first ten chips")

    def test_sbas_codes_two_ways(self):
        for prn in g.SBAS_PRNS:
            self.assertEqual(g.g2_start(prn), g.SBAS_G2_INIT_OCTAL[prn], f"PRN {prn}: initial G2")
            self.assertTrue((g.ca_code(prn) == g.ca_code_by_delay(prn)).all(), f"PRN {prn}: delay")

    def test_gold_correlation_values(self):
        # Every periodic auto- (off peak) and cross-correlation of the 51 codes is -1, -65 or 63.
        codes = [1 - 2 * g.ca_code(p).astype(np.int64) for p in g.ALL_PRNS]
        spectra = [np.fft.fft(c) for c in codes]
        for i in range(len(codes)):
            for j in range(i, len(codes)):
                corr = np.rint(np.fft.ifft(spectra[i] * np.conj(spectra[j])).real).astype(np.int64)
                if i == j:
                    self.assertEqual(corr[0], 1023)
                    corr = corr[1:]
                self.assertTrue(set(np.unique(corr)) <= {-1, -65, 63},
                                f"PRNs {g.ALL_PRNS[i]}, {g.ALL_PRNS[j]}")


class TestAcquisition(unittest.TestCase):
    FS = 4e6

    def test_finds_satellites_and_rejects_absent_ones(self):
        sats = [
            dict(prn=3,  code_phase=12.25,  doppler=1750.0,  cn0=45),
            dict(prn=17, code_phase=600.6,  doppler=-3240.0, cn0=42),
            dict(prn=29, code_phase=1001.9, doppler=420.0,   cn0=40),
        ]
        x = g.synthesize(self.FS, int(self.FS * 0.010), sats, seed=5)
        for s in sats:
            r = g.acquire(x, self.FS, s["prn"])
            self.assertGreaterEqual(r["metric"], g.detect_metric(10), f"PRN {s['prn']}")
            err = (r["code_phase"] - s["code_phase"] + g.CODE_LENGTH / 2) % g.CODE_LENGTH \
                - g.CODE_LENGTH / 2
            self.assertLess(abs(err), 0.5, f"PRN {s['prn']} code phase")
            self.assertLessEqual(abs(r["doppler"] - s["doppler"]), 500.0, f"PRN {s['prn']} Doppler")
        for prn in (1, 8, 22):
            r = g.acquire(x, self.FS, prn)
            self.assertLess(r["metric"], g.detect_metric(10), f"PRN {prn} is absent")

    def test_finds_sbas_with_2_ms_symbols(self):
        x = g.synthesize(self.FS, int(self.FS * 0.010),
                         [dict(prn=131, code_phase=77.0, doppler=-150.0, cn0=43)], seed=7)
        r = g.acquire(x, self.FS, 131)
        self.assertGreaterEqual(r["metric"], g.detect_metric(10))
        self.assertLess(abs(r["code_phase"] - 77.0), 0.5)
        self.assertLess(g.acquire(x, self.FS, 133)["metric"], g.detect_metric(10))

    def test_epoch_is_the_first_sample_of_chip_0(self):
        x = g.synthesize(self.FS, int(self.FS * 0.010),
                         [dict(prn=9, code_phase=700.0, doppler=0.0, cn0=50)], seed=6)
        r = g.acquire(x, self.FS, 9)
        expected = (g.CODE_LENGTH - 700.0) * self.FS / g.CHIP_RATE
        self.assertLessEqual(abs(r["epoch"] - expected), 1.0)


if __name__ == "__main__":
    unittest.main()
