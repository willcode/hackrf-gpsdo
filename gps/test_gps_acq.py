#!/usr/bin/env python3
#
# SPDX-License-Identifier: BSD-3-Clause

"""The M4's acquisition (gps/app/gps_acq.c) on the host, against the models and a recording.

Checked: every PRN's code and arm values against gps_model and gps_track_model; the FFT against
numpy; a block's box-car bins (read from a ring) and its mixing against numpy's; and on data/l1-45s.iq, a search of the GPS
PRNs from 32 one-millisecond blocks 10 ms apart (as the app takes snapshots) finds the satellites
the model's 32-block search finds strongly enough to survive the binning's worst loss (2.5 dB), at
the model's Doppler and code phase, and nothing the model does not.
"""

import ctypes
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fpga-build"))

from dsp import gps_model as g, gps_track_model as tm  # noqa: E402

N, BLOCK, FS = 1024, 4000, 4e6
REC = os.path.join(HERE, "data", "l1-45s.iq")



class TestGPSAcq(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        lib = os.path.join(cls.tmp, "libacq.so")
        a = os.path.join(HERE, "app")
        subprocess.run(["cc", "-std=c11", "-O2", "-Wall", "-Wextra", "-Werror", "-shared", "-fPIC",
                        "-I", a, "-o", lib, os.path.join(a, "gps_acq.c"),
                        os.path.join(a, "gps_nav.c"), "-lm"], check=True)
        L = cls.L = ctypes.CDLL(lib)
        F = ctypes.POINTER(ctypes.c_float)
        L.gps_acq_code.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_uint8)] + \
            [ctypes.POINTER(ctypes.c_uint32)] * 3
        L.gps_acq_fft.argtypes = [F, ctypes.c_int]
        I16 = ctypes.POINTER(ctypes.c_int16)
        L.gps_acq_code_spectrum.argtypes = [ctypes.POINTER(ctypes.c_uint8), F, I16]
        L.gps_acq_fft_stage.argtypes = [F, ctypes.c_int, ctypes.c_int]
        L.gps_acq_block.argtypes = [ctypes.c_void_p, ctypes.c_uint32, ctypes.c_uint32,
                                    ctypes.c_float, F]
        L.gps_acq_correlate.argtypes = [F, I16, F, F, F, ctypes.c_int, ctypes.c_int]
        L.gps_acq_shift.argtypes = [ctypes.c_uint32, ctypes.c_float]
        L.gps_acq_skip.argtypes = [ctypes.c_uint32, ctypes.c_float]
        L.gps_acq_skip.restype = ctypes.c_uint32
        L.gps_acq_product.argtypes = [F, I16, F, ctypes.c_int, ctypes.c_int]
        L.gps_acq_accumulate.argtypes = [F, F, F, ctypes.c_int, ctypes.c_int]
        L.gps_acq_peak.argtypes = [F, F, ctypes.POINTER(ctypes.c_int)]
        L.gps_acq_peak.restype = ctypes.c_float
        L.gps_acq_init()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    @staticmethod
    def fp(a):
        return a.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    @staticmethod
    def ip(a):
        return a.ctypes.data_as(ctypes.POINTER(ctypes.c_int16))

    def code(self, prn):
        c = (ctypes.c_uint8 * 1023)()
        v = [ctypes.c_uint32() for _ in range(3)]
        self.assertEqual(self.L.gps_acq_code(prn, c, *[ctypes.byref(x) for x in v]), 0)
        return np.array(c), [x.value for x in v]

    def test_codes_and_arm_values(self):
        for prn in g.ALL_PRNS:
            code, (g1, g2, win) = self.code(prn)
            np.testing.assert_array_equal(code, np.array(g.ca_code(prn)) & 1, err_msg=f"PRN {prn}")
            want = tm.arm_values(prn)
            self.assertEqual((g1, g2, win), tuple(int(v) for v in want), f"PRN {prn}")
        self.assertEqual(self.L.gps_acq_code(33, (ctypes.c_uint8 * 1023)(), *[
            ctypes.byref(ctypes.c_uint32()) for _ in range(3)]), -1)

    def test_fft(self):
        rng = np.random.default_rng(1)
        z = rng.standard_normal(N) + 1j * rng.standard_normal(N)
        for inverse, ref in ((0, np.fft.fft(z)), (1, np.fft.ifft(z) * N)):
            x = np.empty(2 * N, dtype=np.float32)
            x[0::2], x[1::2] = z.real, z.imag
            self.L.gps_acq_fft(self.fp(x), inverse)
            got = x[0::2] + 1j * x[1::2]
            self.assertLess(np.max(np.abs(got - ref)) / np.max(np.abs(ref)), 1e-4)

    @staticmethod
    def ring(words):
        """words as the app's ring: a power of two in size."""
        ring = np.zeros(1 << (len(words) - 1).bit_length(), dtype=np.uint16)
        ring[:len(words)] = words
        return ring.view(np.uint8)

    def block(self, ring, doppler, at=0):
        """The block from word at of ring, binned and mixed, in the stages' order."""
        out = np.empty(2 * N, dtype=np.float32)
        self.L.gps_acq_block(ring.ctypes.data, ring.nbytes, 2 * at, doppler, self.fp(out))
        return out

    def stages(self, x, inverse):
        for s in range(5):
            self.L.gps_acq_fft_stage(self.fp(x), inverse, s)

    def correlate(self, spec, cs, work, acc, acc1, shift):
        for half in (0, 1):
            self.L.gps_acq_correlate(self.fp(spec), self.ip(cs), self.fp(work), self.fp(acc),
                                     self.fp(acc1), shift, half)

    def test_bins_and_mix(self):
        rng = np.random.default_rng(2)
        i, q = rng.integers(-128, 128, 8192), rng.integers(-128, 128, 8192)
        words = ((i & 0xFF) | ((q & 0xFF) << 8)).astype(np.uint16)
        k = np.arange(BLOCK)
        b = k * N // BLOCK
        ring = self.ring(words)
        # Bin j lands at j's base-4 digits reversed.
        rev = np.zeros(N, dtype=int)
        for d in range(5):
            rev = rev * 4 + (np.arange(N) >> (2 * d)) % 4
        # 6000: across the ring's end; 6001: from a word's middle.
        for at in (0, 6000, 6001, 8191):
            idx = (at + k) % 8192
            ref = np.bincount(b, weights=i[idx], minlength=N) + \
                1j * np.bincount(b, weights=q[idx], minlength=N)
            for fd in (0.0, -3125.0, 4870.3):
                out = self.block(ring, fd, at)
                got = (out[0::2] + 1j * out[1::2])[rev]
                if fd == 0:
                    np.testing.assert_array_equal(got, ref)
                t = (np.arange(N) + 0.5) * BLOCK / N / FS
                want = ref * np.exp(-2j * np.pi * fd * t)
                self.assertLess(np.max(np.abs(got - want)) / np.max(np.abs(want)), 1e-4)

    def test_half_bins(self):
        """A code at 100, 100.25, 100.5 and 100.75 bins from the block's start: the best cell is
        the nearest half bin, and the half-bin row keeps the metric up."""
        code, _ = self.code(7)
        cs = np.empty(2 * N, dtype=np.int16)
        self.L.gps_acq_code_spectrum(code.astype(np.uint8).ctypes.data_as(
            ctypes.POINTER(ctypes.c_uint8)), self.fp(np.empty(2 * N, dtype=np.float32)), self.ip(cs))
        rng = np.random.default_rng(3)
        for frac, want in ((0.0, 200), (0.25, 200), (0.5, 201), (0.75, 202)):
            epoch = (100 + frac) * BLOCK / N                       # samples
            x = 40 * np.array(g.code_samples(7, FS, BLOCK, code_phase=-epoch * g.CHIP_RATE / FS))
            x = x + 20 * (rng.standard_normal(BLOCK) + 1j * rng.standard_normal(BLOCK))
            i, q = np.clip(np.round(x.real), -127, 127).astype(int), \
                np.clip(np.round(x.imag), -127, 127).astype(int)
            words = ((i & 0xFF) | ((q & 0xFF) << 8)).astype(np.uint16)
            work = np.empty(2 * N, dtype=np.float32)
            acc, acc1 = np.zeros(N, dtype=np.float32), np.zeros(N, dtype=np.float32)
            spec = self.block(self.ring(words), 0.0)
            self.stages(spec, 0)
            self.correlate(spec, cs, work, acc, acc1, 0)
            cell = ctypes.c_int()
            m = self.L.gps_acq_peak(self.fp(acc), self.fp(acc1), ctypes.byref(cell))
            self.assertIn(cell.value, (want, want + (1 if frac == 0.25 else 0)), f"{frac}")
            m0 = self.L.gps_acq_peak(self.fp(acc), None, ctypes.byref(ctypes.c_int()))
            if frac == 0.5:
                self.assertGreater(m, 1.5 * m0)                    # the half bin recovers it

    def code_spectrum(self, prn):
        code, _ = self.code(prn)
        cs = np.empty(2 * N, dtype=np.int16)
        self.L.gps_acq_code_spectrum(code.astype(np.uint8).ctypes.data_as(
            ctypes.POINTER(ctypes.c_uint8)), self.fp(np.empty(2 * N, dtype=np.float32)), self.ip(cs))
        return cs

    def test_rotation(self):
        """The spectrum mixed at a Doppler and taken a bin on correlates as the block mixed 1 kHz
        on: the cold pass's pairs share a block's bins and FFT."""
        rng = np.random.default_rng(4)
        i, q = rng.integers(-60, 60, 8192), rng.integers(-60, 60, 8192)
        ring = self.ring(((i & 0xFF) | ((q & 0xFF) << 8)).astype(np.uint16))
        cs = self.code_spectrum(11)
        for fd in (-4500.0, 1250.0):
            rows = []
            for f, rot in ((fd, 1), (fd + 1000, 0)):
                spec = self.block(ring, f)
                self.stages(spec, 0)
                work = np.empty(2 * N, dtype=np.float32)
                self.L.gps_acq_product(self.fp(spec), self.ip(cs), self.fp(work), 0, rot)
                self.stages(work, 1)
                rows.append(work[0::2] ** 2 + work[1::2] ** 2)
            self.assertLess(np.max(np.abs(rows[0] - rows[1])) / np.max(rows[1]), 1e-3)

    def test_skip_to_whole_bins(self):
        """A block started gps_acq_skip samples on starts within half a sample of the first's
        bins, so its shift is even; with no half-bin row an odd shift goes to the next bin."""
        rng = np.random.default_rng(5)
        for _ in range(2000):
            since, fd = int(rng.integers(0, 600000)), float(rng.uniform(-5000, 5000))
            d = self.L.gps_acq_skip(since, fd)
            self.assertLess(d, 4)
            p = BLOCK / (1 + fd / 1575.42e6)
            r = (since + d) % p
            e = (r * N / BLOCK) % 1
            self.assertLess(min(e, 1 - e) * BLOCK / N, 0.5)
            self.assertEqual(self.L.gps_acq_shift(since + d, fd) % 2, 0)
        work = np.zeros(2 * N, dtype=np.float32)
        work[2 * 10] = 3.0
        acc = np.zeros(N, dtype=np.float32)
        self.L.gps_acq_accumulate(self.fp(work), self.fp(acc), None, 7, 0)
        self.assertEqual(int(np.argmax(acc)), 14)
        self.assertEqual(acc[14], 9.0)

    @unittest.skipUnless(os.path.exists(REC), "no recording in data/")
    def test_search_recording(self):
        raw = np.fromfile(REC, dtype=np.int8, count=2 * BLOCK * 330, offset=2 * BLOCK * 100)
        i, q = raw[0::2].astype(np.int64), raw[1::2].astype(np.int64)
        words = ((i & 0xFF) | ((q & 0xFF) << 8)).astype(np.uint16)
        starts = [40000 * k for k in range(32)]     # 32 blocks, 10 ms apart
        dopplers = np.arange(-5000, 5001, 250.0)
        code_spec = np.empty(2 * N, dtype=np.int16)
        work = np.empty(2 * N, dtype=np.float32)
        rings = [self.ring(words[s:s + BLOCK]) for s in starts]
        found = {}
        for prn in range(1, 33):
            code, _ = self.code(prn)
            self.L.gps_acq_code_spectrum(code.astype(np.uint8).ctypes.data_as(
                ctypes.POINTER(ctypes.c_uint8)), self.fp(work), self.ip(code_spec))
            best = None
            for fd in dopplers:
                acc = np.zeros(N, dtype=np.float32)
                acc1 = np.zeros(N, dtype=np.float32)
                for s, ring in zip(starts, rings):
                    spec = self.block(ring, fd)
                    self.stages(spec, 0)
                    self.correlate(spec, code_spec, work, acc, acc1, self.L.gps_acq_shift(s, fd))
                cell = ctypes.c_int()
                m = self.L.gps_acq_peak(self.fp(acc), self.fp(acc1), ctypes.byref(cell))
                if best is None or m > best[0]:
                    best = (m, fd, cell.value)
            found[prn] = best
        # The model's search on 32 contiguous blocks from the same start.
        x = (i + 1j * q)[:32 * BLOCK]
        ref = {f["prn"]: f for f in g.acquire_all(x, FS, prns=range(1, 33), blocks=32,
                                                  doppler_max=5000, doppler_step=250)
               if f["detected"]}
        self.assertGreaterEqual(len(ref), 4)
        threshold = 1 + 10 / np.sqrt(len(starts))
        # A code a quarter bin from both grids keeps (1 - 1/4)^2 of its excess over the noise.
        required = [p for p, f in ref.items() if 1 + 0.5625 * (f["metric"] - 1) > threshold]
        self.assertGreaterEqual(len(required), 3)
        for prn in required:
            m, fd, cell = found[prn]
            f = ref[prn]
            self.assertGreater(m, threshold, f"PRN {prn}")
            self.assertLessEqual(abs(fd - f["doppler"]), 250, f"PRN {prn}")
            want = 2 * f["epoch"] * N / BLOCK                     # in half bins
            d = (cell - want + N) % (2 * N) - N
            self.assertLessEqual(abs(d), 1.5, f"PRN {prn}: half bin {cell}, model {want:.1f}")
        # Nothing over the threshold that the model does not detect, but a strong satellite's
        # cross-correlation.
        top = max(v[0] for v in found.values())
        for prn, (m, fd, cell) in found.items():
            if prn not in ref:
                self.assertTrue(m < threshold or m - 1 < 0.02 * (top - 1), f"PRN {prn}: {m:.2f}")


if __name__ == "__main__":
    unittest.main()
