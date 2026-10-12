#
# SPDX-License-Identifier: BSD-3-Clause

"""Reference model of GPS and SBAS L1 C/A reception, for the correlator image to be checked against.

This floating-point reference precedes the gateware's integer arithmetic. It holds:

1. The C/A codes: GPS PRN 1-32 (IS-GPS-200) and SBAS PRN 120-138 (RTCA DO-229). Each is G1 xor G2 of
   two 10-stage generators, G1 = 1 + x^3 + x^10 and G2 = 1 + x^2 + x^3 + x^6 + x^8 + x^9 + x^10, G1
   starting at all ones. A PRN differs only in where G2 starts, so one generator with a 10-bit G2 start
   value makes every code: that is how the gateware makes them (`ca_code`). The spec's other
   definitions are kept so the tests can check it: for GPS, the G2 stage pairs, the G2 delays and the
   first ten chips; for SBAS, the G2 delays and the initial G2 settings.
2. A signal synthesizer: complex baseband at the sample rate, any number of satellites with code
   phase, Doppler (carrier and code), C/N0 and data symbols (GPS: 20 code periods a symbol; SBAS: 2),
   plus noise, quantized to the 8-bit I/Q the HackRF delivers.
3. Acquisition by parallel code-phase search: per Doppler bin, 1 ms blocks with the carrier wiped off
   are correlated with the code through FFTs, and the squared magnitudes are summed over blocks.

Conventions. Chip c(k) is 0 or 1 and is sent as 1 - 2 c(k). A satellite's code phase is the chip,
fractional, being sent at sample 0. Acquisition reports the code epoch: the first sample at which
chip 0 begins, in [0, samples per ms).
"""

from functools import lru_cache

import numpy as np

L1_HZ        = 1575.42e6
CHIP_RATE    = 1.023e6
CODE_LENGTH  = 1023
GPS_PRNS     = range(1, 33)
SBAS_PRNS    = range(120, 139)
ALL_PRNS     = list(GPS_PRNS) + list(SBAS_PRNS)


def codes_per_symbol(prn):
    """Code periods per data symbol: GPS 50 bit/s; SBAS 500 symbol/s (250 bit/s, rate 1/2 coded)."""
    return 2 if prn in SBAS_PRNS else 20


# IS-GPS-200 Table 3-Ia: per PRN, the two G2 stages (1-based) summed for the code, the G2 delay in
# chips that the same sum equals, and the first ten chips in octal (the first chip is the leading 1).
G2_TAPS = {
     1: (2, 6),   2: (3, 7),   3: (4, 8),   4: (5, 9),   5: (1, 9),   6: (2, 10),  7: (1, 8),
     8: (2, 9),   9: (3, 10), 10: (2, 3),  11: (3, 4),  12: (5, 6),  13: (6, 7),  14: (7, 8),
    15: (8, 9),  16: (9, 10), 17: (1, 4),  18: (2, 5),  19: (3, 6),  20: (4, 7),  21: (5, 8),
    22: (6, 9),  23: (1, 3),  24: (4, 6),  25: (5, 7),  26: (6, 8),  27: (7, 9),  28: (8, 10),
    29: (1, 6),  30: (2, 7),  31: (3, 8),  32: (4, 9),
}
# RTCA DO-229 for PRN 120-138, which also agrees with gnss-sdr's table.
G2_DELAY = {
     1: 5,    2: 6,    3: 7,    4: 8,    5: 17,   6: 18,   7: 139,  8: 140,  9: 141, 10: 251,
    11: 252, 12: 254, 13: 255, 14: 256, 15: 257, 16: 258, 17: 469, 18: 470, 19: 471, 20: 472,
    21: 473, 22: 474, 23: 509, 24: 512, 25: 513, 26: 514, 27: 515, 28: 516, 29: 859, 30: 860,
    31: 861, 32: 862,
   120: 145, 121: 175, 122: 52,  123: 21,  124: 237, 125: 235, 126: 886, 127: 657, 128: 634,
   129: 762, 130: 355, 131: 1012, 132: 176, 133: 603, 134: 130, 135: 359, 136: 595, 137: 68,
   138: 386,
}
FIRST_CHIPS_OCTAL = {
     1: 0o1440,  2: 0o1620,  3: 0o1710,  4: 0o1744,  5: 0o1133,  6: 0o1455,  7: 0o1131,
     8: 0o1454,  9: 0o1626, 10: 0o1504, 11: 0o1642, 12: 0o1750, 13: 0o1764, 14: 0o1772,
    15: 0o1775, 16: 0o1776, 17: 0o1156, 18: 0o1467, 19: 0o1633, 20: 0o1715, 21: 0o1746,
    22: 0o1763, 23: 0o1063, 24: 0o1706, 25: 0o1743, 26: 0o1761, 27: 0o1770, 28: 0o1774,
    29: 0o1127, 30: 0o1453, 31: 0o1625, 32: 0o1712,
}
# RTCA DO-229: the initial G2 setting of each SBAS PRN in octal, the first G2 output chip leading.
SBAS_G2_INIT_OCTAL = {
    120: 0o1106, 121: 0o1241, 122: 0o0267, 123: 0o0232, 124: 0o1617, 125: 0o1076, 126: 0o1764,
    127: 0o0717, 128: 0o1532, 129: 0o1250, 130: 0o0341, 131: 0o0551, 132: 0o0520, 133: 0o1731,
    134: 0o0706, 135: 0o1216, 136: 0o0740, 137: 0o1007, 138: 0o0450,
}

G1_TAPS = (3, 10)
G2_FEEDBACK = (2, 3, 6, 8, 9, 10)


def lfsr_out(taps, start=(1,) * 10, n=CODE_LENGTH):
    """n outputs of a 10-stage generator: its stage 10, then a shift with stage 1 = xor of the taps.

    start[i] is stage i + 1 at the first step.
    """
    reg = list(start)
    out = np.empty(n, dtype=np.uint8)
    for k in range(n):
        out[k] = reg[9]
        f = 0
        for t in taps:
            f ^= reg[t - 1]
        reg = [f] + reg[:9]
    return out


def g2_start(prn):
    """The PRN's 10-bit G2 start value, as the gateware takes it: bit 9 is stage 10, the first output.

    G2's first ten outputs are stages 10, 9, ..., 1 at the start, so the start value is G2's first ten
    chips from the PRN's delay, read with the first chip as the most significant bit.
    """
    first = np.roll(lfsr_out(G2_FEEDBACK), G2_DELAY[prn])[:10]
    return int("".join(str(c) for c in first), 2)


def ca_code(prn):
    """The 1023 chips (0 or 1) of the PRN's code, as the gateware makes them: G1 from all ones, G2
    from the PRN's start value."""
    return _ca_code(prn).copy()


@lru_cache(maxsize=None)
def _ca_code(prn):
    v = g2_start(prn)
    start = [(v >> i) & 1 for i in range(10)]      # stage i + 1 = bit i
    return lfsr_out(G1_TAPS) ^ lfsr_out(G2_FEEDBACK, start=start)


def ca_code_by_taps(prn):
    """A GPS PRN's code as IS-GPS-200 draws it: G1 xor the sum of two G2 stages, both from all ones."""
    a, b = G2_TAPS[prn]
    g1 = [1] * 10
    g2 = [1] * 10
    chips = np.empty(CODE_LENGTH, dtype=np.uint8)
    for k in range(CODE_LENGTH):
        chips[k] = g1[9] ^ g2[a - 1] ^ g2[b - 1]
        f1 = g1[2] ^ g1[9]
        f2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9]
        g1 = [f1] + g1[:9]
        g2 = [f2] + g2[:9]
    return chips


def ca_code_by_delay(prn):
    """The code as G1 xor G2 delayed by the PRN's G2 delay."""
    return lfsr_out(G1_TAPS) ^ np.roll(lfsr_out(G2_FEEDBACK), G2_DELAY[prn])


def code_samples(prn, fs, n, code_phase=0.0, code_rate=CHIP_RATE):
    """n samples of the code as +1/-1, chip code_phase at sample 0, at code_rate chips per second."""
    chips = 1.0 - 2.0 * ca_code(prn)
    idx = np.floor(code_phase + np.arange(n) * (code_rate / fs)).astype(np.int64) % CODE_LENGTH
    return chips[idx]


def synthesize(fs, n, sats, noise_rms=20.0, seed=0, bits=8):
    """n complex samples quantized to signed `bits`-bit I and Q, as a complex128 array.

    `sats` is a list of dicts: prn, code_phase (chips at sample 0), doppler (Hz), cn0 (dB-Hz), and
    optionally phase (rad) and data (array of +1/-1 symbols, each lasting codes_per_symbol(prn) code
    periods; random when missing). The code is
    Doppler-stretched by doppler / 1540, the ratio of L1 to the chip rate. noise_rms is per component,
    in LSB.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(n) / fs
    noise_power = 2 * noise_rms ** 2                    # complex, over the band fs
    x = noise_rms * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    for s in sats:
        amp = np.sqrt(10 ** (s["cn0"] / 10) * noise_power / fs)
        code_rate = CHIP_RATE * (1 + s["doppler"] / L1_HZ)
        chip_pos = s["code_phase"] + t * code_rate      # chips since the code's chip 0, unwrapped
        chips = 1.0 - 2.0 * ca_code(s["prn"])
        code = chips[np.floor(chip_pos).astype(np.int64) % CODE_LENGTH]
        bit_index = np.floor(chip_pos / (CODE_LENGTH * codes_per_symbol(s["prn"]))).astype(np.int64)
        data = s.get("data")
        if data is None:
            data = rng.choice([-1.0, 1.0], size=bit_index.max() + 2)
        bit = np.asarray(data)[bit_index - bit_index.min()]
        phase = 2 * np.pi * s["doppler"] * t + s.get("phase", 0.0)
        x += amp * bit * code * np.exp(1j * phase)
    lim = 1 << (bits - 1)
    i = np.clip(np.round(x.real), -lim, lim - 1)
    q = np.clip(np.round(x.imag), -lim, lim - 1)
    return i + 1j * q


def read_hackrf_iq(path, count=None, offset=0):
    """Interleaved signed 8-bit I/Q (hackrf_transfer -r) as complex128."""
    raw = np.fromfile(path, dtype=np.int8, count=-1 if count is None else 2 * count,
                      offset=2 * offset)
    return raw[0::2].astype(np.float64) + 1j * raw[1::2].astype(np.float64)


class RawIQ:
    """A hackrf_transfer file read lazily: len() in samples, a slice gives complex128."""

    def __init__(self, path):
        self.raw = np.memmap(path, dtype=np.int8, mode="r")

    def __len__(self):
        return len(self.raw) // 2

    def __getitem__(self, sl):
        start, stop, _ = sl.indices(len(self))
        r = np.asarray(self.raw[2 * start:2 * stop], dtype=np.float64)
        return r[0::2] + 1j * r[1::2]


def acquire(x, fs, prn, doppler_max=6000.0, doppler_step=500.0, blocks=10):
    """Parallel code-phase search for one PRN over `blocks` 1 ms blocks summed noncoherently.

    Returns a dict: epoch (first sample of chip 0, in [0, samples per ms)), code_phase (chip being sent
    at sample 0), doppler (Hz, the bin's center), metric (the peak over the mean of the cells more than
    one chip from it) and the full search grid.
    """
    n = int(round(fs / 1000))
    if abs(n - fs / 1000) > 1e-9:
        raise ValueError("the sample rate must be a whole number of samples per millisecond")
    if len(x) < n * blocks:
        raise ValueError("not enough samples for the blocks asked for")
    local = np.conj(np.fft.fft(code_samples(prn, fs, n)))
    dopplers = np.arange(-doppler_max, doppler_max + doppler_step / 2, doppler_step)
    t = np.arange(n * blocks) / fs
    grid = np.zeros((len(dopplers), n))
    for di, fd in enumerate(dopplers):
        y = (x[:n * blocks] * np.exp(-2j * np.pi * fd * t)).reshape(blocks, n)
        corr = np.fft.ifft(np.fft.fft(y, axis=1) * local, axis=1)
        grid[di] = np.sum(np.abs(corr) ** 2, axis=0)
    di, epoch = np.unravel_index(np.argmax(grid), grid.shape)
    guard = int(np.ceil(fs / CHIP_RATE)) + 1
    dist = np.abs((np.arange(n) - epoch + n // 2) % n - n // 2)
    noise = grid[di][dist > guard].mean()
    code_phase = (-epoch * CHIP_RATE / fs) % CODE_LENGTH
    return dict(prn=prn, epoch=int(epoch), code_phase=code_phase, doppler=float(dopplers[di]),
                metric=float(grid[di, epoch] / noise), grid=grid, dopplers=dopplers)


# The metric of a noise cell averages 1 and its spread falls as 1 / sqrt(blocks), so a peak is a
# satellite above 1 + DETECT_SPREAD / sqrt(blocks). The largest of the 100,000 cells of a 4 Msps,
# +-6 kHz search in pure noise sat at 1 + 6.5 to 8.1 / sqrt(blocks) for 10, 40 and 80 blocks.
DETECT_SPREAD = 10.0


def detect_metric(blocks):
    return 1 + DETECT_SPREAD / np.sqrt(blocks)


# A weak peak is taken for a strong satellite's cross-correlation when its excess over the noise
# (metric - 1) is below this share of the strongest detection's. In a 200-block search beside a
# 44 dB-Hz satellite (metric 52), a dozen absent PRNs reached 1.7-1.9 at Dopplers spread around it,
# while real satellites in the same search stood at 3.3 and above.
CROSS_SHARE = 0.02


def acquire_all(x, fs, prns=ALL_PRNS, blocks=10, **kw):
    """acquire() for every PRN; the results sorted by metric, strongest first.

    detected: above detect_metric(blocks) and not explained as a cross-correlation of the strongest
    detection (`cross` marks those).
    """
    found = sorted((acquire(x, fs, prn, blocks=blocks, **kw) for prn in prns),
                   key=lambda f: -f["metric"])
    top = found[0]["metric"] if found else 1.0
    for f in found:
        above = f["metric"] >= detect_metric(blocks)
        f["cross"] = above and f is not found[0] and f["metric"] - 1 < CROSS_SHARE * (top - 1)
        f["detected"] = above and not f["cross"]
    return found
