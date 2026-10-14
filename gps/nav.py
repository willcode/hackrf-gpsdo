#
# SPDX-License-Identifier: BSD-3-Clause

"""GPS navigation from 1 ms prompt sums (IS-GPS-200 20.3): bit sync, bits, parity, subframes,
ephemeris and almanac, satellite positions, the position solve, DOP and geodetic conversion.

The reference the C (gps/app/gps_nav.c) is checked against.
"""

import numpy as np

PREAMBLE = (1, 0, 0, 0, 1, 0, 1, 1)
BIT_MS = 20


def bit_sync(prompt):
    """The ms offset (0-19) of bit edges in a complex prompt series, and the share of sign changes there.

    Differential: a data edge flips the product of two neighbors' prompts negative, whatever the
    carrier phase.
    """
    d = np.real(prompt[1:] * np.conj(prompt[:-1]))
    flips = np.zeros(BIT_MS)
    for k in np.nonzero(d < 0)[0]:
        flips[(k + 1) % BIT_MS] += -d[k]
    edge = int(np.argmax(flips))
    return edge, float(flips[edge] / flips.sum()) if flips.sum() else 0.0


def bits(prompt, edge, coherent=True):
    """One soft value per whole bit from the bit edge on.

    coherent: the real part of each bit's 20-sum (needs carrier phase lock). Otherwise differential: the
    real part of each bit's sum times the conjugate of the one before, whose sign says whether the bit
    changed.
    """
    n = (len(prompt) - edge) // BIT_MS
    s = prompt[edge:edge + n * BIT_MS].reshape(n, BIT_MS).sum(axis=1)
    if coherent:
        return np.real(s)
    return np.real(s[1:] * np.conj(s[:-1]))


def parity_ok(word, d29, d30):
    """IS-GPS-200 parity of a 30-bit word (list of 0/1, D1 first) given the last two bits before it.

    Returns the 24 data bits with the D30* inversion removed, or None.
    """
    d = [b ^ d30 for b in word[:24]]

    def x(*idx):
        v = 0
        for i in idx:
            v ^= d[i - 1]
        return v

    p = [
        d29 ^ x(1, 2, 3, 5, 6, 10, 11, 12, 13, 14, 17, 18, 20, 23),
        d30 ^ x(2, 3, 4, 6, 7, 11, 12, 13, 14, 15, 18, 19, 21, 24),
        d29 ^ x(1, 3, 4, 5, 7, 8, 12, 13, 14, 15, 16, 19, 20, 22),
        d30 ^ x(2, 4, 5, 6, 8, 9, 13, 14, 15, 16, 17, 20, 21, 23),
        d30 ^ x(1, 3, 5, 6, 7, 9, 10, 14, 15, 16, 17, 18, 21, 22, 24),
        d29 ^ x(3, 5, 6, 8, 9, 10, 11, 13, 15, 19, 22, 23, 24),
    ]
    return d if p == list(word[24:30]) else None


def subframes(hard):
    """Every subframe start in a 0/1 bit list whose TLM and HOW pass parity, in either polarity.

    Returns dicts: bit (index of the preamble's first bit), inverted, tow (HOW's 17-bit count: the GPS
    time of week of the next subframe's start, in 6 s units) and id (1-5).
    """
    found = []
    for inv in (0, 1):
        b = [v ^ inv for v in hard]
        for i in range(2, len(b) - 60 + 1):
            if tuple(b[i:i + 8]) != PREAMBLE:
                continue
            tlm = parity_ok(b[i:i + 30], b[i - 2], b[i - 1])
            if tlm is None:
                continue
            how = parity_ok(b[i + 30:i + 60], b[i + 28], b[i + 29])
            if how is None:
                continue
            tow = int("".join(map(str, how[:17])), 2)
            sid = int("".join(map(str, how[19:22])), 2)
            if 1 <= sid <= 5:
                found.append(dict(bit=i, inverted=bool(inv), tow=tow, id=sid))
    return sorted(found, key=lambda f: f["bit"])


# Ephemeris and the position solution (IS-GPS-200 20.3.3 and Table 20-IV).

C        = 299792458.0
MU       = 3.986005e14
OMEGA_E  = 7.2921151467e-5
GPS_PI   = 3.1415926535898
REL_F    = -4.442807633e-10
WEEK_S   = 604800
# The orbit a GPS satellite can have (gps_nav.h): sqrt(A) (m^0.5; 5,153.6 nominal) and eccentricity.
SQRTA_MIN, SQRTA_MAX, E_MAX = 5000.0, 5300.0, 0.03


def words(hard, start):
    """The ten words of the subframe at bit `start` as 24-bit lists, or None if any fails parity."""
    out = []
    for w in range(10):
        i = start + 30 * w
        if i + 30 > len(hard) or i < 2:
            return None
        d = parity_ok(hard[i:i + 30], hard[i - 2], hard[i - 1])
        if d is None:
            return None
        out.append(d)
    return out


def field(ws, *parts, signed=False):
    """A field made of (word, first bit, last bit) parts, most significant first; words and bits from 1."""
    v, n = 0, 0
    for w, a, b in parts:
        for bit in ws[w - 1][a - 1:b]:
            v = (v << 1) | bit
            n += 1
    if signed and v >> (n - 1):
        v -= 1 << n
    return v


def decode_subframe(ws):
    """The ID, the TOW count, and the clock or orbit fields of subframes 1-3, scaled to SI and radians."""
    sid = field(ws, (2, 20, 22))
    out = dict(id=sid, tow=field(ws, (2, 1, 17)))
    sc = GPS_PI
    if sid == 1:
        out.update(
            week=field(ws, (3, 1, 10)), health=field(ws, (3, 17, 22)),
            iodc=field(ws, (3, 23, 24), (8, 1, 8)),
            tgd=field(ws, (7, 17, 24), signed=True) * 2.0 ** -31,
            toc=field(ws, (8, 9, 24)) * 16.0,
            af2=field(ws, (9, 1, 8), signed=True) * 2.0 ** -55,
            af1=field(ws, (9, 9, 24), signed=True) * 2.0 ** -43,
            af0=field(ws, (10, 1, 22), signed=True) * 2.0 ** -31)
    elif sid == 2:
        out.update(
            iode=field(ws, (3, 1, 8)),
            crs=field(ws, (3, 9, 24), signed=True) * 2.0 ** -5,
            dn=field(ws, (4, 1, 16), signed=True) * 2.0 ** -43 * sc,
            m0=field(ws, (4, 17, 24), (5, 1, 24), signed=True) * 2.0 ** -31 * sc,
            cuc=field(ws, (6, 1, 16), signed=True) * 2.0 ** -29,
            e=field(ws, (6, 17, 24), (7, 1, 24)) * 2.0 ** -33,
            cus=field(ws, (8, 1, 16), signed=True) * 2.0 ** -29,
            sqrta=field(ws, (8, 17, 24), (9, 1, 24)) * 2.0 ** -19,
            toe=field(ws, (10, 1, 16)) * 16.0)
    elif sid == 3:
        out.update(
            cic=field(ws, (3, 1, 16), signed=True) * 2.0 ** -29,
            omega0=field(ws, (3, 17, 24), (4, 1, 24), signed=True) * 2.0 ** -31 * sc,
            cis=field(ws, (5, 1, 16), signed=True) * 2.0 ** -29,
            i0=field(ws, (5, 17, 24), (6, 1, 24), signed=True) * 2.0 ** -31 * sc,
            crc=field(ws, (7, 1, 16), signed=True) * 2.0 ** -5,
            omega=field(ws, (7, 17, 24), (8, 1, 24), signed=True) * 2.0 ** -31 * sc,
            omegadot=field(ws, (9, 1, 24), signed=True) * 2.0 ** -43 * sc,
            iode=field(ws, (10, 1, 8)),
            idot=field(ws, (10, 9, 22), signed=True) * 2.0 ** -43 * sc)
    return out


def plausible(orbit):
    """An ephemeris's or almanac entry's orbit within SQRTA_MIN-SQRTA_MAX and E_MAX, as gps_nav.c
    takes them: any parity-valid subframe can be sent on the air, and sqrt(A) = 0 leaves
    sat_position dividing by zero."""
    return SQRTA_MIN < orbit["sqrta"] < SQRTA_MAX and orbit["e"] < E_MAX


def decode_almanac(ws):
    """A subframe 4 or 5 page's almanac (IS-GPS-200 20.3.3.5.1.2, Table 20-VI) as an ephemeris for
    sat_position, with the PRN and health; None for a page that holds none (data ID not 1, or an SV
    ID outside 1-32: the pages of other data and the dummy SV 0). The terms an almanac lacks are 0,
    its reference time stands for toe and toc, and i0 is 0.3 semicircles plus its offset."""
    sid = field(ws, (2, 20, 22))
    prn = field(ws, (3, 3, 8))
    if sid not in (4, 5) or field(ws, (3, 1, 2)) != 1 or not 1 <= prn <= 32:
        return None
    sc = GPS_PI
    toa = field(ws, (4, 1, 8)) * 4096.0
    return dict(
        prn=prn, health=field(ws, (5, 17, 24)), toa=toa, toe=toa, toc=toa,
        e=field(ws, (3, 9, 24)) * 2.0 ** -21,
        i0=(0.3 + field(ws, (4, 9, 24), signed=True) * 2.0 ** -19) * sc,
        omegadot=field(ws, (5, 1, 16), signed=True) * 2.0 ** -38 * sc,
        sqrta=field(ws, (6, 1, 24)) * 2.0 ** -11,
        omega0=field(ws, (7, 1, 24), signed=True) * 2.0 ** -23 * sc,
        omega=field(ws, (8, 1, 24), signed=True) * 2.0 ** -23 * sc,
        m0=field(ws, (9, 1, 24), signed=True) * 2.0 ** -23 * sc,
        af0=field(ws, (10, 1, 8), (10, 20, 22), signed=True) * 2.0 ** -20,
        af1=field(ws, (10, 9, 19), signed=True) * 2.0 ** -38,
        af2=0.0, tgd=0.0, dn=0.0, idot=0.0, cuc=0.0, cus=0.0, crc=0.0, crs=0.0, cic=0.0, cis=0.0)


def _wrap(t):
    return (t + WEEK_S / 2) % WEEK_S - WEEK_S / 2


def sat_position(eph, t):
    """ECEF position (m) at GPS time of week t, and the clock correction (s), from one ephemeris."""
    a = eph["sqrta"] ** 2
    tk = _wrap(t - eph["toe"])
    n = np.sqrt(MU / a ** 3) + eph["dn"]
    m = eph["m0"] + n * tk
    e_anom = m
    for _ in range(10):
        e_anom = m + eph["e"] * np.sin(e_anom)
    nu = np.arctan2(np.sqrt(1 - eph["e"] ** 2) * np.sin(e_anom), np.cos(e_anom) - eph["e"])
    phi = nu + eph["omega"]
    s2, c2 = np.sin(2 * phi), np.cos(2 * phi)
    u = phi + eph["cus"] * s2 + eph["cuc"] * c2
    r = a * (1 - eph["e"] * np.cos(e_anom)) + eph["crs"] * s2 + eph["crc"] * c2
    i = eph["i0"] + eph["idot"] * tk + eph["cis"] * s2 + eph["cic"] * c2
    xp, yp = r * np.cos(u), r * np.sin(u)
    om = eph["omega0"] + (eph["omegadot"] - OMEGA_E) * tk - OMEGA_E * eph["toe"]
    pos = np.array([xp * np.cos(om) - yp * np.cos(i) * np.sin(om),
                    xp * np.sin(om) + yp * np.cos(i) * np.cos(om),
                    yp * np.sin(i)])
    dt = _wrap(t - eph["toc"])
    clock = (eph["af0"] + eph["af1"] * dt + eph["af2"] * dt * dt
             + REL_F * eph["e"] * eph["sqrta"] * np.sin(e_anom) - eph["tgd"])
    return pos, clock


def ecef_to_geodetic(p):
    """WGS 84 latitude and longitude (degrees) and height (m)."""
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    lon = np.arctan2(p[1], p[0])
    rho = np.hypot(p[0], p[1])
    lat = np.arctan2(p[2], rho * (1 - e2))
    for _ in range(10):
        n = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        h = rho / np.cos(lat) - n
        lat = np.arctan2(p[2], rho * (1 - e2 * n / (n + h)))
    return np.degrees(lat), np.degrees(lon), h


def elevation(rx, sat):
    lat, lon, _ = ecef_to_geodetic(rx)
    lat, lon = np.radians(lat), np.radians(lon)
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    d = sat - rx
    return np.arcsin(d @ up / np.linalg.norm(d))


def dop(rx, sats):
    """Geometric, position, horizontal, vertical and time dilution of precision at ECEF rx for
    satellites at the ECEF positions sats; infinite for a geometry that fixes nothing (fewer than
    four, or in a plane through rx)."""
    if len(sats) < 4:
        return (np.inf,) * 5
    h = np.array([np.r_[-(s - rx) / np.linalg.norm(s - rx), 1.0] for s in sats])
    try:
        q = np.linalg.inv(h.T @ h)
    except np.linalg.LinAlgError:
        return (np.inf,) * 5
    lat, lon, _ = ecef_to_geodetic(rx)
    la, lo = np.radians(lat), np.radians(lon)
    enu = np.array([[-np.sin(lo), np.cos(lo), 0],
                    [-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)],
                    [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]])
    qe = enu @ q[:3, :3] @ enu.T
    return (np.sqrt(np.trace(q)), np.sqrt(np.trace(q[:3, :3])), np.sqrt(qe[0, 0] + qe[1, 1]),
            np.sqrt(qe[2, 2]), np.sqrt(q[3, 3]))


def solve(t_tx, ephs, t_rx, iterations=10):
    """Least-squares position and clock from transmit times at one receive instant.

    t_tx: per satellite, the GPS time of week at which the signal received at t_rx left it (from the code
    phase). t_rx: the receiver's time of week for that instant, any guess within tens of ms. Returns the
    ECEF position, the receiver clock bias (s, true = t_rx - bias), the residuals (m) and the elevations
    (None when the solve never came near the Earth); a NaN position and bias when a range is not
    finite.
    """
    x = np.zeros(4)
    for it in range(iterations):
        rows, res, els = [], [], []
        for t, eph in zip(t_tx, ephs):
            pos, clk = sat_position(eph, t)
            pos2, clk = sat_position(eph, t - clk)
            rho = C * (t_rx - (t - clk))
            rng = np.linalg.norm(pos2 - x[:3])
            rot = OMEGA_E * rng / C
            sat = np.array([pos2[0] * np.cos(rot) + pos2[1] * np.sin(rot),
                            -pos2[0] * np.sin(rot) + pos2[1] * np.cos(rot), pos2[2]])
            rng = np.linalg.norm(sat - x[:3])
            tropo = 0.0
            if it > 2 and np.linalg.norm(x[:3]) > 6e6:
                el = elevation(x[:3], sat)
                tropo = 2.47 / (np.sin(el) + 0.0121)
                els.append(el)
            res.append(rho - (rng + x[3] + tropo))
            rows.append(np.r_[-(sat - x[:3]) / rng, 1.0])
        h, r = np.array(rows), np.array(res)
        if not (np.all(np.isfinite(h)) and np.all(np.isfinite(r))):
            return np.full(3, np.nan), np.nan, r, None
        dx = np.linalg.lstsq(h, r, rcond=None)[0]
        x += dx
        if it > 4 and np.linalg.norm(dx) < 1e-4:
            break
    return x[:3], x[3] / C, r, np.degrees(els) if els else None


def solve_clock(t_tx, ephs, t_rx, rx):
    """The receiver clock bias (s) with the position held at rx, as a timing receiver does once it has
    surveyed its position; the mean over satellites of each one's range residual. Returns the bias and
    the per-satellite residuals about it (m)."""
    res = []
    for t, eph in zip(t_tx, ephs):
        _, clk = sat_position(eph, t)
        pos, clk = sat_position(eph, t - clk)
        rho = C * (t_rx - (t - clk))
        rot = OMEGA_E * np.linalg.norm(pos - rx) / C
        sat = np.array([pos[0] * np.cos(rot) + pos[1] * np.sin(rot),
                        -pos[0] * np.sin(rot) + pos[1] * np.cos(rot), pos[2]])
        el = elevation(rx, sat)
        res.append(rho - np.linalg.norm(sat - rx) - 2.47 / (np.sin(el) + 0.0121))
    res = np.array(res)
    return res.mean() / C, res - res.mean()


def consistent(subs):
    """The subframes whose TOW count agrees with the majority's count per 300-bit step.

    A preamble and two parity checks can pass by chance in noise (about once in several thousand
    positions); such a header's TOW does not fit the sequence and is dropped."""
    if not subs:
        return subs
    keys = [s["tow"] - s["bit"] // 300 for s in subs]
    vals, counts = np.unique(keys, return_counts=True)
    best = vals[np.argmax(counts)]
    return [s for s, k in zip(subs, keys) if k == best and s["bit"] % 300 == subs[keys.index(best)]["bit"] % 300]
