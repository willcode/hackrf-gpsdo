#!/usr/bin/env python3
#
# SPDX-License-Identifier: BSD-3-Clause

"""SBAS L1 messages from 1 ms prompt sums (RTCA DO-229): symbols, Viterbi, CRC-24Q.

An SBAS message is 250 bits at 250 bit/s: an 8-bit preamble (0x53, 0x9A, 0xC6 in turn), a 6-bit
message type, 212 data bits and a 24-bit CRC-24Q. It is sent rate 1/2 convolutionally coded,
constraint length 7, as 500 symbols a second (2 ms each, two code periods).

The decoder does not assume the conventions it cannot check: which symbol of a pair comes first, whether
the G2 symbol is inverted, and the carrier's sign (Costas). It tries each, and only messages whose
CRC passes are reported. A CRC-24Q pass in a 250-bit window happens by chance about once in 16 million,
so passing messages are real.

Usage: sbas.py dumps.npz PRN
"""

import sys

import numpy as np

G1, G2 = 0o171, 0o133
PREAMBLES = (0x53, 0x9A, 0xC6)
CRC24Q = 0x1864CFB


def _parity(v):
    return bin(v).count("1") & 1


# Encoder outputs per (state, input): the 7-bit register is the input bit above the six previous ones.
_OUT = np.zeros((64, 2, 2), dtype=np.int8)
for _s in range(64):
    for _u in range(2):
        _r = (_u << 6) | _s
        _OUT[_s, _u] = (_parity(_r & G1), _parity(_r & G2))


def viterbi(soft):
    """Bits from soft symbol pairs (positive meaning bit 0 sent as +1), G1 symbol first."""
    n = len(soft) // 2
    pairs = np.asarray(soft[:2 * n], dtype=float).reshape(n, 2)
    metric = np.zeros(64)               # the start state is unknown mid-stream
    prev = np.zeros((n, 64), dtype=np.int16)
    sign = 1 - 2 * _OUT.astype(float)   # (state, input, symbol) as +-1
    states = np.arange(64)
    for k in range(n):
        # Next state = (input << 5) | (state >> 1); every next state has two predecessors.
        branch = sign @ pairs[k]        # (state, input)
        cand = np.full((64, 2), -1e18)
        src = np.zeros((64, 2), dtype=np.int16)
        for u in range(2):
            nxt = (u << 5) | (states >> 1)
            lsb = states & 1
            cand[nxt, lsb] = metric + branch[:, u]
            src[nxt, lsb] = states
        pick = np.argmax(cand, axis=1)
        metric = cand[states, pick]
        prev[k] = src[states, pick]
    bits = np.zeros(n, dtype=np.uint8)
    s = int(np.argmax(metric))
    for k in range(n - 1, -1, -1):
        bits[k] = s >> 5
        s = int(prev[k, s])
    return bits


def crc24q(bits):
    reg = 0
    for b in bits:
        reg ^= int(b) << 23
        reg = ((reg << 1) ^ CRC24Q) & 0xFFFFFF if reg & 0x800000 else (reg << 1) & 0xFFFFFF
    return reg


def messages(bits):
    """Every 250-bit window starting with a preamble whose CRC-24Q is zero over all 250 bits."""
    out = []
    for i in range(len(bits) - 250 + 1):
        pre = int("".join(map(str, bits[i:i + 8])), 2)
        if pre in PREAMBLES and crc24q(bits[i:i + 250]) == 0:
            mtype = int("".join(map(str, bits[i + 8:i + 14])), 2)
            out.append(dict(bit=i, preamble=pre, type=mtype))
    return out


def decode(prompt, edge):
    """Try both symbol pairings, G2 inverted or not, and both carrier signs; return the variant with
    the most CRC-passing messages, as (messages, variant description)."""
    n = (len(prompt) - edge) // 2
    sym = np.real(prompt[edge:edge + 2 * n].reshape(n, 2).sum(axis=1))
    sym = sym / np.median(np.abs(sym))
    best = ([], None)
    for shift in (0, 1):
        for inv2 in (False, True):
            s = sym[shift:].copy()
            s = s[:len(s) // 2 * 2]
            if inv2:
                s[1::2] = -s[1::2]
            bits = viterbi(s)
            for flip in (0, 1):
                found = messages(bits ^ flip)
                if len(found) > len(best[0]):
                    best = (found, f"pair offset {shift}, G2 {'inverted' if inv2 else 'as is'}, "
                                   f"bits {'inverted' if flip else 'as is'}")
    return best


def main():
    z = np.load(sys.argv[1])
    prn = int(sys.argv[2])
    p = z[f"{prn}_IP"] + 1j * z[f"{prn}_QP"]
    from_ms = 2000
    found, how, edge = [], None, from_ms
    for e in (from_ms, from_ms + 1):         # the symbol edge, if the dumps do not carry it
        f, h = decode(p, e)
        if len(f) > len(found):
            found, how, edge = f, h, e
    print(f"PRN {prn}: {len(found)} messages pass CRC-24Q in {(len(p) - edge) / 1000:.0f} s"
          f" of symbols ({how})")
    if found:
        types = [m["type"] for m in found]
        print("  types:", " ".join(str(t) for t in types))
        steps = np.diff([m["bit"] for m in found])
        print("  spacing in bits:", sorted(set(int(v) for v in steps)))


if __name__ == "__main__":
    main()
