#
# SPDX-License-Identifier: BSD-3-Clause

"""The register-6 command port: bytes written as the host will write them come out as engine requests
and settings, with the toggle protocol and resynchronization."""

import unittest

from amaranth.sim import Simulator

from dsp.gps_cmd import GPSCommandPort, LENGTHS


def encode(op, ch, value_bits=0):
    """The register-6 bytes of one command, without the toggle bit."""
    out = [0x40 | (op << 4) | ch]
    for i in range(LENGTHS[op]):
        out.append((value_bits >> (6 * i)) & 0x3F)
    return out


def arm_bits(W, Wc, g1, g2, win, start):
    return (W | (Wc << 32) | (g1 << 64) | (g2 << 74) | (win << 84) | (start << 87))


class TestCommandPort(unittest.TestCase):
    def run_bytes(self, seq, hold_cycles=50):
        dut = GPSCommandPort()
        sim = Simulator(dut)
        sim.add_clock(1 / 32e6)
        got, sets = [], []
        self.dropped = None

        async def bench(ctx):
            tog = 0
            for b in seq:
                tog ^= 1
                ctx.set(dut.reg, b | (tog << 7))
                for _ in range(8):
                    if ctx.get(dut.set_valid):
                        sets.append((ctx.get(dut.set_ch), ctx.get(dut.set_which), ctx.get(dut.set_value)))
                    await ctx.tick()
            for _ in range(hold_cycles):
                if ctx.get(dut.cfg_valid):
                    got.append({k: ctx.get(getattr(dut, f"cfg_{k}"))
                                for k in ("kind", "ch", "W", "Wc", "g1", "g2", "win", "start")})
                    ctx.set(dut.cfg_ready, 1)
                    await ctx.tick()
                    ctx.set(dut.cfg_ready, 0)
                await ctx.tick()
            self.dropped = ctx.get(dut.dropped)

        sim.add_testbench(bench)
        sim.run()
        return got, sets

    def test_arm(self):
        f = dict(W=0x8000303A, Wc=0xFFF5C28F, g1=0x3F0, g2=0x155, win=5, start=0xDEADBEEF)
        got, _ = self.run_bytes(encode(0, 6, value_bits=arm_bits(**f)))
        self.assertEqual(len(got), 1)
        self.assertEqual(got[0], dict(kind=0, ch=6, **f))

    def test_pending_after_abandoned_command_and_selftest_bytes(self):
        seq = [0x25, 0x25]                                   # self-test style payload bytes: ignored
        seq += encode(0, 2, value_bits=123)[:5]              # an arm cut short
        seq += encode(1, 3, value_bits=0x1234_5678 | (0x0BAD_F00D << 32))
        got, _ = self.run_bytes(seq)
        self.assertEqual(len(got), 1)
        self.assertEqual((got[0]["kind"], got[0]["ch"], got[0]["W"], got[0]["Wc"]),
                         (1, 3, 0x1234_5678, 0x0BAD_F00D))

    def test_setting_and_stop(self):
        seq = encode(3, 1, value_bits=0 | (0x01 << 4)) + encode(2, 4)
        got, sets = self.run_bytes(seq)
        self.assertEqual(sets, [(1, 0, 1)])
        self.assertEqual([(r["kind"], r["ch"]) for r in got], [(2, 4)])

    def test_command_while_held_is_dropped(self):
        # The engine takes nothing while the bytes arrive: the first request is held, and the
        # second, begun meanwhile, is dropped and counted.
        seq = encode(1, 3, value_bits=7) + encode(2, 5)
        got, _ = self.run_bytes(seq)
        self.assertEqual([(r["kind"], r["ch"], r["W"]) for r in got], [(1, 3, 7)])
        self.assertEqual(self.dropped, 1)

    def test_channel_beyond_the_engine_is_dropped(self):
        # The engine has channels 0-6: requests for 7 and above are dropped, not held, and a
        # setting is unaffected by its channel field.
        seq = encode(0, 7, value_bits=1) + encode(2, 9) + encode(1, 15, value_bits=2)
        seq += encode(3, 7, value_bits=1 | (0x03 << 4)) + encode(2, 6)
        got, sets = self.run_bytes(seq)
        self.assertEqual([(r["kind"], r["ch"]) for r in got], [(2, 6)])
        self.assertEqual(sets, [(7, 1, 3)])
        self.assertEqual(self.dropped, 3)


if __name__ == "__main__":
    unittest.main()
