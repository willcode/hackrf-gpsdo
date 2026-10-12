#
# SPDX-License-Identifier: BSD-3-Clause

"""Register 6 as a byte-stream command port for the GPS image, in the engine's clock domain.

The firmware lets the host write registers 1-6 only, and register 6 is the image's own, so commands
arrive one register write at a time. Each written byte is:
  bit 7  toggle: differs from the previous byte's, so a repeated value is still a new byte
  bit 6  start: this byte begins a command
  5-0    six data bits
The register is 0 after a load, which is no command: the image then behaves as the receive-only
image. The firmware's SPI self-test writes 0xA5 and restores the old value; both are payload bytes
outside a command and are ignored.

A start byte's data is the opcode (bits 5-4) and the channel (bits 3-0). Payload bytes follow, six bits
each, least significant first:
  0 arm      120 bits (20 bytes): W 0-31, Wc 32-63, G1 64-73, G2 74-83, window 84-86, start 87-118
  1 pending   66 bits (11 bytes): W 0-31, Wc 32-63
  2 stop       0 bytes
  3 setting   12 bits (2 bytes): which 0-3, value 4-11. Setting 0 selects the output mode.
A start byte in the middle of a command abandons it. Arm, pending and stop go to the engine's config
port and are held until it takes them; a command whose start byte arrives while one is still held
is dropped and counted in `dropped`, and so is an arm, pending or stop for a channel the engine does
not have (7 and above), which it would never take.

The register value crosses from `sync` by a two-stage synchronizer on bit 7. The other bits are read
two cycles after the toggle is seen; they were written together with it, and the host's writes come
milliseconds apart.
"""

from amaranth               import Module, Signal, Cat
from amaranth.lib           import wiring
from amaranth.lib.cdc       import FFSynchronizer
from amaranth.lib.wiring    import In, Out

from dsp.gps_engine         import CHANNELS

LENGTHS = {0: 20, 1: 11, 2: 0, 3: 2}
PAYLOAD_BITS = 120


class GPSCommandPort(wiring.Component):
    reg:       In(8)       # register 6, from the sync domain

    cfg_valid: Out(1)
    cfg_ready: In(1)
    cfg_kind:  Out(2)
    cfg_ch:    Out(4)
    cfg_W:     Out(32)
    cfg_Wc:    Out(32)
    cfg_g1:    Out(10)
    cfg_g2:    Out(10)
    cfg_win:   Out(3)
    cfg_start: Out(32)

    set_valid: Out(1)      # pulses when a setting command completes
    set_which: Out(4)
    set_value: Out(8)
    set_ch:    Out(4)

    dropped:   Out(8)

    def elaborate(self, platform):
        m = Module()

        tog = Signal()
        m.submodules.sync_tog = FFSynchronizer(self.reg[7], tog)
        tog_seen = Signal()
        new = Signal(2)                   # delay so the data bits are sampled after the toggle
        m.d.sync += new.eq(Cat(tog != tog_seen, new[0]))
        with m.If(tog != tog_seen):
            m.d.sync += tog_seen.eq(tog)
        byte = Signal(7)
        strobe = Signal()
        m.d.comb += strobe.eq(new[1])
        with m.If(new[0]):
            m.d.sync += byte.eq(self.reg[:7])

        op = Signal(2)
        ch = Signal(4)
        remaining = Signal(5)
        pos = Signal(5)
        busy = Signal()
        payload = Signal(PAYLOAD_BITS)
        held = Signal()

        fields = dict(W=(0, 32), Wc=(32, 64), g1=(64, 74), g2=(74, 84), win=(84, 87), start=(87, 119))
        m.d.comb += [
            self.cfg_valid.eq(held),
            self.cfg_kind.eq(op),
            self.cfg_ch.eq(ch),
        ]
        for name, (lo, hi) in fields.items():
            m.d.comb += getattr(self, f"cfg_{name}").eq(payload[lo:hi])
        with m.If(held & self.cfg_ready):
            m.d.sync += held.eq(0)

        finish = Signal()
        with m.If(strobe):
            with m.If(byte[6]):
                # Start byte. While a request is held its fields must not move, so a new command
                # waits for the hold to clear before it may overwrite the payload; one that cannot
                # be accepted now is dropped.
                with m.If(held):
                    m.d.sync += [busy.eq(0), self.dropped.eq(self.dropped + 1)]
                with m.Else():
                    m.d.sync += [op.eq(byte[4:6]), ch.eq(byte[:4]), pos.eq(0)]
                    with m.Switch(byte[4:6]):
                        for code, n in LENGTHS.items():
                            with m.Case(code):
                                m.d.sync += remaining.eq(n)
                                if n == 0:
                                    m.d.comb += finish.eq(1)
                                    m.d.sync += busy.eq(0)
                                else:
                                    m.d.sync += busy.eq(1)
            with m.Elif(busy):
                for i in range(PAYLOAD_BITS // 6):
                    with m.If(pos == i):
                        m.d.sync += payload[6 * i:6 * i + 6].eq(byte[:6])
                m.d.sync += [pos.eq(pos + 1), remaining.eq(remaining - 1)]
                with m.If(remaining == 1):
                    m.d.sync += busy.eq(0)
                    m.d.comb += finish.eq(1)

        # A finished command: settings act at once; the others are held for the engine.
        done = Signal()
        m.d.sync += done.eq(finish)
        with m.If(done):
            with m.If(op == 3):
                m.d.comb += [self.set_valid.eq(1)]
            with m.Elif(held | (ch >= CHANNELS)):
                m.d.sync += self.dropped.eq(self.dropped + 1)
            with m.Else():
                m.d.sync += held.eq(1)
        m.d.comb += [
            self.set_which.eq(payload[0:4]),
            self.set_value.eq(payload[4:12]),
            self.set_ch.eq(ch),
        ]
        return m
