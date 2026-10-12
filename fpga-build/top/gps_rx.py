#
# SPDX-License-Identifier: BSD-3-Clause

"""Receive-only image with the GPS/SBAS correlator (dsp/gps_engine.py), kind 0x43.

The decimation chain keeps half-band stages 1 to 3 (decimation up to 8, the firmware's choice at
4 Msps). The engine and its command port run in `adclk`, which at 4 Msps is 32 MHz: exactly 8 cycles
a sample, locked to the TCXO, unlike `sync` (the iCE40's own oscillator, +-20% on this part).

Register 0x06 is the command port (dsp/gps_cmd.py). It is 0 after every load, which is no command,
and the output is then the plain I/Q stream of the receive-only image. Setting 0 selects the output:
0 for I/Q, 1 for the engine's dump records (16-bit words, dsp/gps_engine.py). Register 0x07 reads
the command port's count of dropped commands (8 bits, wrapping); it is not writable.

The platform resources, their pins and their directions are those of the standard image, so no pin
drives a net that the standard image leaves as an input.
"""

from amaranth               import Elaboratable, Module, Signal, DomainRenamer
from amaranth.lib.cdc       import FFSynchronizer
from amaranth.lib.wiring    import connect

from amaranth_future        import fixed

from board                  import PralinePlatform
from interface              import MAX586xInterface, SGPIOInterface, SPIRegisterInterface
from dsp.fir                import HalfBandDecimator
from dsp.dc_block           import DCBlock
from dsp.quarter_shift      import QuarterShift
from util                   import ClockConverter, IQSample
from util.clocks            import ClockDomainGeneratorGB
from util.identity          import add_identity, KIND_GPS, CAP_RX, CAP_STANDARD_MAP
from dsp.gps_engine         import GPSEngine
from dsp.gps_cmd            import GPSCommandPort


class Top(Elaboratable):

    def elaborate(self, platform):
        m = Module()

        m.submodules.clkgen = ClockDomainGeneratorGB()
        adc_clk = "adclk"
        dac_clk = "daclk"

        # The ADC data pins include 39 to 41, which are also the outputs of the
        # RGB current-sink drivers. They are plain inputs here; no SB_RGBA_DRV,
        # SB_RGB_DRV or SB_LEDDA_IP is instantiated.
        m.submodules.adcdac_intf = adcdac_intf = MAX586xInterface(adc_domain=adc_clk, dac_domain=dac_clk)
        # host_data is shared with the MCU's SGPIO pins. The interface enables
        # its drivers only while the MCU's direction input is low, and releases
        # them when it is high. Its DAC stream stays unconnected, so with the
        # direction high it accepts no samples and drives no clock pulses.
        m.submodules.mcu_intf    = mcu_intf    = SGPIOInterface(sample_width=16, domain="sync")

        m.d.comb += adcdac_intf.q_invert.eq(platform.request("q_invert").i)

        taps = [-2, 0, 7, 0, -18, 0, 41, 0, -92, 0, 320, 512, 320, 0, -92, 0, 41, 0, -18, 0, 7, 0, -2]
        taps = [tap / 1024 for tap in taps]
        taps2 = [3, 0, -16, 0, 77, 128, 77, 0, -16, 0, 3]
        taps2 = [tap / 256 for tap in taps2]
        taps3 = [-9, 0, 73, 128, 73, 0, -9]
        taps3 = [tap / 256 for tap in taps3]

        common_rx_filter_opts = dict(
            data_shape=fixed.SQ(1, 7),
            shape_out=fixed.SQ(1, 7),
            always_ready=True,
            domain=adc_clk,
        )

        rx_chain = {
            "dc_block":      DCBlock(width=8, num_channels=2, domain=adc_clk),
            "quarter_shift": DomainRenamer(adc_clk)(QuarterShift()),
            # Half-band stages 3 to 1: decimation up to 8. Stages 4 and 5 (below 2.5 Msps) are left
            # out to make room for the correlator; a larger setting acts as 8.
            "hbfir3":        HalfBandDecimator(taps3, **common_rx_filter_opts),
            "hbfir2":        HalfBandDecimator(taps2, **common_rx_filter_opts),
            "hbfir1":        HalfBandDecimator(taps, **common_rx_filter_opts),
            "clkconv":       ClockConverter(IQSample(8), 8, adc_clk, "sync"),
        }
        for k, v in rx_chain.items():
            m.submodules[f"rx_{k}"] = v

        last = adcdac_intf.adc_stream
        decimated = None
        for name, block in rx_chain.items():
            if name == "clkconv":
                decimated = last
            connect(m, last, block.input)
            last = block.output

        # ---- correlator, command port and dump output, in adclk ----------------------------
        m.submodules.engine = engine = DomainRenamer(adc_clk)(GPSEngine())
        m.submodules.cmd = cmd = DomainRenamer(adc_clk)(GPSCommandPort())
        m.d.comb += [
            engine.in_valid.eq(decimated.valid),
            # The decimator's payload is an array of two fixed-point values, I then Q, with the
            # same bits as IQSample(8).
            engine.in_i.eq(decimated.p.as_value()[:8].as_signed()),
            engine.in_q.eq(decimated.p.as_value()[8:16].as_signed()),
        ]
        for k in ("valid", "kind", "ch", "W", "Wc", "g1", "g2", "win", "start"):
            m.d.comb += getattr(engine, f"cfg_{k}").eq(getattr(cmd, f"cfg_{k}"))
        m.d.comb += cmd.cfg_ready.eq(engine.cfg_ready)

        # Setting 1 asks the engine for a snapshot of (value + 1) * 1024 input samples.
        m.d.comb += [
            engine.snap_valid.eq(cmd.set_valid & (cmd.set_which == 1)),
            engine.snap_value.eq(cmd.set_value),
        ]

        mode_ad = Signal(2)
        with m.If(cmd.set_valid & (cmd.set_which == 0)):
            m.d[adc_clk] += mode_ad.eq(cmd.set_value[:2])
        mode = Signal(2)
        m.submodules.mode_sync = FFSynchronizer(mode_ad, mode, o_domain="sync")

        # Record words across to sync by a toggle handshake, one word at a time: the engine sends
        # at most one word a sample (eight adclk cycles), and a word's round trip takes about five
        # cycles of either clock. The word is held from its toggle until the acknowledgment, so it
        # needs no synchronizer. (Amaranth's asynchronous FIFO missed 48 MHz in its own pointer
        # logic at this fill.)
        hs_word = Signal(16)
        hs_req = Signal()
        hs_ack = Signal()
        hs_ack_ad = Signal()
        hs_req_s = Signal()
        m.submodules.hs_req_sync = FFSynchronizer(hs_req, hs_req_s, o_domain="sync")
        m.submodules.hs_ack_sync = FFSynchronizer(hs_ack, hs_ack_ad, o_domain=adc_clk)
        m.d.comb += engine.out_ready.eq(hs_req == hs_ack_ad)
        with m.If(engine.out_valid & (hs_req == hs_ack_ad)):
            m.d[adc_clk] += [hs_word.eq(engine.out_data), hs_req.eq(~hs_req)]
        word_new = Signal()
        m.d.comb += word_new.eq(hs_req_s != hs_ack)
        m.d.sync += hs_ack.eq(hs_req_s)

        # The selected stream is registered before the interface: combinationally, the crossing,
        # the mode select and the interface's state machine made one long path. Records are
        # discarded while I/Q is selected.
        sel_valid = Signal()
        sel_p = Signal(16)
        m.d.comb += [
            mcu_intf.adc_stream.valid.eq(sel_valid),
            mcu_intf.adc_stream.p.eq(sel_p),
        ]
        with m.If(mode == 1):
            m.d.sync += [sel_valid.eq(word_new), sel_p.eq(hs_word)]
        with m.Else():
            m.d.sync += [sel_valid.eq(last.valid), sel_p.eq(last.p)]

        # The DAC stays at mid-scale so the transmit path carries no signal.
        m.d.comb += adcdac_intf.dac_stream.valid.eq(0)

        # SPI register interface: one SB_SPI core on the board's SPI pins 14 to 17.
        spi_port = platform.request("spi")
        m.submodules.spi_regs = spi_regs = SPIRegisterInterface(spi_port, registered_read=2)
        ctrl     = spi_regs.add_register(0x01, init=0)
        rx_decim = spi_regs.add_register(0x02, init=0, size=3)
        rx_pstep = spi_regs.add_register(0x03, init=0)
        tx_ctrl  = spi_regs.add_register(0x04, init=0, size=1)
        tx_intrp = spi_regs.add_register(0x05, init=0, size=3)
        tx_pstep = spi_regs.add_register(0x06, init=0)
        add_identity(spi_regs, KIND_GPS, CAP_RX | CAP_STANDARD_MAP)

        m.d.comb += cmd.reg.eq(tx_pstep)

        # The dropped count crosses to sync in Gray code, so a read taken during a step is the old
        # or the new value, never a mix.
        dropped_gray_ad = Signal(8)
        m.d[adc_clk] += dropped_gray_ad.eq(cmd.dropped ^ (cmd.dropped >> 1))
        dropped_gray = Signal(8)
        m.submodules.dropped_sync = FFSynchronizer(dropped_gray_ad, dropped_gray, o_domain="sync")
        dropped = Signal(8)
        for i in range(8):
            m.d.sync += dropped[i].eq(dropped_gray[i:].xor())
        spi_regs.add_read_only_register(0x07, read=dropped)

        m.d.sync += [
            mcu_intf.trigger_en                 .eq(ctrl[7]),
            mcu_intf.prbs                       .eq(ctrl[6]),
            rx_chain["dc_block"].enable         .eq(ctrl[0]),
            rx_chain["quarter_shift"].enable    .eq(rx_pstep[-2]),
            rx_chain["quarter_shift"].up        .eq(rx_pstep[-1]),
        ]
        rx_enables = Signal(5)
        m.d.sync += rx_enables.eq((1 << rx_decim) - 1)
        m.d.comb += [
            rx_chain["hbfir3"].enable           .eq(rx_enables[2]),
            rx_chain["hbfir2"].enable           .eq(rx_enables[1]),
            rx_chain["hbfir1"].enable           .eq(rx_enables[0]),
        ]

        return m


if __name__ == "__main__":
    PralinePlatform().build(Top(), name="14_gps", build_dir="build-gps")
