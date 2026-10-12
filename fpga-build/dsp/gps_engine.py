#
# SPDX-License-Identifier: BSD-3-Clause

"""Time-shared GPS/SBAS correlator: the arithmetic of gps_track_model.Channel for seven channels.

Each sample is worked in eight slots, one a cycle: slots 0 to 6 are the channels, slot 7 reads out
the dump records (and slot 6's turn on the words' memory takes the config writes). A channel's state
lives in block RAM and is never copied into registers: its slot reads the state, and the next two
cycles compute and write it back. Its sums are read and written four and five cycles later, after the
multipliers. A channel is visited every eight cycles, so no read meets an unfinished write.

Per sample, exactly as gps_track_model.Channel:
- prompt chip from the code window c(k-1), c(k), c(k+1); early is c(k+1) when the code fraction's top
  bit is set, else c(k); late is c(k) when it is set, else c(k-1);
- carrier table entry theta >> 26; yi = xi C + xq S, yq = xq C - xi S;
- E, P, L sums of +-yi, +-yq (chip 0 is +1);
- then the fraction advances by W; a carry steps the chip, the window and G1 and G2. The carrier
  phase advances by Wc.
A sample whose chip count is 1023 starts a new period. Pending words, if any, take effect for the
advance after that sample, gps_track_model's "period after next".

Memories (all free of read-write collisions by the slot schedule, so the block RAMs are instantiated
directly; util/ebr.py):
  nco     NCO_FIELDS per channel: read in the slot, written the next cycle (stage F)
  code    CODE_FIELDS per channel: read in the slot, written two cycles later (stage G)
  words   W (31 bits; bit 31 is 0) and Wc (bits 32-63) per channel in two halves, in SPRAM. The
          channel's current half is read in the slot before its own; slot 6 writes the other half
          for a pending request, or the current half for an arm. The visit before a boundary (which
          `at1023` announces) switches halves when a request is pending, so the boundary's own visit
          reads the new words.
  sums    six SUM_BITS sums per channel and bank; a period accumulates in one bank, so the finished
          period's sums stay put in the other until they are read out
  capture the record's timing part, written at the boundary: its sample index, code fraction, W, Wc

Config port: one request at a time, held by the requester (cfg_valid and the fields) until cfg_ready
pulses.
  kind 0 (arm): W (below 2^31), Wc, G1, G2, window, and `cfg_start`, the index of the sample BEFORE
    the channel's first. Slot 6 writes the words, then the channel's own slot writes the phase state:
    inactive, waiting for that index. At that sample's visit the channel is set to chip 0, fraction 0,
    carrier phase 0; the next sample is the first of its period 0, an ordinary advance.
  kind 1 (pending): W, Wc, for the next period boundary; taken in slot 6, once the channel's previous
    pending request has been applied. One to a channel that is not tracking (not yet started, or
    stopped) is taken and discarded, so it cannot wait for a boundary that never comes.
  kind 2 (stop): taken in the channel's slot.
The channel is cfg_ch's low three bits; the command port drops requests to channel 7 and above.

Snapshots: snap_valid with snap_value v sends the next (v + 1) * 1024 input samples, each as a word
(I low byte, Q high byte) when it enters the slots, then a three-word trailer from which
decode_snapshot recovers the first sample's index. A request during a snapshot waits for its
trailer, and its value is kept apart until its own snapshot begins. Record readout waits meanwhile, so a channel
that reaches two boundaries in that time loses a record (dump_lost). Snapshots serve acquisition: a
search in them gives the code phase at a known sample index, as an arm requires. A record does not
start in its channel's last 15 chips, so one that waited out a snapshot is either sent before the
next boundary or lost to it, never read across it.

Dump records: a boundary marks its channel and captures the new period's start; the readout takes
one marked channel at a time and sends RECORD_WORDS 16-bit words, least significant first, one per
sample (slot 7 reads the capture and the finished bank once per word):
  0       0xD5 | channel << 8
  1, 2    n0: the index of the new period's first sample
  3, 4    F: its code fraction (chip 0)
  5, 6    W and 7, 8 Wc: the words in use from that sample
  9-20    IE, QE, IP, QP, IL, QL of the finished period, each sign-extended to 32 bits
"""

from amaranth               import Module, Signal, Cat, Const, Mux, signed
from amaranth.lib           import wiring
from amaranth.lib.wiring    import In, Out

from dsp.mul16ce            import Mul16CE, LATENCY as MUL_LATENCY
from dsp                    import gps_track_model as tm
from util.ebr               import EBR
from util.spram             import SPRAM

MAGIC = 0xD5
SNAP_MAGIC = 0xD7
CHANNELS = 7
SLOTS = 8
SUM_BITS = 28              # a period's sum of 17-bit terms over at most about 4,100 samples
# The phase state in two memories: the NCOs, written in stage F, and the fields the code carry
# steps, written a cycle later (stage G) from registered values.
NCO_FIELDS = [("frac", 32), ("theta", 32)]
CODE_FIELDS = [("active", 1), ("armed", 1), ("at1023", 1), ("chip", 10), ("g1", 10), ("g2", 10), ("win", 3)]
PHASE_FIELDS = NCO_FIELDS + CODE_FIELDS
SUM_NAMES = ("ie", "qe", "ip", "qp", "il", "ql")
RECORD_WORDS = 21

G1_TAPS = (3, 10)
G2_TAPS = (2, 3, 6, 8, 9, 10)


def _offsets(fields):
    off, out = 0, {}
    for name, width in fields:
        out[name] = (off, width)
        off += width
    return out


NCO_OFF = _offsets(NCO_FIELDS)
CODE_OFF = _offsets(CODE_FIELDS)


def lfsr_step(reg, taps):
    fb = reg[taps[0] - 1]
    for t in taps[1:]:
        fb = fb ^ reg[t - 1]
    return Cat(fb, reg[:9])


def sign_extend(v, width):
    return Cat(v, v[-1].replicate(width - len(v)))


def decode_record(words):
    """A dump record's 16-bit words as a dict with the model's dump fields for the new period's start
    (n_start, F_start), its words (W, Wc) and the finished period's sums."""
    w = [int(v) for v in words]
    assert w[0] & 0xFF == MAGIC, "bad magic"
    u32 = lambda lo: w[lo] | (w[lo + 1] << 16)
    s32 = lambda lo: u32(lo) - (1 << 32) if u32(lo) & (1 << 31) else u32(lo)
    rec = dict(ch=w[0] >> 8, n_start=u32(1), F_start=u32(3), W=u32(5), Wc=s32(7))
    for i, name in enumerate(("IE", "QE", "IP", "QP", "IL", "QL")):
        rec[name] = s32(9 + 2 * i)
    return rec


def decode_snapshot(words):
    """A snapshot's words (its samples, then the three trailer words) as (index of the first
    sample, complex samples)."""
    w = [int(v) for v in words]
    assert w[-3] & 0xFF == SNAP_MAGIC, "bad snapshot magic"
    count = ((w[-3] >> 8) + 1) * 1024
    assert len(w) == count + 3, "snapshot length does not match its trailer"
    lo, hi = w[-2], w[-1]
    second = ((((hi - (lo == 0xFFFF)) & 0xFFFF) << 16) | lo)  # index at the trailer's second word
    last = (second - 2) & 0xFFFFFFFF
    s8 = lambda v: v - 256 if v & 0x80 else v
    x = [complex(s8(v & 0xFF), s8(v >> 8)) for v in w[:count]]
    return (last - count + 1) & 0xFFFFFFFF, x


class GPSEngine(wiring.Component):
    in_valid:  In(1)
    in_i:      In(signed(8))
    in_q:      In(signed(8))

    cfg_valid: In(1)
    cfg_ready: Out(1)
    cfg_kind:  In(2)
    cfg_ch:    In(4)
    cfg_W:     In(32)
    cfg_Wc:    In(32)
    cfg_g1:    In(10)
    cfg_g2:    In(10)
    cfg_win:   In(3)
    cfg_start: In(32)

    out_valid: Out(1)
    out_ready: In(1)
    out_data:  Out(16)

    snap_valid: In(1)       # pulse: send a snapshot of (snap_value + 1) * 1024 input samples
    snap_value: In(8)

    overrun:   Out(1)       # a sample arrived before the last one's slots began
    snap_lost: Out(1)       # a snapshot word found the output still full
    dump_lost: Out(1)       # a boundary came before the previous one's record was read out

    def elaborate(self, platform):
        m = Module()

        # ---- sample intake and slot sequencer ---------------------------------------------
        # Samples come every eight cycles at most, so one register holds the next.
        x_full = Signal()
        x_i, x_q = Signal(signed(8)), Signal(signed(8))
        xi, xq = Signal(signed(8)), Signal(signed(8))
        n = Signal(32, init=(1 << 32) - 1)       # index of the sample in the slots
        slot = Signal(3)
        busy = Signal()                          # slots of a sample are running
        r_valid = Signal()                       # stage R holds a slot this cycle
        m.d.comb += r_valid.eq(busy)
        with m.If(busy):
            m.d.sync += slot.eq(slot + 1)
            with m.If(slot == SLOTS - 1):
                m.d.sync += busy.eq(0)
        start_next = Signal()
        m.d.comb += start_next.eq(x_full & (~busy | (slot == SLOTS - 1)))
        with m.If(start_next):
            m.d.sync += [busy.eq(1), slot.eq(0), n.eq(n + 1), xi.eq(x_i), xq.eq(x_q), x_full.eq(0)]
        with m.If(self.in_valid):
            m.d.sync += [x_i.eq(self.in_i), x_q.eq(self.in_q), x_full.eq(1)]
            with m.If(x_full & ~start_next):
                m.d.sync += self.overrun.eq(1)

        # ---- memories -----------------------------------------------------------------------
        m.submodules.nco = nco = EBR(sum(w for _, w in NCO_FIELDS), SLOTS)
        m.submodules.code = code = EBR(sum(w for _, w in CODE_FIELDS), SLOTS)
        m.submodules.sums = sums = EBR(6 * SUM_BITS, 2 * SLOTS)
        m.submodules.capture = capture = EBR(128, SLOTS)
        words = [SPRAM(sim_depth=2 * SLOTS) for _ in range(4)]
        for i, p in enumerate(words):
            m.submodules[f"words{i}"] = p
        wd_addr = Signal(4)
        wd_we = Signal()
        wd_din = Signal(64)
        wd_dout = Cat(*[p.dout for p in words])
        wd_reg = Signal(64)                      # read a slot ahead, so registered on the way out
        m.d.sync += wd_reg.eq(wd_dout)
        for i, p in enumerate(words):
            m.d.comb += [p.addr.eq(wd_addr), p.we.eq(wd_we), p.din.eq(wd_din[16 * i:16 * i + 16])]

        half = Signal(CHANNELS)                  # each channel's current half of the words
        pend = Signal(CHANNELS)                  # a pending request waits in the other half
        par = Signal(CHANNELS)                   # each channel's accumulating bank

        # ---- stage R: phase read --------------------------------------------------------------
        for mem in (nco, code):
            m.d.comb += [mem.raddr.eq(slot), mem.re.eq(r_valid & (slot != SLOTS - 1))]

        # ---- stage F: compute and write back the phase state ---------------------------------
        # The sample and its index travel with the slot: the next sample is taken while this one's
        # last slots are still in the pipeline.
        F_valid, F_slot = Signal(), Signal(3)
        F_n, F_n1, F_xi, F_xq = Signal(32), Signal(32), Signal(signed(8)), Signal(signed(8))
        m.d.sync += [F_valid.eq(r_valid), F_slot.eq(slot), F_n.eq(n), F_n1.eq(n + 1),
                     F_xi.eq(xi), F_xq.eq(xq)]
        st = {name: nco.rdata[o:o + w] for name, (o, w) in NCO_OFF.items()}
        st.update({name: code.rdata[o:o + w] for name, (o, w) in CODE_OFF.items()})
        W_e, Wc_e = wd_reg[:32], wd_reg[32:]
        is_ch = F_valid & (F_slot != SLOTS - 1)
        act = Signal()
        pre = Signal()                           # the visit before the first sample
        # The start index is compared a visit ahead, so the 32-bit compare ends in a register: in
        # two 16-bit halves, then their AND a cycle later, well before the channel's next visit.
        go = Signal(CHANNELS)
        # The halves' slot is kept one-hot, so the write needs no decode (a binary copy merges with
        # stage G's slot, and its decode then sits on the pending-request paths).
        # Gated by F_valid: between samples the slot rests at 0, and an ungated latch would overwrite
        # channel 0's result with a stale compare, so a start just after an input gap was missed.
        go_lo, go_hi, go_at = Signal(), Signal(), Signal(CHANNELS)
        m.d.sync += [go_lo.eq(F_n1[:16] == st["frac"][:16]), go_hi.eq(F_n1[16:] == st["frac"][16:])]
        for c in range(CHANNELS):
            m.d.sync += go_at[c].eq(F_valid & (F_slot == c))
            with m.If(go_at[c]):
                m.d.sync += go[c].eq(go_lo & go_hi)
        bnd = Signal()
        m.d.comb += [
            act.eq(is_ch & st["active"]),
            pre.eq(is_ch & ~st["active"] & st["armed"] & go.bit_select(F_slot, 1)),
            bnd.eq(act & st["at1023"]),
        ]

        # An arm (after its words are written) or a stop is taken in the target channel's own slot,
        # in place of its computed state.
        arm_w_done = Signal()
        cfg_phase = Signal()
        m.d.comb += cfg_phase.eq(is_ch & self.cfg_valid & (self.cfg_kind != 1)
                                 & (self.cfg_ch[:3] == F_slot)
                                 & ((self.cfg_kind == 2) | arm_w_done))

        nxt = Signal(33)
        m.d.comb += nxt.eq(st["frac"] + W_e)
        m.d.comb += [
            nco.waddr.eq(F_slot),
            nco.wdata.eq(Cat(Mux(cfg_phase, self.cfg_start, Mux(act, nxt[:32], 0)),
                             Mux(act, st["theta"] + Wc_e, 0))),
            nco.we.eq(act | pre | cfg_phase),
        ]

        # ---- stage G: the fields the code carry steps ------------------------------------------
        # The carry lands in a register; the chip, the generators and the window follow from it a
        # cycle later. The config request is still held here, and is released at this stage.
        G = dict(slot=Signal(3), we=Signal(), cfg=Signal(), act=Signal(), pre=Signal(), step=Signal(),
                 last=Signal())
        Gst = {name: Signal(w, name=f"G_{name}") for name, w in CODE_FIELDS}
        m.d.sync += [G["slot"].eq(F_slot), G["we"].eq(act | pre | cfg_phase), G["cfg"].eq(cfg_phase),
                     G["act"].eq(act), G["pre"].eq(pre), G["step"].eq(act & nxt[32]),
                     G["last"].eq(~st["at1023"] & (st["chip"] == 1022))]
        m.d.sync += [Gst[name].eq(st[name]) for name, _ in CODE_FIELDS]
        step = G["step"]
        chip_base = Mux(Gst["at1023"], 0, Gst["chip"])
        gw = Gst["win"]
        out_bit = Gst["g1"][9] ^ Gst["g2"][9]
        new = {name: Signal(w, name=f"new_{name}") for name, w in CODE_FIELDS}
        m.d.comb += [
            new["active"].eq(~G["cfg"]),
            new["armed"].eq(Mux(G["cfg"], self.cfg_kind == 0, ~G["act"])),
            new["at1023"].eq(~G["cfg"] & (G["pre"] | (step & G["last"]))),
            new["chip"].eq(Mux(step, chip_base + 1, chip_base)),
            new["g1"].eq(Mux(G["cfg"], self.cfg_g1, Mux(step, lfsr_step(Gst["g1"], G1_TAPS), Gst["g1"]))),
            new["g2"].eq(Mux(G["cfg"], self.cfg_g2, Mux(step, lfsr_step(Gst["g2"], G2_TAPS), Gst["g2"]))),
            new["win"].eq(Mux(G["cfg"], self.cfg_win, Mux(step, Cat(gw[1], gw[2], out_bit), gw))),
        ]
        m.d.comb += [
            code.waddr.eq(G["slot"]),
            code.wdata.eq(Cat(*[new[name] for name, _ in CODE_FIELDS])),
            code.we.eq(G["we"]),
        ]

        # The visit before a boundary switches to a pending request's half, at stage G (the
        # channel's next read is seven cycles away).
        flip = Signal()
        m.d.comb += flip.eq(G["we"] & new["at1023"] & pend.bit_select(G["slot"], 1))
        with m.If(flip):
            m.d.sync += [half.bit_select(G["slot"], 1).eq(~half.bit_select(G["slot"], 1)),
                         pend.bit_select(G["slot"], 1).eq(0)]
        # A bit per tracking channel, equal to the code memory's `active` after stage G: set at a
        # channel's first sample, cleared by an arm or a stop.
        trk = Signal(CHANNELS)
        with m.If(G["pre"] | G["cfg"]):
            m.d.sync += trk.bit_select(G["slot"], 1).eq(~G["cfg"])

        # The words' port: the SPRAM sits at the die's corner, so its address comes from a register
        # and its data goes through one, and each channel's words are read in the slot before its
        # own (slot 7 or an idle cycle reads channel 0). Slot 6 has no read; it writes the words of
        # an arm (current half) or of a pending request (other half), decided the cycle before from
        # registers alone: a pending request waits while the channel's previous one is unapplied,
        # one to a channel not tracking is taken and discarded (`w_skip`), and an arm first drops
        # any pending one. A channel without a pending request cannot switch halves, so the half
        # read here stays the one written.
        cfg_ch = self.cfg_ch[:3]
        ch_pend = pend.bit_select(cfg_ch, 1)
        ch_trk = trk.bit_select(cfg_ch, 1)
        nslot = Signal(3)                        # the slot of the next cycle; 7 before a sample
        rd_next = Signal(3)
        m.d.comb += [nslot.eq(Mux(start_next, 0, Mux(busy & (slot != SLOTS - 1), slot + 1, SLOTS - 1))),
                     rd_next.eq(nslot + 1)]
        to6 = busy & (slot == SLOTS - 3) & self.cfg_valid
        w_pend, w_arm, w_skip = Signal(), Signal(), Signal()
        m.d.sync += [wd_we.eq(0), w_pend.eq(0), w_arm.eq(0), w_skip.eq(0),
                     wd_addr.eq(Cat(rd_next, half.bit_select(rd_next, 1)))]
        with m.If(to6 & ~ch_trk & (self.cfg_kind == 1)):
            m.d.sync += w_skip.eq(1)
        with m.If(to6 & ch_trk & ~ch_pend & (self.cfg_kind == 1)):
            m.d.sync += [wd_we.eq(1), w_pend.eq(1), wd_addr.eq(Cat(cfg_ch, ~half.bit_select(cfg_ch, 1)))]
        with m.If(to6 & ~arm_w_done & (self.cfg_kind == 0)):
            with m.If(ch_pend):
                m.d.sync += pend.bit_select(cfg_ch, 1).eq(0)
            with m.Else():
                m.d.sync += [wd_we.eq(1), w_arm.eq(1), wd_addr.eq(Cat(cfg_ch, half.bit_select(cfg_ch, 1)))]
        m.d.sync += wd_din.eq(Cat(self.cfg_W[:31], Const(0, 1), self.cfg_Wc))
        with m.If(w_pend):
            m.d.sync += pend.bit_select(cfg_ch, 1).eq(1)
        with m.If(w_arm):
            m.d.sync += arm_w_done.eq(1)
        with m.If(G["cfg"]):
            m.d.sync += arm_w_done.eq(0)
        m.d.comb += self.cfg_ready.eq(G["cfg"] | w_pend | w_skip)

        # The boundary: switch banks, capture the new period's start, mark the record (not for the
        # first sample after an arm, which ends no period).
        bank = Signal()
        m.d.comb += bank.eq(par.bit_select(F_slot, 1) ^ bnd)
        with m.If(bnd):
            m.d.sync += par.bit_select(F_slot, 1).eq(bank)
        m.d.comb += [
            capture.waddr.eq(F_slot),
            capture.wdata.eq(Cat(F_n, st["frac"], W_e, Wc_e)),
            capture.we.eq(bnd),
        ]
        dump_req = Signal(CHANNELS)
        sending = Signal()
        rd_ch = Signal(3)
        rd_bank = Signal()
        rd_cnt = Signal(range(RECORD_WORDS + 1))
        rd_sel = Signal(3)
        # A channel in its last 15 chips (about 58 samples: the readout's 21 and room for the
        # output's stalls) does not start a record: its boundary, landing during the readout, would clear the
        # bank being read and overwrite the start, giving the next period's start with a few
        # samples' sums. Only a record that waited out a snapshot can be that late. It waits for the
        # boundary and sends that record, the older one lost.
        near = Signal(CHANNELS)
        with m.If(G["we"]):
            m.d.sync += near.bit_select(G["slot"], 1).eq(
                G["act"] & ~Gst["at1023"] & Gst["chip"][4:].all())
        # Registered, so the request loop (request, priority, start, request) spans two cycles; a
        # start reads it a cycle late, harmless while the cycle after a start is sending.
        ready = Signal(CHANNELS)
        m.d.sync += ready.eq(dump_req & ~near)
        for c in reversed(range(CHANNELS)):
            with m.If(ready[c]):
                m.d.comb += rd_sel.eq(c)
        snap_busy = Signal()
        begin = ~sending & (ready != 0) & ~snap_busy
        with m.If(begin):
            m.d.sync += [sending.eq(1), rd_ch.eq(rd_sel), rd_cnt.eq(0),
                         rd_bank.eq(~par.bit_select(rd_sel, 1))]
        # A boundary marks its channel a cycle later, off the boundary's path; the channel is in
        # its last chips until then, so no record of it starts meanwhile.
        mark, mark_slot = Signal(), Signal(3)
        m.d.sync += [mark.eq(bnd & ~st["armed"]), mark_slot.eq(F_slot)]
        for c in range(CHANNELS):
            with m.If(begin & (rd_sel == c)):
                m.d.sync += dump_req[c].eq(0)
            with m.If(mark & (mark_slot == c)):
                with m.If(dump_req[c] | (sending & (rd_ch == c))):
                    m.d.sync += self.dump_lost.eq(1)
                m.d.sync += dump_req[c].eq(1)

        # Correlation bits, table entry, sample, to the multiplier stage.
        M = dict(valid=Signal(), slot=Signal(3), bank=Signal(), clear=Signal(), e=Signal(), p=Signal(),
                 l=Signal(), svc=Signal())
        M_xi, M_xq = Signal(signed(8)), Signal(signed(8))
        m.d.sync += [M_xi.eq(F_xi), M_xq.eq(F_xq)]
        # The carrier table as logic: as a 64-entry read-only memory, yosys built 1,024 flip-flops.
        cs = Signal(16)
        with m.Switch(st["theta"][26:32]):
            for i, (c, s_) in enumerate(zip(tm.COS, tm.SIN)):
                with m.Case(i):
                    m.d.comb += cs.eq((int(c) & 0xFF) | ((int(s_) & 0xFF) << 8))
        M_cs = Signal(16)
        m.d.sync += M_cs.eq(cs)
        frac31, sw = st["frac"][31], st["win"]
        m.d.sync += [
            M["valid"].eq(act), M["slot"].eq(F_slot), M["bank"].eq(bank), M["clear"].eq(st["at1023"]),
            M["e"].eq(Mux(frac31, sw[2], sw[1])), M["p"].eq(sw[1]),
            M["l"].eq(Mux(frac31, sw[1], sw[0])),
            M["svc"].eq(F_valid & (F_slot == SLOTS - 1)),
        ]

        # ---- multipliers ---------------------------------------------------------------------
        cos_v = M_cs[:8].as_signed()
        sin_v = M_cs[8:].as_signed()
        muls = [Mul16CE() for _ in range(4)]
        for i, mu in enumerate(muls):
            m.submodules[f"mul{i}"] = mu
            m.d.comb += mu.ce.eq(1)
        for mu, a, b in ((muls[0], M_xi, cos_v), (muls[1], M_xq, sin_v),
                         (muls[2], M_xq, cos_v), (muls[3], M_xi, sin_v)):
            m.d.comb += [mu.a.eq(a), mu.b.eq(b)]
        prod = [mu.o.as_signed() for mu in muls]
        assert MUL_LATENCY == 2

        # ---- stage Y: the products' sums yi, yq; issue the sums read (slot 7: the record's next
        # word). The two additions sit in separate cycles: in one, the multiplier outputs, the
        # 17-bit and the 28-bit carry chains made the longest path.
        M2 = {k: Signal.like(v, name=f"M2_{k}") for k, v in M.items()}
        Y = {k: Signal.like(v, name=f"Y_{k}") for k, v in M.items()}
        for k in M:
            m.d.sync += [M2[k].eq(M[k]), Y[k].eq(M2[k])]
        yi = Signal(signed(17))
        yq = Signal(signed(17))
        m.d.sync += [yi.eq(prod[0] + prod[1]), yq.eq(prod[2] - prod[3])]
        out_full = Signal()
        issue = Signal()
        m.d.comb += [
            issue.eq(Y["svc"] & sending & ~out_full),
            sums.raddr.eq(Mux(issue, Cat(rd_bank, rd_ch), Cat(Y["bank"], Y["slot"]))),
            sums.re.eq(1),
            capture.raddr.eq(rd_ch),
            capture.re.eq(issue),
        ]

        # ---- stage A: accumulate, or send a record word ----------------------------------------
        A = {k: Signal.like(v, name=f"A_{k}") for k, v in Y.items()}
        for k in Y:
            m.d.sync += A[k].eq(Y[k])
        A_issue = Signal()
        m.d.sync += A_issue.eq(issue)
        old = {name: sums.rdata[i * SUM_BITS:(i + 1) * SUM_BITS] for i, name in enumerate(SUM_NAMES)}
        acc = {}
        for (si, sq, bit) in (("ie", "qe", "e"), ("ip", "qp", "p"), ("il", "ql", "l")):
            neg = A[bit]
            for s, y in ((si, yi), (sq, yq)):
                base = Mux(A["clear"], 0, old[s])
                ys = sign_extend(y.as_unsigned(), SUM_BITS) ^ neg.replicate(SUM_BITS)
                acc[s] = Signal(SUM_BITS, name=f"acc_{s}")
                m.d.comb += acc[s].eq(base + ys + neg)
        m.d.comb += [sums.waddr.eq(Cat(A["bank"], A["slot"])),
                     sums.wdata.eq(Cat(*[acc[s] for s in SUM_NAMES])), sums.we.eq(A["valid"])]

        # ---- record output ----------------------------------------------------------------------
        record = [Cat(Const(MAGIC, 8), rd_ch, Const(0, 5))]
        record += [capture.rdata[16 * i:16 * i + 16] for i in range(8)]
        for i in range(6):
            v = sign_extend(sums.rdata[i * SUM_BITS:(i + 1) * SUM_BITS], 32)
            record += [v[:16], v[16:]]
        assert len(record) == RECORD_WORDS
        word = Signal(16)
        with m.Switch(rd_cnt):
            for i, wv in enumerate(record):
                with m.Case(i):
                    m.d.comb += word.eq(wv)
        out_data = Signal(16)
        m.d.comb += [self.out_valid.eq(out_full), self.out_data.eq(out_data)]
        with m.If(out_full & self.out_ready):
            m.d.sync += out_full.eq(0)
        # The record word is registered once more before the output register, so the snapshot
        # selection below does not follow the 21-way selection in the same cycle.
        rec_word, A2_issue = Signal(16), Signal()
        m.d.sync += [rec_word.eq(word), A2_issue.eq(A_issue)]
        with m.If(A_issue):
            m.d.sync += rd_cnt.eq(rd_cnt + 1)
        with m.If(A2_issue):
            m.d.sync += [out_data.eq(rec_word), out_full.eq(1)]
            with m.If(rd_cnt == RECORD_WORDS):
                m.d.sync += sending.eq(0)

        # ---- snapshots ----------------------------------------------------------------------
        # Raw input samples, one word per sample as each enters the slots, while record readout
        # waits; then a trailer of three words, one per following sample: SNAP_MAGIC | length
        # code << 8, then the index of the sample that entered at the second word, low half, and at
        # the third word, high half (decode_snapshot reassembles it).
        snap_req, snap_on = Signal(), Signal()
        snap_len = Signal(8)                     # the running snapshot's, which its trailer reports
        snap_next = Signal(8)                    # the waiting request's
        snap_left = Signal(18)
        trail = Signal(2)
        m.d.comb += snap_busy.eq(snap_req | snap_on | (trail != 0))
        # A cycle after a sample enters the slots, xi and xq hold it and n is its index.
        tick = Signal()
        m.d.sync += tick.eq(start_next)
        # A request during a snapshot waits for its trailer.
        first = tick & snap_req & ~sending & ~snap_on & (trail == 0)
        with m.If(first | (tick & snap_on)):
            with m.If(out_full & ~self.out_ready):
                m.d.sync += self.snap_lost.eq(1)
            m.d.sync += [out_data.eq(Cat(xi, xq)), out_full.eq(1)]
            with m.If(first):
                m.d.sync += [snap_req.eq(0), snap_on.eq(1), snap_len.eq(snap_next),
                             snap_left.eq(Cat(Const(1022, 10), snap_next))]
            with m.Elif(snap_left == 0):
                m.d.sync += [snap_on.eq(0), trail.eq(1)]
            with m.Else():
                m.d.sync += snap_left.eq(snap_left - 1)
        with m.Elif(tick & (trail != 0)):
            with m.If(out_full & ~self.out_ready):
                m.d.sync += self.snap_lost.eq(1)
            with m.Switch(trail):
                with m.Case(1):
                    m.d.sync += out_data.eq(Cat(Const(SNAP_MAGIC, 8), snap_len))
                with m.Case(2):
                    m.d.sync += out_data.eq(n[:16])
                with m.Case(3):
                    m.d.sync += out_data.eq(n[16:])
            m.d.sync += [out_full.eq(1), trail.eq(trail + 1)]
        # After the snapshot's start, so a request in the cycle a snapshot begins is kept.
        with m.If(self.snap_valid):
            m.d.sync += [snap_req.eq(1), snap_next.eq(self.snap_value)]

        return m
