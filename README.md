# hackrf-gpsdo

> **Disclaimer.** This is independent work. It is not a product of, endorsed by, or supported by
> Great Scott Gadgets. It is provided "as is", without warranty of any kind, express or implied,
> including any warranty of merchantability, fitness for a particular purpose or accuracy. You use
> it entirely at your own risk. Installing modified firmware, loading custom FPGA images and
> changing the board's clock synthesizer can leave a device unusable until it is recovered, can
> affect your warranty, and can produce incorrect frequencies or unintended emissions. The authors
> accept no liability for any damage or loss arising from its use. You are responsible for operating
> your equipment in accordance with the laws and regulations that apply where you are.

A GPS-disciplined frequency reference for the HackRF Pro. A GPS L1 C/A receiver runs on the board
itself: the FPGA correlates, the LPC4330's Cortex-M4 searches, tracks, navigates and measures the
board's 25 MHz TCXO against the satellites' carriers, and it steers the Si5351C's PLL B so that
CLKOUT (10 MHz) and the radio's LO references follow GPS. The host only loads the image and the app,
sets up the radio and prints what the receiver reports.

## What it does

- **Gateware:** a receive-only FPGA image (decimation to 4 Msps) with a seven-channel GPS/SBAS
  correlator. For each channel and code period it sends the early, prompt and late sums with the
  exact sample timing, and on request it sends snapshots of raw samples for the search.
- **The M4 app:** the satellite search (FFT-based, from the almanac when there is one), the tracking
  loops (FLL pull-in, Costas PLL and DLL), the navigation message, position fixes and a survey, the
  frequency estimate, and the steering loop. It runs as an unprivileged app under the MPU, loaded
  into RAM beside the FPGA image.
- **The firmware** (a fork of the HackRF firmware, the `firmware/hackrf` submodule) loads FPGA
  images staged in RAM, runs the app, and gives it a queue of FPGA register writes, messages to and
  from the host, and control of PLL B's feedback ratio. Nothing is written to flash at run time.
- **`hackrf_gps`**, the host tool, switches the board to GPS, passes the settings and a warm-start
  state, prints the receiver's status each second, and on exit puts PLL B back and returns the
  board to the stock receive image.

PLL A, which clocks the ADC and the samples, stays on the raw TCXO; only PLL B is steered, within
5 ppm of nominal. CLKOUT and anything tuned with the HackRF's synthesizers then carry GPS accuracy
while the TCXO drifts underneath.

This is a frequency reference, not a time reference: there is no pulse per second and no UTC. With
a known position one satellite is enough for the frequency estimate; four are needed for a fix.

**What to expect.** With an active antenna and a clear sky, the receiver's own measure of the
steered CLKOUT is about 0.06–0.08 ppb rms each second, with an Allan deviation near 5e-11 at 1 s and
below 1e-12 at 1,000 s. Against an independent GPSDO-referenced receiver over two hours, the mean
offset was within 0.2 ppb and the scatter about 0.3 ppb rms each second (that measurement's own
noise included). From a cold start with the position given, the first satellites are tracked in
about 15 s, the first estimate follows about 30 s later (once their ephemeris is read), and the
steering catches the TCXO about 40 s after that. A warm start (a saved almanac and position) tracks
four satellites in about 10 s.

To check CLKOUT independently, record it (or one of its harmonics) with another SDR locked to a
GPSDO, for example a USRP B2xx with its GPSDO option, and compare its phase over time.

## Repository layout

| Path | Contents |
|---|---|
| `firmware/hackrf` | submodule: the HackRF firmware and host tools with RAM image loading and app support |
| `fpga-build/` | the gateware (Amaranth): the GPS image `top/gps_rx.py`, the correlator `dsp/gps_engine.py` and its command port `dsp/gps_cmd.py`, the models they are tested against, the image gate (`check_image.py`) and packer (`pack.py`) |
| `fpga-build/packs/gps_rx.pack` | the built, gated GPS image |
| `apps/` | the app SDK (`sdk/`), the GPS app (`gps/gps_app.c`) and sandbox test apps |
| `gps/app/` | the app's loops, navigation, frequency estimate and search, shared with the host tests |
| `gps/tools/hackrf_gps.c` | the host tool |
| `gps/*.py`, `gps/test_*.py`, `fpga-build/test/` | reference code and the test suites |
| `prebuilt/` | the firmware (`hackrf_usb.bin` for flash, `hackrf_usb_fpga.dfu` for RAM), `gps_app.bin`, and `SHA256SUMS` |
| `build.sh` | builds everything from source and compares it with the committed binaries |

## Requirements

- A HackRF Pro.
- A GPS L1 antenna with a view of the sky. An active antenna needs the HackRF's antenna port power
  (`-p 1`). A passive antenna runs with port power off (`-p 0`); with no amplifier at the antenna,
  only strong satellites are seen, so a directional antenna or a saved position helps.
- A Linux host. To build: `git`, `cmake`, `make`, a C compiler, `arm-none-eabi-gcc` with newlib,
  `pkg-config`, the development files of `libusb-1.0`, `fftw3` (single precision) and `libhackrf`,
  and Python 3.12 or newer with `venv`. To install the firmware: `dfu-util` or `hackrf_spiflash`.

## Getting the source

```bash
git clone --recurse-submodules https://github.com/willcode/hackrf-gpsdo.git
```

In an existing clone, `git submodule update --init --recursive` fetches the firmware and its
libopencm3 submodule.

## Building

For everyday use, build the app and the host tool in place, and use the committed GPS image and the
prebuilt firmware:

```bash
make -C apps gps-app
```
```bash
make -C gps/tools
```

`gps/tools/build/hackrf_gps` then finds `fpga-build/packs/gps_rx.pack` and `apps/build/gps_app.bin`
by its own location.

To build everything from source and check it against the committed binaries:

```bash
./build.sh out
```

`out` must be empty or absent. The script clones the repository at `HEAD` into `out/src`, builds
the firmware, the fork's libhackrf and tools, the app, `hackrf_gps` and the GPS image, and reports
whether each matches `prebuilt/` and `fpga-build/packs/gps_rx.pack`. `--tests` also runs the test
suites; `--venv DIR` reuses a Python environment that has `fpga-build/requirements-lock.txt`
installed; `--only firmware,host,app,tools,gateware` builds a subset. The gateware environment is
made with `pip` from the lock file and needs network access; the pinned YoWASP Yosys and nextpnr
versions make the image build deterministic. A full build takes a few minutes.

The firmware build is deterministic for a given compiler; the prebuilt files were built with
`arm-none-eabi-gcc` 15.2.0.

## Installing the firmware

The receiver needs the firmware from `firmware/hackrf` on the board. The flash firmware keeps the
stock receive, transmit and sweep functions; the RAM firmware serves the receiver only.

**From RAM, nothing written:** put the board in DFU mode (hold DFU, press and release RESET, release
DFU) and load the RAM firmware:

```bash
dfu-util --device 1fc9:000c --alt 0 --download prebuilt/hackrf_usb_fpga.dfu
```

It runs until the board is reset or unpowered. A RAM firmware has no FPGA images of its own in
flash, so at exit `hackrf_gps` cannot return the board to the stock receive image; reset the board
for other use.

**In flash:** first save the current flash, then write the new firmware and reset:

```bash
hackrf_spiflash -r hackrf-pro-flash-backup.bin
```
```bash
hackrf_spiflash -w prebuilt/hackrf_usb.bin -R
```

`hackrf_spiflash -w` erases the whole flash before writing. To go back, write the backup the same
way. If the board does not start, the boot ROM's DFU mode always works: load a RAM firmware as above,
then write the backup.

## Running

Turn CLKOUT on (the steering moves it whether it is on or not):

```bash
hackrf_clock -o 1
```

Then, with an active antenna:

```bash
gps/tools/build/hackrf_gps -p 1 -t 0 -D 20
```

This powers the antenna, runs until Ctrl-C, and turns the steering loop on 20 s into tracking. Keep
`hackrf_gps` running: the app runs only while the host stays attached, and a USB reconfiguration
unloads it. On exit, PLL B goes back to nominal.

| Option | Meaning | Default |
|---|---|---|
| `-t S` | seconds of tracking; 0 runs until Ctrl-C | 60 |
| `-D S` | turn the steering loop on at second S of tracking | off |
| `-p 0\|1`, `-a 0\|1` | antenna port power; RF amplifier | 0; 1 |
| `-l DB`, `-g DB` | LNA gain (0–40, 8 dB steps); VGA gain (0–62, 2 dB steps) | 32; 30 |
| `-P X,Y,Z` | hold this ECEF position, in meters | survey |
| `-P RUN.fix` | hold the mean of a previous run's fixes | |
| `-P GRID[,HEIGHT]` | start from a Maidenhead square's center (4–10 characters); the survey replaces it | |
| `-E HEIGHT` | antenna height above the WGS 84 ellipsoid, held until the first full fix | 0 |
| `-W FILE`, `-W none` | the warm-start state, read at the start and written at the end | `gps/data/hackrf_gps.state` |
| `-o FILE` | record the stream to FILE, the fixes to FILE.fix and the estimates to FILE.freq | none |
| `-G FILE`, `-A FILE` | the GPS image and the app to load | `fpga-build/packs/gps_rx.pack`, `apps/build/gps_app.bin` |

`hackrf_gps --help` lists every option. With the prebuilt app, add `-A prebuilt/gps_app.bin`.

**What it prints:** the switch to GPS and the tuner's offset; the first search's progress; then each
second the channels (PRN, carrier, prompt power, lock), a fix line (satellites, GDOP, GPS time of
week) and a frequency line:

```
frequency: TCXO -97.426 ppb from 4 satellites (rms 0.017 ppb), position surveyed (60 fixes); PLL B steered +97.380 ppb, CLKOUT -0.046 ppb
```

"TCXO" is the board oscillator's error, "PLL B steered" the correction in effect, and "CLKOUT" their
sum: the steered output's remaining error by the receiver's own measure.

**The first run** without a saved state searches cold and surveys its position from 60 fixes; give
a grid square (`-P FN20xr`) to start nearer. Later runs start warm from `gps/data/hackrf_gps.state`.
That file, and any `.fix` file, holds the antenna's position; `hackrf_gps` prints the position
nowhere else.

**Board selection:** `hackrf_gps` runs on the one HackRF Pro whose firmware supports apps, and
leaves any other HackRF alone (another HackRF clocked from its CLKOUT, for example). With two such
boards attached it refuses to start.

## Tests

```bash
./build.sh --tests out
```

runs the gateware simulations against their models (each correlator field, the command port, the
packer and the image gate) and the host tests (the app's loops, search and navigation against the
Python references, and the app itself in a host harness). Some host tests replay recordings in
`gps/data/` and are skipped when those are absent. The sandbox test apps (`make -C apps check`)
run on a board with the firmware installed and need `hackrf_app_test` from the fork's host tools.

## License

The code in this repository is under the BSD 3-Clause license (`LICENSE`). Files from or based on
the HackRF gateware keep their Great Scott Gadgets copyright notices under the same license. The
firmware in `firmware/hackrf` is a fork of the HackRF repository and keeps its licenses: GPL 2.0 or
later for the firmware and tools, BSD 3-Clause for libhackrf.

HackRF is a trademark of Great Scott Gadgets.
