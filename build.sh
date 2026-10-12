#!/bin/sh
# SPDX-License-Identifier: BSD-3-Clause
#
# Build everything from committed source into a new directory: the firmware, the host library and
# tools, the M4 app, hackrf_gps and the GPS gateware image, then compare the results with the
# committed ones.
#
# usage: build.sh [--venv DIR] [--tests] [--only STEPS] OUT
#
# OUT must be empty or absent. Nothing outside OUT is written. The sources are a clone of this
# repository at HEAD (uncommitted changes are not built), with the firmware submodule and its
# libopencm3 submodule taken from this checkout:
#   OUT/src              the clone
#   OUT/build/firmware   hackrf_usb.bin (flash image) and hackrf_usb_fpga.dfu (RAM firmware)
#   OUT/build/host       libhackrf and the hackrf tools of the firmware fork
#   OUT/build/app        gps_app.bin
#   OUT/build/tools      hackrf_gps
#   OUT/build/gateware   gps_rx.pack and its build directory
#   OUT/logs             one log per step
#
# STEPS is a comma-separated list of firmware, host, app, tools and gateware (default: all). --venv
# uses an existing gateware environment instead of making OUT/venv from
# fpga-build/requirements-lock.txt, which needs network access. --tests also runs the gateware and
# host test suites. The firmware and the app are compared with prebuilt/SHA256SUMS, and the pack with
# the committed fpga-build/packs/gps_rx.pack.
set -eu

here=$(cd "$(dirname "$0")" && pwd)
venv= tests= steps=firmware,host,app,tools,gateware
while [ $# -gt 1 ]; do
    case $1 in
        --venv) venv=$(cd "$2" && pwd); shift 2 ;;
        --tests) tests=1; shift ;;
        --only) steps=$2; shift 2 ;;
        *) break ;;
    esac
done
[ $# -eq 1 ] || { sed -n '8p' "$0" | sed 's/^# //'; exit 2; }
mkdir -p "$1"
out=$(cd "$1" && pwd)
[ -z "$(ls -A "$out")" ] || { echo "build.sh: $out is not empty"; exit 2; }
want() { case ",$steps," in *",$1,"*) return 0 ;; esac; return 1; }

src=$out/src build=$out/build logs=$out/logs
mkdir -p "$build" "$logs"
step() { echo "== $1 (log: $logs/$2.log)"; }
fail() { echo "build.sh: $1 failed; see $logs/$2.log"; tail -5 "$logs/$2.log"; exit 1; }
# Compare a built file with its line in prebuilt/SHA256SUMS.
check() {
    want_sum=$(awk -v f="$2" '$2 == f { print $1 }' "$src/prebuilt/SHA256SUMS")
    got_sum=$(sha256sum "$1" | cut -d' ' -f1)
    echo "   $got_sum  $2"
    if [ "$want_sum" = "$got_sum" ]; then
        echo "   $2 is identical to prebuilt/$2"
    else
        echo "   $2 differs from prebuilt/$2"
    fi
}

step "sources" sources
fw=$here/firmware/hackrf
[ -e "$fw/.git" ] || { echo "build.sh: firmware/hackrf is not checked out; run: git submodule update --init --recursive"; exit 2; }
{
    git clone "$here" "$src" &&
    git -C "$src" submodule init &&
    git -C "$src" config submodule.firmware/hackrf.url "$fw" &&
    git -C "$src" -c protocol.file.allow=always submodule update &&
    if [ -e "$fw/firmware/libopencm3/.git" ]; then
        git -C "$src/firmware/hackrf" submodule init firmware/libopencm3 &&
        git -C "$src/firmware/hackrf" config submodule.firmware/libopencm3.url "$fw/firmware/libopencm3"
    fi &&
    git -C "$src/firmware/hackrf" -c protocol.file.allow=always submodule update --init firmware/libopencm3
} >"$logs/sources.log" 2>&1 || fail sources sources
echo "   $(git -C "$src" rev-parse --short HEAD), firmware $(git -C "$src/firmware/hackrf" describe --always --dirty)"

if want firmware; then
    step "firmware" firmware
    {
        cmake -S "$src/firmware/hackrf/firmware/hackrf_usb" -B "$build/firmware" -DBOARD=PRALINE &&
        cmake --build "$build/firmware" -j6 --target hackrf_usb.bin hackrf_usb_fpga.dfu
    } >"$logs/firmware.log" 2>&1 || fail firmware firmware
    check "$build/firmware/hackrf_usb.bin" hackrf_usb.bin
    check "$build/firmware/hackrf_usb_fpga.dfu" hackrf_usb_fpga.dfu
fi

if want host; then
    step "host library and tools" host
    {
        cmake -S "$src/firmware/hackrf/host" -B "$build/host" &&
        cmake --build "$build/host" -j6
    } >"$logs/host.log" 2>&1 || fail host host
fi

if want app; then
    step "M4 app" app
    {
        make -C "$src/apps" gps-app &&
        mkdir -p "$build/app" && cp "$src/apps/build/gps_app.bin" "$build/app/"
    } >"$logs/app.log" 2>&1 || fail app app
    check "$build/app/gps_app.bin" gps_app.bin
fi

if want tools; then
    step "hackrf_gps" tools
    {
        make -C "$src/gps/tools" &&
        mkdir -p "$build/tools" && cp "$src/gps/tools/build/hackrf_gps" "$build/tools/"
    } >"$logs/tools.log" 2>&1 || fail hackrf_gps tools
fi

if want gateware || [ -n "$tests" ]; then
    if [ -z "$venv" ]; then
        step "gateware environment" venv
        venv=$out/venv
        {
            python3 -m venv "$venv" &&
            "$venv/bin/pip" install -r "$src/fpga-build/requirements-lock.txt"
        } >"$logs/venv.log" 2>&1 || fail "gateware environment" venv
    fi
fi

if want gateware; then
    step "GPS image" gateware
    gw=$src/fpga-build
    g=$build/gateware/build-gps
    (
        export PATH="$venv/bin:$PATH" PYTHONPATH="$gw"
        export YOSYS=yowasp-yosys NEXTPNR_ICE40=yowasp-nextpnr-ice40 ICEPACK=yowasp-icepack
        mkdir -p "$build/gateware" && cd "$gw" &&
        rm -rf build-gps && mkdir -p build-gps && touch build-gps/.build-start &&
        python top/gps_rx.py &&
        python record_versions.py --since build-gps/.build-start build-gps 14_gps &&
        python check_image.py build-gps 14_gps &&
        python pack.py --gate build-gps/14_gps.gate.json build-gps/14_gps.bin "$build/gateware/gps_rx.pack" &&
        cp -r build-gps "$g"
    ) >"$logs/gateware.log" 2>&1 || fail "GPS image" gateware
    grep -h "Max frequency for clock" "$g/14_gps.tim" | tail -3 | sed 's/^.*clock /   /'
    echo "   $(sha256sum "$build/gateware/gps_rx.pack" | cut -d' ' -f1)  gps_rx.pack"
    if cmp -s "$build/gateware/gps_rx.pack" "$src/fpga-build/packs/gps_rx.pack"; then
        echo "   gps_rx.pack is identical to fpga-build/packs/gps_rx.pack"
    else
        echo "   gps_rx.pack differs from fpga-build/packs/gps_rx.pack"
    fi
fi

if [ -n "$tests" ]; then
    step "gateware tests" tests-gateware
    (cd "$src/fpga-build" && PYTHONPATH=. "$venv/bin/python" -m unittest discover -s test) \
        >"$logs/tests-gateware.log" 2>&1 || fail "gateware tests" tests-gateware
    grep -hE "^(Ran|OK|FAILED)" "$logs/tests-gateware.log" | sed 's/^/   /'
    step "host tests" tests-host
    (cd "$src/gps" && "$venv/bin/python" -m unittest discover -p 'test_*.py') \
        >"$logs/tests-host.log" 2>&1 || fail "host tests" tests-host
    grep -hE "^(Ran|OK|FAILED)" "$logs/tests-host.log" | sed 's/^/   /'
fi

echo "== done: $out"
