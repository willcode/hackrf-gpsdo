# SPDX-License-Identifier: BSD-3-Clause
#
# Runtime apps for the HackRF Pro's M4 (firmware with app support).
# make test-apps    builds build/test_<n>.bin, one per sandbox test
# make gps-app      builds build/gps_app.bin
# make check        runs each on the device (needs the RAM_FPGA build
#                   with app support running, and hackrf_app_test from
#                   its host tools; TOOL=path selects another)

CC      = arm-none-eabi-gcc
OBJCOPY = arm-none-eabi-objcopy
CFLAGS  = -mcpu=cortex-m4 -mthumb -mfloat-abi=hard -mfpu=fpv4-sp-d16 -Os -g \
          -ffreestanding -fno-common -ffunction-sections -fdata-sections \
          -Wall -Wextra -Werror -Isdk
LDFLAGS = -nostdlib -Wl,--gc-sections -Wl,--no-warn-rwx-segments -T sdk/app.ld
TOOL   ?= hackrf_app_test

# test number : the cause the firmware must report (app_sandbox.h);
# 0 runs, and is then stopped (cause 1).
TESTS = 0:1 1:2 2:2 3:6 4:7 5:2 6:2 7:4 8:10 9:2 10:3 11:1 12:1

build/test_%.elf: test/test_apps.c sdk/app_entry.S sdk/app.ld sdk/app.h
	@mkdir -p build
	$(CC) $(CFLAGS) -DTEST=$* $(LDFLAGS) -o $@ sdk/app_entry.S test/test_apps.c -lgcc

build/%.bin: build/%.elf
	$(OBJCOPY) -O binary $< $@

# The GPS app: the tracking loops, navigation, the frequency estimate and acquisition for the GPS
# image (gps/app/gps_loops.c, gps_nav.c, gps_freq.c and gps_acq.c are shared with the host's
# tests). Budget 2 ms per call. Linked whole (-flto) and once-called functions kept apart: about
# 760 bytes fewer, in an image region of 28 KB. No errno (there is none here), so a single-precision
# square root is the FPU's instruction rather than a library call.
GPS_CFLAGS = -flto -fno-inline-functions-called-once -fno-math-errno \
	'-DGPS_ACQ_STATE=__attribute__((section(".bss.appdata")))'
build/gps_app.elf: gps/gps_app.c ../gps/app/gps_loops.c ../gps/app/gps_loops.h \
		../gps/app/gps_nav.c ../gps/app/gps_nav.h ../gps/app/gps_freq.c ../gps/app/gps_freq.h \
		../gps/app/gps_acq.c ../gps/app/gps_acq.h sdk/app_lib.c sdk/app_entry.S sdk/app.ld sdk/app.h
	@mkdir -p build
	$(CC) $(CFLAGS) $(GPS_CFLAGS) -I../gps/app $(LDFLAGS) -Wl,--defsym=APP_BUDGET_US=2000 -o $@ \
		sdk/app_entry.S gps/gps_app.c ../gps/app/gps_loops.c ../gps/app/gps_nav.c \
		../gps/app/gps_freq.c ../gps/app/gps_acq.c sdk/app_lib.c -lgcc

gps-app: build/gps_app.bin

test-apps: $(foreach t,$(TESTS),build/test_$(firstword $(subst :, ,$(t))).bin)

check: test-apps
	@for t in $(TESTS); do \
		n=$${t%%:*}; c=$${t##*:}; \
		echo "== test $$n, expect cause $$c"; \
		if [ $$n = 0 ]; then \
			$(TOOL) -l build/test_0.bin -w 2 && $(TOOL) -m 01020304 -r 0404030201 && \
			$(TOOL) -s -e 1 || exit 1; \
		elif [ $$c = 1 ]; then \
			$(TOOL) -l build/test_$$n.bin -w 2 && $(TOOL) -s -e 1 || exit 1; \
		else \
			$(TOOL) -l build/test_$$n.bin -w 3 -e $$c || exit 1; \
		fi; \
	done

clean:
	rm -rf build

.PHONY: test-apps gps-app check clean
.SECONDARY:
