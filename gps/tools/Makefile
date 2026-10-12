# SPDX-License-Identifier: BSD-3-Clause
#
# hackrf_gps: needs libhackrf, libusb-1.0 and fftw3f (found by pkg-config).

CFLAGS = -std=gnu11 -O2 -Wall -Wextra -Werror $(shell pkg-config --cflags libhackrf libusb-1.0)
LIBS   = $(shell pkg-config --libs libhackrf libusb-1.0 fftw3f) -lm -lpthread

all: build/hackrf_gps

build/hackrf_gps: hackrf_gps.c
	@mkdir -p build
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

clean:
	rm -rf build

.PHONY: all clean
