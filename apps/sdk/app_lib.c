/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * The few C library functions the compiler may call on its own (structure
 * copies and clears); apps link no C library. Kept (used) under -flto, where
 * the calls the compiler makes itself appear only after the link has dropped
 * what nothing named.
 */

#include <stddef.h>

__attribute__((used)) void* memset(void* dst, int c, size_t n)
{
	unsigned char* d = dst;
	while (n--) {
		*d++ = (unsigned char) c;
	}
	return dst;
}

__attribute__((used)) void* memcpy(void* dst, const void* src, size_t n)
{
	unsigned char* d = dst;
	const unsigned char* s = src;
	while (n--) {
		*d++ = *s++;
	}
	return dst;
}
