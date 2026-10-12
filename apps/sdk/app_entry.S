/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * The one entry point the firmware calls, unprivileged on the app's stack:
 * r0 the operation, r1 its argument. app_main's result goes back with
 * SVC 0, the only way out; the firmware set lr to an address that faults.
 */
	.syntax unified
	.thumb
	.section .text.app_entry, "ax"
	.global app_entry
	.thumb_func
app_entry:
	bl app_main
	svc #0
1:	b 1b
