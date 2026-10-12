/*
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * A runtime app's image: linked to run where the firmware copies it
 * (APP_IMAGE_BASE in firmware/common/app_sandbox.h), header first. The
 * header's fields are written here, so the image needs no build step after
 * the link beyond objcopy. .data stays in place (the image is loaded at its
 * run address); the firmware zeroes .bss, the stack below the image and the
 * data region at 0x10088000. Variables marked APP_STATE (app.h) live in the
 * data region instead of after the image: 40 kB less the 256 bytes at its
 * top that carry host messages.
 *
 * An app may set its time budget per call with --defsym=APP_BUDGET_US=n.
 */

MEMORY
{
	app (rwx) : ORIGIN = 0x10019000, LENGTH = 28K
	data (rw) : ORIGIN = 0x10088000, LENGTH = 0xA000 - 256
}

ENTRY(app_entry)

SECTIONS
{
	.header : {
		LONG(0x50504148)                      /* magic "HAPP" */
		SHORT(1)                              /* API version */
		SHORT(32)                             /* header size */
		LONG(__image_end - ORIGIN(app))       /* image size */
		LONG(__bss_end - __image_end)         /* .bss size, from the image's end */
		LONG(app_entry | 1)                   /* entry, Thumb */
		LONG(DEFINED(APP_BUDGET_US) ? APP_BUDGET_US : 0)
		LONG(0)
		LONG(0)
	} > app

	.text : {
		*(.text.app_entry)
		*(.text*)
		*(.rodata*)
	} > app

	.data : {
		*(.data*)
		. = ALIGN(4);
	} > app

	__image_end = .;

	/* Before .bss, whose *(.bss*) would otherwise take it. */
	.appdata (NOLOAD) : {
		*(.bss.appdata*)
	} > data

	.bss (NOLOAD) : {
		__bss_start = .;
		*(.bss*)
		*(COMMON)
		. = ALIGN(4);
		__bss_end = .;
	} > app

	/DISCARD/ : { *(.ARM.exidx*) *(.comment) }
}
