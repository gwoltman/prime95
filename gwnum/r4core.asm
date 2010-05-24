; Copyright 2009-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4 FFTs optimized for the Core 2 / Core i7
; architecture.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

xfft_type TEXTEQU <r4>
xarch	TEXTEQU	<CORE>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE r4.mac
INCLUDE r4pass1sc.mac
INCLUDE r4pass2.mac

EXTRN	pass1_aux_entry_point_return:PROC
EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC

_TEXT SEGMENT

;; Some small FFTs fit completely in a large L2 cache and thus do
;; not need to be prefetched.

PREFETCHING = 0

	xpass2gen 8
	xpass2gen 768
	xpass2gen 10
	xpass2gen 1280

	r4_pass1sc32 8K, 8, 4
	r4_pass1sc64 16K, 8, 4
	r4_pass1sc128 32K, 8, 4
	r4_pass1sc256 64K, 8, 4
	r4_pass1sc32 24K, 768, 4
	r4_pass1sc64 48K, 768, 4
	r4_pass1sc128 96K, 768, 4
	r4_pass1sc256 192K, 768, 4
	r4_pass1sc32 32K, 10, 4
	r4_pass1sc64 64K, 10, 4
	r4_pass1sc32 40K, 1280, 4
	r4_pass1sc64 80K, 1280, 4

	r4_pass1sc32ac 8K, 8, 4
	r4_pass1sc64ac 16K, 8, 4
	r4_pass1sc128ac 32K, 8, 4
	r4_pass1sc256ac 64K, 8, 4
	r4_pass1sc32ac 24K, 768, 4
	r4_pass1sc64ac 48K, 768, 4
	r4_pass1sc128ac 96K, 768, 4
	r4_pass1sc256ac 192K, 768, 4
	r4_pass1sc32ac 32K, 10, 4
	r4_pass1sc64ac 64K, 10, 4
	r4_pass1sc32ac 40K, 1280, 4
	r4_pass1sc64ac 80K, 1280, 4

;; Define names and features for this architecture and FFT type

PREFETCHING = 64
TLB_PRIMING = 0

;; Generate pass 2 routines optmized for this architecture

	xpass2gen 8
	xpass2gen 768
	xpass2gen 10
	xpass2gen 1280

;; Routines for many FFT sizes

	r4_pass1sc32 8K, 8, 4
	r4_pass1sc64 16K, 8, 4
	r4_pass1sc128 32K, 8, 4
	r4_pass1sc256 64K, 8, 4
	r4_pass1sc32 24K, 768, 4
	r4_pass1sc64 48K, 768, 4
	r4_pass1sc128 96K, 768, 4
	r4_pass1sc256 192K, 768, 4
	r4_pass1sc32 32K, 10, 4
	r4_pass1sc64 64K, 10, 4
	r4_pass1sc32 40K, 1280, 4
	r4_pass1sc64 80K, 1280, 4

	r4_pass1sc32ac 8K, 8, 4
	r4_pass1sc64ac 16K, 8, 4
	r4_pass1sc128ac 32K, 8, 4
	r4_pass1sc256ac 64K, 8, 4
	r4_pass1sc32ac 24K, 768, 4
	r4_pass1sc64ac 48K, 768, 4
	r4_pass1sc128ac 96K, 768, 4
	r4_pass1sc256ac 192K, 768, 4
	r4_pass1sc32ac 32K, 10, 4
	r4_pass1sc64ac 64K, 10, 4
	r4_pass1sc32ac 40K, 1280, 4
	r4_pass1sc64ac 80K, 1280, 4

_TEXT	ENDS
END
