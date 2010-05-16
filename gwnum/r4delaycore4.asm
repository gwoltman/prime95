; Copyright 2009-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These FFTs were split out of r4delaycore.asm because of MASM limitations.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

xfft_type TEXTEQU <r4delay>
xarch	TEXTEQU	<CORE>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE r4.mac
INCLUDE r4delaypass1sc.mac
INCLUDE r4delaypass2.mac

EXTRN	pass1_aux_entry_point_return:PROC
EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC

_TEXT SEGMENT

;; Define names and features for this architecture and FFT type

PREFETCHING = 64
TLB_PRIMING = 0

;; Generate pass 2 routines optmized for this architecture

	xpass2gen 4608
	xpass2gen 6144
	xpass2gen 7680
	xpass2gen 13
	xpass2gen 10240
	xpass2gen 12800

;; Routines for many FFT sizes

	; The 13 levels variants (4608, 6144, 7680, 8192, 10240, 12800)

	r4delay_pass1sc256 1152K, 4608, 4
	r4delay_pass1sc256 1536K, 6144, 4
	r4delay_pass1sc256 1920K, 7680, 4
	r4delay_pass1sc256 2M, 13, 4
	r4delay_pass1sc256 2560K, 10240, 4
	r4delay_pass1sc256 3200K, 12800, 4
	r4delay_pass1sc512 4M, 13, 4
	r4delay_pass1sc2560 20M, 13, 2
	r4delay_pass1sc4096 18M, 4608, 2

	; The all-complex 13 levels variants (4608, 6144, 7680, 8192, 10240, 12800)

	r4delay_pass1sc256ac 1152K, 4608, 4
	r4delay_pass1sc256ac 1536K, 6144, 4
	r4delay_pass1sc256ac 1920K, 7680, 4
	r4delay_pass1sc256ac 2560K, 10240, 4
	r4delay_pass1sc256ac 3200K, 12800, 4
	r4delay_pass1sc512ac 3M, 6144, 4
	r4delay_pass1sc512ac 4M, 13, 4

_TEXT	ENDS
END
