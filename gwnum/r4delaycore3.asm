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

;; Generate pass 2 routines optimized for this architecture

	xpass2gen 2304
	xpass2gen 3072
	xpass2gen 3840

;; Routines for many FFT sizes

	; The 12 levels variants (2304, 3072, 3840)

	r4delay_pass1sc256 576K, 2304, 4
	r4delay_pass1sc256 768K, 3072, 4
	r4delay_pass1sc256 960K, 3840, 4

	r4delay_pass1sc320 720K, 2304, 4
	r4delay_pass1sc320 960K, 3072, 4
	r4delay_pass1sc320 1200K, 3840, 4

	r4delay_pass1sc384 864K, 2304, 4
	r4delay_pass1sc384 1152K, 3072, 4
	r4delay_pass1sc384 1440K, 3840, 4

	r4delay_pass1sc448 1008K, 2304, 4
	r4delay_pass1sc448 1344K, 3072, 4
	r4delay_pass1sc448 1680K, 3840, 4

	r4delay_pass1sc512 1152K, 2304, 4
	r4delay_pass1sc512 1536K, 3072, 4
	r4delay_pass1sc512 1920K, 3840, 4

	r4delay_pass1sc640 1440K, 2304, 4
	r4delay_pass1sc640 1920K, 3072, 4
	r4delay_pass1sc640 2400K, 3840, 4

	r4delay_pass1sc768 1728K, 2304, 4
	r4delay_pass1sc768 2304K, 3072, 4
	r4delay_pass1sc768 2880K, 3840, 4

	r4delay_pass1sc896 2016K, 2304, 4
	r4delay_pass1sc896 2688K, 3072, 4
	r4delay_pass1sc896 3360K, 3840, 4

	r4delay_pass1sc1024 2304K, 2304, 4
	r4delay_pass1sc1024 3M, 3072, 4
	r4delay_pass1sc1024 3840K, 3840, 4

	r4delay_pass1sc1280 2880K, 2304, 4
	r4delay_pass1sc1280 3840K, 3072, 4
	r4delay_pass1sc1280 4800K, 3840, 4

	; The all-complex 12 levels variants (2304, 3072, 3840)

	r4delay_pass1sc256ac 576K, 2304, 4
	r4delay_pass1sc256ac 768K, 3072, 4
	r4delay_pass1sc256ac 960K, 3840, 4

	r4delay_pass1sc384ac 864K, 2304, 4
	r4delay_pass1sc384ac 1152K, 3072, 4
	r4delay_pass1sc384ac 1440K, 3840, 4

	r4delay_pass1sc512ac 1152K, 2304, 4
	r4delay_pass1sc512ac 1536K, 3072, 4
	r4delay_pass1sc512ac 1920K, 3840, 4

	r4delay_pass1sc640ac 1440K, 2304, 4
	r4delay_pass1sc640ac 1920K, 3072, 4
	r4delay_pass1sc640ac 2400K, 3840, 4

	r4delay_pass1sc768ac 1728K, 2304, 4
	r4delay_pass1sc768ac 2304K, 3072, 4
	r4delay_pass1sc768ac 2880K, 3840, 4

	r4delay_pass1sc1024ac 2304K, 2304, 4
	r4delay_pass1sc1024ac 3M, 3072, 4
	r4delay_pass1sc1024ac 3840K, 3840, 4

	r4delay_pass1sc1280ac 2880K, 2304, 4
	r4delay_pass1sc1280ac 3840K, 3072, 4
	r4delay_pass1sc1280ac 4800K, 3840, 4

	r4delay_pass1sc5120ac 11520K, 2304, 4

_TEXT	ENDS
END
