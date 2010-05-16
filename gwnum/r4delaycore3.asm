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

	xpass2gen 2304
	xpass2gen 3072
	xpass2gen 3840
	xpass2gen 12
	xpass2gen 5120
	xpass2gen 6400

;; Routines for many FFT sizes

	; The 12 levels variants (2304, 3072, 3840, 4096, 5120, 6400)

	r4delay_pass1sc256 576K, 2304, 4
	r4delay_pass1sc256 768K, 3072, 4
	r4delay_pass1sc256 960K, 3840, 4
	r4delay_pass1sc256 1M, 12, 4
	r4delay_pass1sc256 1280K, 5120, 4
	r4delay_pass1sc256 1600K, 6400, 4

	r4delay_pass1sc320 720K, 2304, 4
	r4delay_pass1sc320 960K, 3072, 4
	r4delay_pass1sc320 1200K, 3840, 4
	r4delay_pass1sc320 1280K, 12, 4
	r4delay_pass1sc320 1600K, 5120, 4
	r4delay_pass1sc320 2000K, 6400, 4

	r4delay_pass1sc384 864K, 2304, 4
	r4delay_pass1sc384 1152K, 3072, 4
	r4delay_pass1sc384 1440K, 3840, 4
	r4delay_pass1sc384 1536K, 12, 4
	r4delay_pass1sc384 1920K, 5120, 4
	r4delay_pass1sc384 2400K, 6400, 4

	r4delay_pass1sc448 1008K, 2304, 4
	r4delay_pass1sc448 1344K, 3072, 4
	r4delay_pass1sc448 1680K, 3840, 4
	r4delay_pass1sc448 1792K, 12, 4
	r4delay_pass1sc448 2240K, 5120, 4
	r4delay_pass1sc448 2800K, 6400, 4

	r4delay_pass1sc512 1152K, 2304, 4
	r4delay_pass1sc512 1536K, 3072, 4
	r4delay_pass1sc512 1920K, 3840, 4
	r4delay_pass1sc512 2M, 12, 4
	r4delay_pass1sc512 2560K, 5120, 4
	r4delay_pass1sc512 3200K, 6400, 4

	r4delay_pass1sc640 1440K, 2304, 4
	r4delay_pass1sc640 1920K, 3072, 4
	r4delay_pass1sc640 2400K, 3840, 4
	r4delay_pass1sc640 2560K, 12, 4
	r4delay_pass1sc640 3200K, 5120, 4
	r4delay_pass1sc640 4000K, 6400, 4

	r4delay_pass1sc768 1728K, 2304, 4
	r4delay_pass1sc768 2304K, 3072, 4
	r4delay_pass1sc768 2880K, 3840, 4
	r4delay_pass1sc768 3M, 12, 4
	r4delay_pass1sc768 3840K, 5120, 4
	r4delay_pass1sc768 4800K, 6400, 4

	r4delay_pass1sc896 2016K, 2304, 4
	r4delay_pass1sc896 2688K, 3072, 4
	r4delay_pass1sc896 3360K, 3840, 4
	r4delay_pass1sc896 3584K, 12, 4
	r4delay_pass1sc896 4480K, 5120, 4
	r4delay_pass1sc896 5600K, 6400, 4

	r4delay_pass1sc1024 2304K, 2304, 4
	r4delay_pass1sc1024 3M, 3072, 4
	r4delay_pass1sc1024 3840K, 3840, 4
	r4delay_pass1sc1024 4M, 12, 4
	r4delay_pass1sc1024 5M, 5120, 4
	r4delay_pass1sc1024 6400K, 6400, 4

	r4delay_pass1sc1280 2880K, 2304, 4
	r4delay_pass1sc1280 3840K, 3072, 4
	r4delay_pass1sc1280 4800K, 3840, 4
	r4delay_pass1sc1280 5M, 12, 4
	r4delay_pass1sc1280 6400K, 5120, 4
	r4delay_pass1sc1280 8000K, 6400, 4

	r4delay_pass1sc1536 6M, 12, 4

	r4delay_pass1sc1792 7M, 12, 4

	r4delay_pass1sc2048 8M, 12, 4

	r4delay_pass1sc2560 10M, 12, 4

	r4delay_pass1sc3072 12M, 12, 4

	r4delay_pass1sc3584 14M, 12, 4

	r4delay_pass1sc4096 16M, 12, 4
	r4delay_pass1sc4096 16M, 12, 2

	; The 12 levels variants (2304, 3072, 3840, 4096, 5120, 6400)

	r4delay_pass1sc256ac 576K, 2304, 4
	r4delay_pass1sc256ac 1M, 12, 4
	r4delay_pass1sc384ac 864K, 2304, 4
	r4delay_pass1sc384ac 1152K, 3072, 4
	r4delay_pass1sc384ac 1536K, 12, 4
	r4delay_pass1sc512ac 1536K, 3072, 4
	r4delay_pass1sc512ac 1920K, 3840, 4
	r4delay_pass1sc512ac 2M, 12, 4
	r4delay_pass1sc640ac 2560K, 12, 4
	r4delay_pass1sc640ac 3200K, 5120, 4
	r4delay_pass1sc640ac 4000K, 6400, 4
	r4delay_pass1sc768ac 1728K, 2304, 4
	r4delay_pass1sc768ac 2304K, 3072, 4
	r4delay_pass1sc768ac 3M, 12, 4
	r4delay_pass1sc1024ac 3M, 3072, 4
	r4delay_pass1sc1024ac 4M, 12, 4
	r4delay_pass1sc1280ac 5M, 12, 4
	r4delay_pass1sc1536ac 6M, 12, 4
	r4delay_pass1sc2048ac 8M, 12, 4
	r4delay_pass1sc2560ac 10M, 12, 4
	r4delay_pass1sc3072ac 12M, 12, 4
	r4delay_pass1sc4096ac 16M, 12, 4
	r4delay_pass1sc5120ac 11520K, 2304, 4
	r4delay_pass1sc5120ac 20M, 12, 4

_TEXT	ENDS
END
