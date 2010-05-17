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

	r4delay_pass1sc320 1440K, 4608, 4
	r4delay_pass1sc320 1920K, 6144, 4
	r4delay_pass1sc320 2400K, 7680, 4
	r4delay_pass1sc320 2560K, 13, 4
	r4delay_pass1sc320 3200K, 10240, 4
	r4delay_pass1sc320 4000K, 12800, 4

	r4delay_pass1sc384 1728K, 4608, 4
	r4delay_pass1sc384 2304K, 6144, 4
	r4delay_pass1sc384 2880K, 7680, 4
	r4delay_pass1sc384 3M, 13, 4
	r4delay_pass1sc384 3840K, 10240, 4
	r4delay_pass1sc384 4800K, 12800, 4

	r4delay_pass1sc448 2016K, 4608, 4
	r4delay_pass1sc448 2688K, 6144, 4
	r4delay_pass1sc448 3360K, 7680, 4
	r4delay_pass1sc448 3584K, 13, 4
	r4delay_pass1sc448 4480K, 10240, 4
	r4delay_pass1sc448 5600K, 12800, 4

	r4delay_pass1sc512 2304K, 4608, 4
	r4delay_pass1sc512 3M, 6144, 4
	r4delay_pass1sc512 3840K, 7680, 4
	r4delay_pass1sc512 4M, 13, 4
	r4delay_pass1sc512 4M, 13, 2
	r4delay_pass1sc512 4M, 13, 1
	r4delay_pass1sc512 5M, 10240, 4
	r4delay_pass1sc512 6400K, 12800, 4

	r4delay_pass1sc640 2880K, 4608, 4
	r4delay_pass1sc640 3840K, 6144, 4
	r4delay_pass1sc640 4800K, 7680, 4
	r4delay_pass1sc640 5M, 13, 4
	r4delay_pass1sc640 6400K, 10240, 4
	r4delay_pass1sc640 8000K, 12800, 4

	r4delay_pass1sc768 3456K, 4608, 4
	r4delay_pass1sc768 4608K, 6144, 4
	r4delay_pass1sc768 5760K, 7680, 4
	r4delay_pass1sc768 6M, 13, 4
	r4delay_pass1sc768 7680K, 10240, 4
	r4delay_pass1sc768 9600K, 12800, 4

	r4delay_pass1sc896 4032K, 4608, 4
	r4delay_pass1sc896 5376K, 6144, 4
	r4delay_pass1sc896 6720K, 7680, 4
	r4delay_pass1sc896 7M, 13, 4
	r4delay_pass1sc896 8960K, 10240, 4
	r4delay_pass1sc896 11200K, 12800, 4

	r4delay_pass1sc1024 4608K, 4608, 4
	r4delay_pass1sc1024 6M, 6144, 4
	r4delay_pass1sc1024 7680K, 7680, 4
	r4delay_pass1sc1024 8M, 13, 4
	r4delay_pass1sc1024 10M, 10240, 4
	r4delay_pass1sc1024 12800K, 12800, 4

	r4delay_pass1sc2048 9M, 4608, 4
	r4delay_pass1sc2048 12M, 6144, 4
	r4delay_pass1sc2048 15M, 7680, 4
	r4delay_pass1sc2048 16M, 13, 4
	r4delay_pass1sc2048 16M, 13, 2
	r4delay_pass1sc2048 16M, 13, 1
	r4delay_pass1sc2048 20M, 10240, 4
	r4delay_pass1sc2048 20M, 10240, 2
	r4delay_pass1sc2048 25M, 12800, 4

	r4delay_pass1sc2560 20M, 13, 4
	r4delay_pass1sc2560 20M, 13, 2

	r4delay_pass1sc4096 18M, 4608, 4
	r4delay_pass1sc4096 18M, 4608, 2
	r4delay_pass1sc4096 24M, 6144, 2
	r4delay_pass1sc4096 30M, 7680, 2
	r4delay_pass1sc4096 32M, 13, 2
	;r4delay_pass1sc4096 40M, 10240, 2
	;r4delay_pass1sc4096 50M, 12800, 2

	; The all-complex 13 levels variants (4608, 6144, 7680, 8192, 10240, 12800)

	r4delay_pass1sc256ac 1152K, 4608, 4
	r4delay_pass1sc256ac 1536K, 6144, 4
	r4delay_pass1sc256ac 1920K, 7680, 4
	r4delay_pass1sc256ac 2560K, 10240, 4
	r4delay_pass1sc256ac 3200K, 12800, 4

	r4delay_pass1sc512ac 3M, 6144, 4
	r4delay_pass1sc512ac 4M, 13, 4
	r4delay_pass1sc512ac 4M, 13, 2
	r4delay_pass1sc512ac 4M, 13, 1

	r4delay_pass1sc2048ac 16M, 13, 4

_TEXT	ENDS
END
