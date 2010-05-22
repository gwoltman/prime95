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

;; Routines for many FFT sizes

	; The 13 levels variants (4608, 6144, 7680)

	r4delay_pass1sc256 1152K, 4608, 4
	r4delay_pass1sc256 1536K, 6144, 4
	r4delay_pass1sc256 1920K, 7680, 4

	r4delay_pass1sc320 1440K, 4608, 4
	r4delay_pass1sc320 1920K, 6144, 4
	r4delay_pass1sc320 2400K, 7680, 4

	r4delay_pass1sc384 1728K, 4608, 4
	r4delay_pass1sc384 2304K, 6144, 4
	r4delay_pass1sc384 2880K, 7680, 4

toonear	r4delay_pass1sc448 2016K, 4608, 4
	r4delay_pass1sc448 2688K, 6144, 4
	r4delay_pass1sc448 3360K, 7680, 4

	r4delay_pass1sc512 2304K, 4608, 4
	r4delay_pass1sc512 3M, 6144, 4
	r4delay_pass1sc512 3840K, 7680, 4

	r4delay_pass1sc640 2880K, 4608, 4
	r4delay_pass1sc640 3840K, 6144, 4
	r4delay_pass1sc640 4800K, 7680, 4

	r4delay_pass1sc768 3456K, 4608, 4
	r4delay_pass1sc768 4608K, 6144, 4
	r4delay_pass1sc768 5760K, 7680, 4

toonear	r4delay_pass1sc896 4032K, 4608, 4
	r4delay_pass1sc896 5376K, 6144, 4
	r4delay_pass1sc896 6720K, 7680, 4

	r4delay_pass1sc1024 4608K, 4608, 4
	r4delay_pass1sc1024 6M, 6144, 4
	r4delay_pass1sc1024 7680K, 7680, 4

	r4delay_pass1sc1280 5760K, 4608, 4
	r4delay_pass1sc1280 7680K, 6144, 4
	r4delay_pass1sc1280 9600K, 7680, 4

	r4delay_pass1sc1536 6912K, 4608, 4
	r4delay_pass1sc1536 9M, 6144, 4
	r4delay_pass1sc1536 11520K, 7680, 4

toonear	r4delay_pass1sc1792 8064K, 4608, 4
	r4delay_pass1sc1792 10752K, 6144, 4
	r4delay_pass1sc1792 13440K, 7680, 4

	r4delay_pass1sc2048 9M, 4608, 4
	r4delay_pass1sc2048 12M, 6144, 4
	r4delay_pass1sc2048 15M, 7680, 4

	r4delay_pass1sc2560 11520K, 4608, 4
	r4delay_pass1sc2560 15M, 6144, 4
	r4delay_pass1sc2560 19200K, 7680, 4

	r4delay_pass1sc3072 13824K, 4608, 4
	r4delay_pass1sc3072 18M, 6144, 4
	r4delay_pass1sc3072 23040K, 7680, 4

toonear	r4delay_pass1sc3584 16128K, 4608, 4
	r4delay_pass1sc3584 21M, 6144, 4
	r4delay_pass1sc3584 26880K, 7680, 4

	r4delay_pass1sc4096 18M, 4608, 4
	r4delay_pass1sc4096 18M, 4608, 2
	r4delay_pass1sc4096 24M, 6144, 4
	r4delay_pass1sc4096 30M, 7680, 4

	; The all-complex 13 levels variants (4608, 6144, 7680)

	r4delay_pass1sc256ac 1152K, 4608, 4
	r4delay_pass1sc256ac 1536K, 6144, 4
	r4delay_pass1sc256ac 1920K, 7680, 4

	r4delay_pass1sc384ac 1728K, 4608, 4
	r4delay_pass1sc384ac 2304K, 6144, 4
	r4delay_pass1sc384ac 2880K, 7680, 4

	r4delay_pass1sc512ac 2304K, 4608, 4
	r4delay_pass1sc512ac 3M, 6144, 4
	r4delay_pass1sc512ac 3840K, 7680, 4

	r4delay_pass1sc640ac 2880K, 4608, 4
	r4delay_pass1sc640ac 3840K, 6144, 4
	r4delay_pass1sc640ac 4800K, 7680, 4

	r4delay_pass1sc768ac 3456K, 4608, 4
	r4delay_pass1sc768ac 4608K, 6144, 4
	r4delay_pass1sc768ac 5760K, 7680, 4

	r4delay_pass1sc1024ac 4608K, 4608, 4
	r4delay_pass1sc1024ac 6M, 6144, 4
	r4delay_pass1sc1024ac 7680K, 7680, 4

	r4delay_pass1sc1280ac 5760K, 4608, 4
	r4delay_pass1sc1280ac 7680K, 6144, 4
	r4delay_pass1sc1280ac 9600K, 7680, 4

	r4delay_pass1sc1536ac 6912K, 4608, 4
	r4delay_pass1sc1536ac 9M, 6144, 4
	r4delay_pass1sc1536ac 11520K, 7680, 4

	r4delay_pass1sc2048ac 9M, 4608, 4
	r4delay_pass1sc2048ac 12M, 6144, 4
	r4delay_pass1sc2048ac 15M, 7680, 4

	r4delay_pass1sc2560ac 11520K, 4608, 4
	r4delay_pass1sc2560ac 15M, 6144, 4
	r4delay_pass1sc2560ac 19200K, 7680, 4

	r4delay_pass1sc3072ac 13824K, 4608, 4
	r4delay_pass1sc3072ac 18M, 6144, 4
	r4delay_pass1sc3072ac 23040K, 7680, 4

	r4delay_pass1sc4096ac 18M, 4608, 4
	r4delay_pass1sc4096ac 24M, 6144, 4
	r4delay_pass1sc4096ac 30M, 7680, 4

	r4delay_pass1sc5120ac 23040K, 4608, 4
	r4delay_pass1sc5120ac 30M, 6144, 4
toobig	r4delay_pass1sc5120ac 38400K, 7680, 4

_TEXT	ENDS
END
