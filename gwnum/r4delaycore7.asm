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

	xpass2gen 14
	xpass2gen 20480
	xpass2gen 25600

;; Routines for many FFT sizes

	; The 14 levels variants (9216, 12288, 15360, 16384, 20480, 25600)

	r4delay_pass1sc256 4M, 14, 4
	r4delay_pass1sc256 5M, 20480, 4
	r4delay_pass1sc256 6400K, 25600, 4

	r4delay_pass1sc320 5M, 14, 4
	r4delay_pass1sc320 6400K, 20480, 4
	r4delay_pass1sc320 8000K, 25600, 4

	r4delay_pass1sc384 6M, 14, 4
	r4delay_pass1sc384 7680K, 20480, 4
	r4delay_pass1sc384 9600K, 25600, 4

	r4delay_pass1sc448 7M, 14, 4
	r4delay_pass1sc448 8960K, 20480, 4
	r4delay_pass1sc448 11200K, 25600, 4

	r4delay_pass1sc512 8M, 14, 4
	r4delay_pass1sc512 10M, 20480, 4
	r4delay_pass1sc512 12800K, 25600, 4

	r4delay_pass1sc640 10M, 14, 4
	r4delay_pass1sc640 12800K, 20480, 4
	r4delay_pass1sc640 16000K, 25600, 4

	r4delay_pass1sc768 12M, 14, 4
	r4delay_pass1sc768 15M, 20480, 4
	r4delay_pass1sc768 19200K, 25600, 4

	r4delay_pass1sc896 14M, 14, 4
	r4delay_pass1sc896 17920K, 20480, 4
	r4delay_pass1sc896 22400K, 25600, 4

	r4delay_pass1sc1024 16M, 14, 4
	r4delay_pass1sc1024 20M, 20480, 4
	r4delay_pass1sc1024 25M, 25600, 4

	r4delay_pass1sc1280 20M, 14, 4
	r4delay_pass1sc1280 25M, 20480, 4
	r4delay_pass1sc1280 32000K, 25600, 4

	r4delay_pass1sc1536 24M, 14, 4
	r4delay_pass1sc1536 30M, 20480, 4
toobig	r4delay_pass1sc1536 38400K, 25600, 4

	r4delay_pass1sc1792 28M, 14, 4
toobig	r4delay_pass1sc1792 35840K, 20480, 4
toobig	r4delay_pass1sc1792 44800K, 25600, 4

	r4delay_pass1sc2048 32M, 14, 4
toobig	r4delay_pass1sc2048 40M, 20480, 4
toobig	r4delay_pass1sc2048 50M, 25600, 4

toobig	r4delay_pass1sc2560 40M, 14, 4
toobig	r4delay_pass1sc2560 50M, 20480, 4
toobig	r4delay_pass1sc2560 64000K, 25600, 4

toobig	r4delay_pass1sc3072 48M, 14, 4
toobig	r4delay_pass1sc3072 60M, 20480, 4
toobig	r4delay_pass1sc3072 75M, 25600, 4

toobig	r4delay_pass1sc3584 56M, 14, 4
toobig	r4delay_pass1sc3584 70M, 20480, 4
toobig	r4delay_pass1sc3584 89600K, 25600, 4

toobig	r4delay_pass1sc4096 64M, 14, 4
toobig	r4delay_pass1sc4096 80M, 20480, 4
toobig	r4delay_pass1sc4096 100M, 25600, 4

	; The all-complex 14 levels variants (9216, 12288, 15360, 16384, 20480, 25600)

	r4delay_pass1sc256ac 4M, 14, 4
	r4delay_pass1sc256ac 5M, 20480, 4
	r4delay_pass1sc256ac 6400K, 25600, 4

	r4delay_pass1sc384ac 6M, 14, 4
	r4delay_pass1sc384ac 7680K, 20480, 4
	r4delay_pass1sc384ac 9600K, 25600, 4

	r4delay_pass1sc512ac 8M, 14, 4
	r4delay_pass1sc512ac 10M, 20480, 4
	r4delay_pass1sc512ac 12800K, 25600, 4

	r4delay_pass1sc640ac 10M, 14, 4
	r4delay_pass1sc640ac 12800K, 20480, 4
	r4delay_pass1sc640ac 16000K, 25600, 4

	r4delay_pass1sc768ac 12M, 14, 4
	r4delay_pass1sc768ac 15M, 20480, 4
	r4delay_pass1sc768ac 19200K, 25600, 4

	r4delay_pass1sc1024ac 16M, 14, 4
	r4delay_pass1sc1024ac 20M, 20480, 4
	r4delay_pass1sc1024ac 25M, 25600, 4

	r4delay_pass1sc1280ac 20M, 14, 4
	r4delay_pass1sc1280ac 25M, 20480, 4
	r4delay_pass1sc1280ac 32000K, 25600, 4

	r4delay_pass1sc1536ac 24M, 14, 4
	r4delay_pass1sc1536ac 30M, 20480, 4
toobig	r4delay_pass1sc1536ac 38400K, 25600, 4

	r4delay_pass1sc2048ac 32M, 14, 4
toobig	r4delay_pass1sc2048ac 40M, 20480, 4
toobig	r4delay_pass1sc2048ac 50M, 25600, 4

toobig	r4delay_pass1sc2560ac 40M, 14, 4
toobig	r4delay_pass1sc2560ac 50M, 20480, 4
toobig	r4delay_pass1sc2560ac 64000K, 25600, 4

toobig	r4delay_pass1sc3072ac 48M, 14, 4
toobig	r4delay_pass1sc3072ac 60M, 20480, 4
toobig	r4delay_pass1sc3072ac 75M, 25600, 4

toobig	r4delay_pass1sc4096ac 64M, 14, 4
toobig	r4delay_pass1sc4096ac 80M, 20480, 4
toobig	r4delay_pass1sc4096ac 100M, 25600, 4

toobig	r4delay_pass1sc5120ac 80M, 14, 4
toobig	r4delay_pass1sc5120ac 100M, 20480, 4
toobig	r4delay_pass1sc5120ac 125M, 25600, 4

_TEXT	ENDS
END
