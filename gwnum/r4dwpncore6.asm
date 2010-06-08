; Copyright 2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These FFTs were split out of r4dwpncore.asm because of MASM limitations.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

xfft_type TEXTEQU <r4dwpn>
xarch	TEXTEQU	<CORE>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE r4.mac
INCLUDE r4dwpnpass1sc.mac
INCLUDE r4dwpnpass2.mac

EXTRN	pass1_aux_entry_point_return:PROC
EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries_wpn:PROC
EXTRN	xgw_finish_mult_wpn:PROC

_TEXT SEGMENT

;; Define names and features for this architecture and FFT type

PREFETCHING = 64
TLB_PRIMING = 0

;; Generate pass 2 routines optimized for this architecture

	xpass2gen 13
	xpass2gen 10240
	xpass2gen 12800

;; Routines for many FFT sizes

	; The 13 levels variants (8192, 10240, 12800)

	r4dwpn_pass1sc256 2M, 13, 4
	r4dwpn_pass1sc256 2560K, 10240, 4
	r4dwpn_pass1sc256 3200K, 12800, 4

	r4dwpn_pass1sc320 2560K, 13, 4
	r4dwpn_pass1sc320 3200K, 10240, 4
	r4dwpn_pass1sc320 4000K, 12800, 4

	r4dwpn_pass1sc384 3M, 13, 4
	r4dwpn_pass1sc384 3840K, 10240, 4
	r4dwpn_pass1sc384 4800K, 12800, 4

	r4dwpn_pass1sc448 3584K, 13, 4
	r4dwpn_pass1sc448 4480K, 10240, 4
	r4dwpn_pass1sc448 5600K, 12800, 4

	r4dwpn_pass1sc512 4M, 13, 4
	r4dwpn_pass1sc512 4M, 13, 2
	r4dwpn_pass1sc512 4M, 13, 1
	r4dwpn_pass1sc512 5M, 10240, 4
	r4dwpn_pass1sc512 6400K, 12800, 4

	r4dwpn_pass1sc640 5M, 13, 4
	r4dwpn_pass1sc640 6400K, 10240, 4
	r4dwpn_pass1sc640 8000K, 12800, 4

	r4dwpn_pass1sc768 6M, 13, 4
	r4dwpn_pass1sc768 7680K, 10240, 4
	r4dwpn_pass1sc768 9600K, 12800, 4

	r4dwpn_pass1sc896 7M, 13, 4
	r4dwpn_pass1sc896 8960K, 10240, 4
	r4dwpn_pass1sc896 11200K, 12800, 4

	r4dwpn_pass1sc1024 8M, 13, 4
	r4dwpn_pass1sc1024 10M, 10240, 4
	r4dwpn_pass1sc1024 12800K, 12800, 4

	r4dwpn_pass1sc1280 10M, 13, 4
	r4dwpn_pass1sc1280 12800K, 10240, 4
	r4dwpn_pass1sc1280 16000K, 12800, 4

	r4dwpn_pass1sc1536 12M, 13, 4
	r4dwpn_pass1sc1536 15M, 10240, 4
	r4dwpn_pass1sc1536 19200K, 12800, 4

	r4dwpn_pass1sc1792 14M, 13, 4
	r4dwpn_pass1sc1792 17920K, 10240, 4
	r4dwpn_pass1sc1792 22400K, 12800, 4

	r4dwpn_pass1sc2048 16M, 13, 4
	r4dwpn_pass1sc2048 16M, 13, 2
	r4dwpn_pass1sc2048 16M, 13, 1
	r4dwpn_pass1sc2048 20M, 10240, 4
	r4dwpn_pass1sc2048 20M, 10240, 2
	r4dwpn_pass1sc2048 25M, 12800, 4

	r4dwpn_pass1sc2560 20M, 13, 4
	r4dwpn_pass1sc2560 25M, 10240, 4
	r4dwpn_pass1sc2560 32000K, 12800, 4

	r4dwpn_pass1sc3072 24M, 13, 4
	r4dwpn_pass1sc3072 30M, 10240, 4
toobig	r4dwpn_pass1sc3072 38400K, 12800, 4

	r4dwpn_pass1sc3584 28M, 13, 4
toobig	r4dwpn_pass1sc3584 35840K, 10240, 4
toobig	r4dwpn_pass1sc3584 44800K, 12800, 4

	r4dwpn_pass1sc4096 32M, 13, 4
	r4dwpn_pass1sc4096 32M, 13, 2
	r4dwpn_pass1sc4096 32M, 13, 1
toobig	r4dwpn_pass1sc4096 40M, 10240, 4
toobig	r4dwpn_pass1sc4096 50M, 12800, 4

	; The all-complex 13 levels variants (8192, 10240, 12800)

	r4dwpn_pass1sc256ac 2M, 13, 4
	r4dwpn_pass1sc256ac 2560K, 10240, 4
	r4dwpn_pass1sc256ac 3200K, 12800, 4

	r4dwpn_pass1sc384ac 3M, 13, 4
	r4dwpn_pass1sc384ac 3840K, 10240, 4
	r4dwpn_pass1sc384ac 4800K, 12800, 4

	r4dwpn_pass1sc512ac 4M, 13, 4
	r4dwpn_pass1sc512ac 4M, 13, 2
	r4dwpn_pass1sc512ac 4M, 13, 1
	r4dwpn_pass1sc512ac 5M, 10240, 4
	r4dwpn_pass1sc512ac 6400K, 12800, 4

	r4dwpn_pass1sc640ac 5M, 13, 4
	r4dwpn_pass1sc640ac 6400K, 10240, 4
	r4dwpn_pass1sc640ac 8000K, 12800, 4

	r4dwpn_pass1sc768ac 6M, 13, 4
	r4dwpn_pass1sc768ac 7680K, 10240, 4
	r4dwpn_pass1sc768ac 9600K, 12800, 4

	r4dwpn_pass1sc1024ac 8M, 13, 4
	r4dwpn_pass1sc1024ac 10M, 10240, 4
	r4dwpn_pass1sc1024ac 12800K, 12800, 4

	r4dwpn_pass1sc1280ac 10M, 13, 4
	r4dwpn_pass1sc1280ac 12800K, 10240, 4
	r4dwpn_pass1sc1280ac 16000K, 12800, 4

	r4dwpn_pass1sc1536ac 12M, 13, 4
	r4dwpn_pass1sc1536ac 15M, 10240, 4
	r4dwpn_pass1sc1536ac 19200K, 12800, 4

	r4dwpn_pass1sc2048ac 16M, 13, 4
	r4dwpn_pass1sc2048ac 16M, 13, 2
	r4dwpn_pass1sc2048ac 16M, 13, 1
	r4dwpn_pass1sc2048ac 20M, 10240, 4
	r4dwpn_pass1sc2048ac 20M, 10240, 2
	r4dwpn_pass1sc2048ac 25M, 12800, 4

	r4dwpn_pass1sc2560ac 20M, 13, 4
	r4dwpn_pass1sc2560ac 25M, 10240, 4
	r4dwpn_pass1sc2560ac 32000K, 12800, 4

	r4dwpn_pass1sc3072ac 24M, 13, 4
	r4dwpn_pass1sc3072ac 30M, 10240, 4
toobig	r4dwpn_pass1sc3072ac 38400K, 12800, 4

	r4dwpn_pass1sc4096ac 32M, 13, 4
	r4dwpn_pass1sc4096ac 32M, 13, 2
	r4dwpn_pass1sc4096ac 32M, 13, 1
toobig	r4dwpn_pass1sc4096ac 40M, 10240, 4
toobig	r4dwpn_pass1sc4096ac 50M, 12800, 4

toobig	r4dwpn_pass1sc5120ac 40M, 13, 4
toobig	r4dwpn_pass1sc5120ac 50M, 10240, 4
toobig	r4dwpn_pass1sc5120ac 64000K, 12800, 4

_TEXT	ENDS
END
