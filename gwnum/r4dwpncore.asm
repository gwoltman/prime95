; Copyright 2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4/8 DJB FFT with delayed sin/cos multiplies and partial normalization
; optimized for the Core 2 / Core i7 architecture.  This is identical to the r4delay FFTs
; except that the some of the normalization multipliers have been combined with the FFT
; sin/cos data.  Saves lots of multiplies during normalization.
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

;; Some small FFTs fit completely in a large L2 cache and thus do
;; not need to be prefetched.

PREFETCHING = 0

	xpass2gen 8
	xpass2gen 768
	xpass2gen 10
	xpass2gen 1280

	r4dwpn_pass1sc128 32K, 8, 4
	r4dwpn_pass1sc256 64K, 8, 4
	r4dwpn_pass1sc320 80K, 8, 4
	r4dwpn_pass1sc384 96K, 8, 4
	r4dwpn_pass1sc448 112K, 8, 4
	r4dwpn_pass1sc512 128K, 8, 4
	r4dwpn_pass1sc640 160K, 8, 4
	r4dwpn_pass1sc768 192K, 8, 4
	r4dwpn_pass1sc896 224K, 8, 4
	r4dwpn_pass1sc1024 256K, 8, 4
	r4dwpn_pass1sc128 96K, 768, 4
	r4dwpn_pass1sc256 192K, 768, 4
	r4dwpn_pass1sc320 240K, 768, 4
	r4dwpn_pass1sc128 128K, 10, 4
	r4dwpn_pass1sc256 256K, 10, 4
	r4dwpn_pass1sc128 160K, 1280, 4

	r4dwpn_pass1sc128ac 32K, 8, 4
	r4dwpn_pass1sc256ac 64K, 8, 4
	r4dwpn_pass1sc384ac 96K, 8, 4
	r4dwpn_pass1sc512ac 128K, 8, 4
	r4dwpn_pass1sc640ac 160K, 8, 4
	r4dwpn_pass1sc768ac 192K, 8, 4
	r4dwpn_pass1sc1024ac 256K, 8, 4
	r4dwpn_pass1sc128ac 96K, 768, 4
	r4dwpn_pass1sc256ac 192K, 768, 4
	r4dwpn_pass1sc128ac 128K, 10, 4
	r4dwpn_pass1sc256ac 256K, 10, 4
	r4dwpn_pass1sc128ac 160K, 1280, 4

;; Define names and features for this architecture and FFT type

PREFETCHING = 64
TLB_PRIMING = 0

;; Generate pass 2 routines optimized for this architecture

	xpass2gen 8

	xpass2gen 768
	xpass2gen 10
	xpass2gen 1280

;; Routines for many FFT sizes

	; The 8 levels variants

	r4dwpn_pass1sc128 32K, 8, 4
	r4dwpn_pass1sc256 64K, 8, 4
	r4dwpn_pass1sc320 80K, 8, 4
	r4dwpn_pass1sc384 96K, 8, 4
	r4dwpn_pass1sc448 112K, 8, 4
	r4dwpn_pass1sc512 128K, 8, 4
	r4dwpn_pass1sc640 160K, 8, 4
	r4dwpn_pass1sc768 192K, 8, 4
	r4dwpn_pass1sc896 224K, 8, 4
	r4dwpn_pass1sc1024 256K, 8, 4

	; The 10 levels variants (768, 1024, 1280)

	r4dwpn_pass1sc128 96K, 768, 4
	r4dwpn_pass1sc128 128K, 10, 4
	r4dwpn_pass1sc128 160K, 1280, 4

	r4dwpn_pass1sc256 192K, 768, 4
	r4dwpn_pass1sc256 256K, 10, 4
	r4dwpn_pass1sc256 320K, 1280, 4

	r4dwpn_pass1sc320 240K, 768, 4
	r4dwpn_pass1sc320 320K, 10, 4
	r4dwpn_pass1sc320 400K, 1280, 4

	r4dwpn_pass1sc384 288K, 768, 4
	r4dwpn_pass1sc384 384K, 10, 4
	r4dwpn_pass1sc384 480K, 1280, 4

	r4dwpn_pass1sc448 336K, 768, 4
	r4dwpn_pass1sc448 448K, 10, 4
	r4dwpn_pass1sc448 560K, 1280, 4

	r4dwpn_pass1sc512 384K, 768, 4
	r4dwpn_pass1sc512 512K, 10, 4
	r4dwpn_pass1sc512 640K, 1280, 4

	r4dwpn_pass1sc640 480K, 768, 4
	r4dwpn_pass1sc640 640K, 10, 4
	r4dwpn_pass1sc640 800K, 1280, 4

	r4dwpn_pass1sc768 576K, 768, 4
	r4dwpn_pass1sc768 768K, 10, 4
	r4dwpn_pass1sc768 960K, 1280, 4

	r4dwpn_pass1sc896 672K, 768, 4
	r4dwpn_pass1sc896 896K, 10, 4
	r4dwpn_pass1sc896 1120K, 1280, 4

	r4dwpn_pass1sc1024 768K, 768, 4
	r4dwpn_pass1sc1024 1M, 10, 4
	r4dwpn_pass1sc1024 1280K, 1280, 4

	r4dwpn_pass1sc1280 960K, 768, 4
	r4dwpn_pass1sc1280 1280K, 10, 4
	r4dwpn_pass1sc1280 1600K, 1280, 4

	r4dwpn_pass1sc1536 1536K, 10, 4

	r4dwpn_pass1sc1792 1792K, 10, 4

	r4dwpn_pass1sc2048 2M, 10, 4

	r4dwpn_pass1sc2560 2560K, 10, 4

	r4dwpn_pass1sc3072 3M, 10, 4

	r4dwpn_pass1sc3584 3584K, 10, 4

	r4dwpn_pass1sc4096 4M, 10, 4
	r4dwpn_pass1sc4096 4M, 10, 2

	; The all-complex 8 levels variants

	r4dwpn_pass1sc128ac 32K, 8, 4
	r4dwpn_pass1sc256ac 64K, 8, 4
	r4dwpn_pass1sc384ac 96K, 8, 4
	r4dwpn_pass1sc512ac 128K, 8, 4
	r4dwpn_pass1sc640ac 160K, 8, 4
	r4dwpn_pass1sc768ac 192K, 8, 4
	r4dwpn_pass1sc1024ac 256K, 8, 4
	r4dwpn_pass1sc1280ac 320K, 8, 4

	; The all-complex 10 levels variants (768, 1024, 1280)

	r4dwpn_pass1sc128ac 96K, 768, 4
	r4dwpn_pass1sc128ac 128K, 10, 4
	r4dwpn_pass1sc128ac 160K, 1280, 4

	r4dwpn_pass1sc256ac 192K, 768, 4
	r4dwpn_pass1sc256ac 256K, 10, 4
	r4dwpn_pass1sc256ac 320K, 1280, 4

	r4dwpn_pass1sc384ac 288K, 768, 4
	r4dwpn_pass1sc384ac 384K, 10, 4
	r4dwpn_pass1sc384ac 480K, 1280, 4

	r4dwpn_pass1sc512ac 384K, 768, 4
	r4dwpn_pass1sc512ac 512K, 10, 4
	r4dwpn_pass1sc512ac 640K, 1280, 4

	r4dwpn_pass1sc640ac 480K, 768, 4
	r4dwpn_pass1sc640ac 640K, 10, 4
	r4dwpn_pass1sc640ac 800K, 1280, 4

	r4dwpn_pass1sc768ac 576K, 768, 4
	r4dwpn_pass1sc768ac 768K, 10, 4
	r4dwpn_pass1sc768ac 960K, 1280, 4

	r4dwpn_pass1sc1024ac 768K, 768, 4
	r4dwpn_pass1sc1024ac 1M, 10, 4
	r4dwpn_pass1sc1024ac 1280K, 1280, 4

	r4dwpn_pass1sc1280ac 960K, 768, 4
	r4dwpn_pass1sc1280ac 1280K, 10, 4
	r4dwpn_pass1sc1280ac 1600K, 1280, 4

	r4dwpn_pass1sc1536ac 1536K, 10, 4

	r4dwpn_pass1sc2048ac 2M, 10, 4

	r4dwpn_pass1sc2560ac 2560K, 10, 4

	r4dwpn_pass1sc3072ac 3M, 10, 4

	r4dwpn_pass1sc4096ac 4M, 10, 4

	r4dwpn_pass1sc5120ac 5M, 10, 4

_TEXT	ENDS
END
