; Copyright 2009-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4/8 DJB FFT with delayed sin/cos multiplies optimized for the Core 2 / Core i7
; architecture.  This is identical to the radix-4/8 DJB split premultiplier FFT
; except that the first levels use common sin/cos data and the last levels use
; extra complex multiplies to repair the introduced error.  Saves lots of sin/cos memory.
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

	xpass2gen 8

	xpass2gen 768
	xpass2gen 10
	xpass2gen 1280

	xpass2gen 1536
	xpass2gen 11
	xpass2gen 2560

;; Routines for many FFT sizes

	; The 8 levels variants

	r4delay_pass1sc256 64K, 8, 4

	; The 10 levels variants (768, 1024, 1280)

	r4delay_pass1sc256 192K, 768, 4
	r4delay_pass1sc256 256K, 10, 4
	r4delay_pass1sc256 320K, 1280, 4

	r4delay_pass1sc320 240K, 768, 4
	r4delay_pass1sc320 320K, 10, 4
	r4delay_pass1sc320 400K, 1280, 4

	r4delay_pass1sc384 288K, 768, 4
	r4delay_pass1sc384 384K, 10, 4
	r4delay_pass1sc384 480K, 1280, 4

	r4delay_pass1sc448 336K, 768, 4
	r4delay_pass1sc448 448K, 10, 4
	r4delay_pass1sc448 560K, 1280, 4

	r4delay_pass1sc512 384K, 768, 4
	r4delay_pass1sc512 512K, 10, 4
	r4delay_pass1sc512 640K, 1280, 4

	r4delay_pass1sc640 480K, 768, 4
	r4delay_pass1sc640 640K, 10, 4
	r4delay_pass1sc640 800K, 1280, 4

	r4delay_pass1sc768 576K, 768, 4
	r4delay_pass1sc768 768K, 10, 4
	r4delay_pass1sc768 960K, 1280, 4

	r4delay_pass1sc896 672K, 768, 4
	r4delay_pass1sc896 896K, 10, 4
	r4delay_pass1sc896 1120K, 1280, 4

	r4delay_pass1sc1024 768K, 768, 4
	r4delay_pass1sc1024 1M, 10, 4
	r4delay_pass1sc1024 1280K, 1280, 4

	r4delay_pass1sc1280 960K, 768, 4
	r4delay_pass1sc1280 1280K, 10, 4
	r4delay_pass1sc1280 1600K, 1280, 4

	r4delay_pass1sc1536 1536K, 10, 4

	r4delay_pass1sc1792 1792K, 10, 4

	r4delay_pass1sc2048 2M, 10, 4

	r4delay_pass1sc2560 2560K, 10, 4

	r4delay_pass1sc3072 3M, 10, 4

	r4delay_pass1sc3584 3584K, 10, 4

	r4delay_pass1sc4096 4M, 10, 4
	r4delay_pass1sc4096 4M, 10, 2

	; The 11 levels variants (1536, 2048, 2560)

	r4delay_pass1sc256 384K, 1536, 4
	r4delay_pass1sc256 512K, 11, 4
	r4delay_pass1sc256 640K, 2560, 4

	r4delay_pass1sc320 480K, 1536, 4
	r4delay_pass1sc320 640K, 11, 4
	r4delay_pass1sc320 800K, 2560, 4

	r4delay_pass1sc384 576K, 1536, 4
	r4delay_pass1sc384 768K, 11, 4
	r4delay_pass1sc384 960K, 2560, 4

	r4delay_pass1sc448 672K, 1536, 4
	r4delay_pass1sc448 896K, 11, 4
	r4delay_pass1sc448 1120K, 2560, 4

	r4delay_pass1sc512 768K, 1536, 4
	r4delay_pass1sc512 1M, 11, 4
	r4delay_pass1sc512 1280K, 2560, 4

	r4delay_pass1sc640 960K, 1536, 4
	r4delay_pass1sc640 1280K, 11, 4
	r4delay_pass1sc640 1600K, 2560, 4

	r4delay_pass1sc768 1152K, 1536, 4
	r4delay_pass1sc768 1536K, 11, 4
	r4delay_pass1sc768 1920K, 2560, 4

	r4delay_pass1sc896 1344K, 1536, 4
	r4delay_pass1sc896 1792K, 11, 4
	r4delay_pass1sc896 2240K, 2560, 4

	r4delay_pass1sc1024 1536K, 1536, 4
	r4delay_pass1sc1024 2M, 11, 4
	r4delay_pass1sc1024 2560K, 2560, 4

	r4delay_pass1sc1280 1920K, 1536, 4
	r4delay_pass1sc1280 2560K, 11, 4
	r4delay_pass1sc1280 3200K, 2560, 4

	r4delay_pass1sc1536 3M, 11, 4

	r4delay_pass1sc1792 3584K, 11, 4

	r4delay_pass1sc2048 4M, 11, 4

	r4delay_pass1sc2560 5M, 11, 4

	r4delay_pass1sc3072 6M, 11, 4

	r4delay_pass1sc3584 7M, 11, 4

	r4delay_pass1sc4096 8M, 11, 4
	r4delay_pass1sc4096 8M, 11, 2

	; The all-complex 8 levels variants

	r4delay_pass1sc256ac 64K, 8, 4

	; The all-complex 10 levels variants (768, 1024, 1280)

	r4delay_pass1sc256ac 192K, 768, 4
	r4delay_pass1sc256ac 256K, 10, 4
	r4delay_pass1sc256ac 320K, 1280, 4

	; The all-complex 11 levels variants (1536, 2048, 2560)

	r4delay_pass1sc256ac 384K, 1536, 4
	r4delay_pass1sc256ac 512K, 11, 4
	r4delay_pass1sc256ac 640K, 2560, 4

_TEXT	ENDS
END
