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

;; Generate pass 2 routines optimized for this architecture

	xpass2gen 1536
	xpass2gen 11
	xpass2gen 2560

;; Routines for many FFT sizes

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

	; The all-complex 11 levels variants (1536, 2048, 2560)

	r4delay_pass1sc256ac 384K, 1536, 4
	r4delay_pass1sc256ac 512K, 11, 4
	r4delay_pass1sc256ac 640K, 2560, 4

	r4delay_pass1sc384ac 576K, 1536, 4
	r4delay_pass1sc384ac 768K, 11, 4
	r4delay_pass1sc384ac 960K, 2560, 4

	r4delay_pass1sc512ac 768K, 1536, 4
	r4delay_pass1sc512ac 1M, 11, 4
	r4delay_pass1sc512ac 1280K, 2560, 4

	r4delay_pass1sc640ac 960K, 1536, 4
	r4delay_pass1sc640ac 1280K, 11, 4
	r4delay_pass1sc640ac 1600K, 2560, 4

	r4delay_pass1sc768ac 1152K, 1536, 4
	r4delay_pass1sc768ac 1536K, 11, 4
	r4delay_pass1sc768ac 1920K, 2560, 4

	r4delay_pass1sc1024ac 1536K, 1536, 4
	r4delay_pass1sc1024ac 2M, 11, 4
	r4delay_pass1sc1024ac 2560K, 2560, 4

	r4delay_pass1sc1280ac 1920K, 1536, 4
	r4delay_pass1sc1280ac 2560K, 11, 4
	r4delay_pass1sc1280ac 3200K, 2560, 4

_TEXT	ENDS
END
