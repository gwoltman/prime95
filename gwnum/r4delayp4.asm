; Copyright 2009-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4/8 DJB FFT with delayed sin/cos multiplies optimized for the Pentium 4
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
xarch	TEXTEQU	<P4>

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

PREFETCHING = 128

;**********************************************************
; Late model P4s do not require TLB priming.  These P4s
; generally (always?) have at least 1MB of L2 cache.
;**********************************************************

TLB_PRIMING = 0

;; Generate pass 2 routines optimized for this architecture

	xpass2gen 12

;; Routines for many FFT sizes

	r4delay_pass1sc384 1536K, 12, 4
	r4delay_pass1sc512 2M, 12, 4
	r4delay_pass1sc1024 4M, 12, 4

	r4delay_pass1sc512ac 2M, 12, 4
	r4delay_pass1sc1024ac 4M, 12, 4

;**********************************************************
; Early model P4s require TLB priming.  These P4s generally
; (always?) have at least 512KB of L2 cache or less.
;**********************************************************

TLB_PRIMING = 1

;; Generate pass 2 routines optimized for this architecture

	xpass2gen 12

;; Routines for many FFT sizes

	r4delay_pass1sc384 1536K, 12, 4
	r4delay_pass1sc384 1536K, 12, 1
	r4delay_pass1sc512 2M, 12, 4
	r4delay_pass1sc512 2M, 12, 2
	r4delay_pass1sc512 2M, 12, 1
	r4delay_pass1sc1024 4M, 12, 2
	r4delay_pass1sc1024 4M, 12, 1

	r4delay_pass1sc512ac 2M, 12, 4
	r4delay_pass1sc512ac 2M, 12, 1
	r4delay_pass1sc1024ac 4M, 12, 2

_TEXT	ENDS
END
