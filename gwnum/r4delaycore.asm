; Copyright 2009 Mersenne Research, Inc.  All rights reserved
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

xarch	TEXTEQU	<CORE>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE r4delay.mac
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

	xpass2gen 3072
	xpass2gen 12

;; Routines for many FFT sizes

	r4delay_pass1sclevels10complex4 3072K, 3072, 4
	r4delay_pass1sclevels10complex4 4096K, 12, 4

_TEXT	ENDS
END
