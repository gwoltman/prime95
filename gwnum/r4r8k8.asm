; Copyright 2009 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4/8 FFTs optimized for the first generation AMD64
; architecture.  This is identical to the radix-4 split premultiplier FFT
; except that any odd FFT levels are done with a radix-8 step instead of a radix-2 step.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

xarch	TEXTEQU	<K8>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE r4r8.mac
INCLUDE r4r8pass1sc.mac
INCLUDE r4r8pass2.mac

EXTRN	pass1_aux_entry_point_return:PROC
EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC

_TEXT SEGMENT

;; Define names and features for this architecture and FFT type

PREFETCHING = 64
TLB_PRIMING = 0

;; Generate pass 2 routines optmized for this architecture

	xpass2gen 12

;; Routines for many FFT sizes

	r4r8_pass1sclevels10complex4 4096K, 12, 1

_TEXT	ENDS
END
