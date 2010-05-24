; Copyright 2009-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4 FFTs optimized for the first generation AMD64
; architecture.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

xfft_type TEXTEQU <r4>
xarch	TEXTEQU	<K8>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE r4.mac
INCLUDE r4pass1sc.mac
INCLUDE r4pass2.mac

EXTRN	pass1_aux_entry_point_return:PROC
EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC

_TEXT SEGMENT

;; Define names and features for this architecture and FFT type

PREFETCHING = 64
TLB_PRIMING = 0

;; Generate pass 2 routines optmized for this architecture

;	xpass2gen 8

;; Routines for many FFT sizes

;	r4_pass1sc32 8K, 8, 4

_TEXT	ENDS
END
