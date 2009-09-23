; Copyright 2009 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4 FFTs optimized for the Pentium 4
; architecture.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

xarch	TEXTEQU	<P4>

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

PREFETCHING = 128

;**********************************************************
; Late model P4s do not require TLB priming.  These P4s
; generally (always?) have at least 1MB of L2 cache.
;**********************************************************

TLB_PRIMING = 0

;; Generate pass 2 routines optmized for this architecture

	xpass2gen 12

;; Routines for many FFT sizes

	r4_pass1sclevels10complex4 4096K, 12, 4

;**********************************************************
; Early model P4s require TLB priming.  These P4s generally
; (always?) have at least 512KB of L2 cache or less.
;**********************************************************

TLB_PRIMING = 1

;; Generate pass 2 routines optmized for this architecture

	xpass2gen 12

;; Routines for many FFT sizes

	r4_pass1sclevels10complex4 4096K, 12, 2

_TEXT	ENDS
END
