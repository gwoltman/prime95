; Copyright 2001-2009 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4 FFTs optimized for the Pentium 4 architecture.
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
INCLUDE hg.mac
INCLUDE hgpass1.mac
INCLUDE hgpass1sc.mac
INCLUDE hgpass2.mac

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

	hg_pass1sclevels10pfa7 3584K, 12, 4
	hg_pass1sclevels10pfa8 4096K, 12, 4
	hg_pass1sclevels11pfa5 5M, 12, 4
	hg_pass1sclevels11pfa6 6M, 12, 4
	hg_pass1sclevels11pfa6 6M, 12, 2
	hg_pass1sclevels11pfa7 7M, 12, 4
	hg_pass1sclevels11pfa7 7M, 12, 2
	hg_pass1sclevels11pfa8 8M, 12, 4
	hg_pass1sclevels11pfa8 8M, 12, 2

	hg_pass1sclevels10complex3 3072K, 12, 4
	hg_pass1sclevels10complex4 4096K, 12, 4
	hg_pass1sclevels11complex3 6M, 12, 4
	hg_pass1sclevels11complex3 6M, 12, 2
	hg_pass1sclevels11complex4 8M, 12, 4
	hg_pass1sclevels11complex4 8M, 12, 2

	;; We've arbitrarily decide P4-class CPUs aren't
	;; suitable for FFTs larger than 8M.

;**********************************************************
; Early model P4s require TLB priming.  These P4s generally
; (always?) have at least 512KB of L2 cache or less.
;**********************************************************

TLB_PRIMING = 1

;; Generate pass 2 routines optmized for this architecture

	xpass2gen 8
	xpass2gen 10
	xpass2gen 11
	xpass2gen 12

;; Routines for many FFT sizes

	hg_pass1levels6pfa5 10K, 8, 4
	hg_pass1levels6pfa6 12K, 8, 4
	hg_pass1levels6pfa7 14K, 8, 4
	hg_pass1levels6pfa8 16K, 8, 4
	hg_pass1levels7pfa5 20K, 8, 4
	hg_pass1levels7pfa6 24K, 8, 4
	hg_pass1levels7pfa7 28K, 8, 4
	hg_pass1levels7pfa8 32K, 8, 4
	hg_pass1levels5pfa5 40K, 11, 4
	hg_pass1levels5pfa6 48K, 11, 4
	hg_pass1levels5pfa7 56K, 11, 4
	hg_pass1levels5pfa8 64K, 11, 4
	hg_pass1levels6pfa5 80K, 11, 4
	hg_pass1levels6pfa6 96K, 11, 4
	hg_pass1levels6pfa7 112K, 11, 4
	hg_pass1levels6pfa8 128K, 11, 4

	hg_pass1sclevels8pfa5 160K, 10, 2
	hg_pass1sclevels8pfa6 192K, 10, 4
	hg_pass1sclevels8pfa6 192K, 10, 1
	hg_pass1sclevels8pfa7 224K, 10, 4
	hg_pass1sclevels8pfa7 224K, 10, 2
	hg_pass1sclevels8pfa8 256K, 10, 4
	hg_pass1sclevels8pfa8 256K, 10, 2
	hg_pass1sclevels9pfa5 320K, 10, 4
	hg_pass1sclevels9pfa5 320K, 10, 1
	hg_pass1sclevels9pfa6 384K, 10, 4
	hg_pass1sclevels9pfa6 384K, 10, 2
	hg_pass1sclevels9pfa7 448K, 10, 4
	hg_pass1sclevels9pfa7 448K, 10, 2
	hg_pass1sclevels9pfa8 512K, 10, 4
	hg_pass1sclevels9pfa8 512K, 10, 2

	hg_pass1sclevels9pfa5 640K, 11, 4
	hg_pass1sclevels9pfa5 640K, 11, 1
	hg_pass1sclevels9pfa6 768K, 11, 4
	hg_pass1sclevels9pfa6 768K, 11, 2
	hg_pass1sclevels9pfa7 896K, 11, 4
	hg_pass1sclevels9pfa7 896K, 11, 2
	hg_pass1sclevels9pfa8 1024K, 11, 4
	hg_pass1sclevels9pfa8 1024K, 11, 2
	hg_pass1sclevels10pfa5 1280K, 11, 2
	hg_pass1sclevels10pfa6 1536K, 11, 1
	hg_pass1sclevels10pfa7 1792K, 11, 1

	hg_pass1sclevels8pfa5 640K, 12, 2
	hg_pass1sclevels8pfa6 768K, 12, 2
	hg_pass1sclevels8pfa7 896K, 12, 2
	hg_pass1sclevels8pfa8 1024K, 12, 1
	hg_pass1sclevels9pfa5 1280K, 12, 4
	hg_pass1sclevels9pfa6 1536K, 12, 4
	hg_pass1sclevels9pfa6 1536K, 12, 1
	hg_pass1sclevels9pfa7 1792K, 12, 4
	hg_pass1sclevels9pfa7 1792K, 12, 1
	hg_pass1sclevels9pfa8 2048K, 12, 4
	hg_pass1sclevels9pfa8 2048K, 12, 2
	hg_pass1sclevels9pfa8 2048K, 12, 1
	hg_pass1sclevels10pfa5 2560K, 12, 4
	hg_pass1sclevels10pfa5 2560K, 12, 1
	hg_pass1sclevels10pfa6 3072K, 12, 4
	hg_pass1sclevels10pfa6 3072K, 12, 1
	hg_pass1sclevels10pfa7 3584K, 12, 2
	hg_pass1sclevels10pfa7 3584K, 12, 1
	hg_pass1sclevels10pfa8 4096K, 12, 2
	hg_pass1sclevels10pfa8 4096K, 12, 1
	hg_pass1sclevels11pfa5 5M, 12, 1
	hg_pass1sclevels11pfa6 6M, 12, 1
	hg_pass1sclevels11pfa7 7M, 12, 1
	hg_pass1sclevels11pfa8 8M, 12, 1

	hg_pass1levels6complex3	12K, 8, 4
	hg_pass1levels6complex4 16K, 8, 4
	hg_pass1levels7complex3 24K, 8, 4
	hg_pass1levels7complex4 32K, 8, 4
	hg_pass1levels5complex3 48K, 11, 4
	hg_pass1levels5complex4 64K, 11, 4
	hg_pass1levels6complex3 96K, 11, 4
	hg_pass1levels6complex4	128K, 11, 4

	hg_pass1sclevels8complex3 192K, 10, 4
	hg_pass1sclevels8complex3 192K, 10, 1
	hg_pass1sclevels8complex4 256K, 10, 4
	hg_pass1sclevels8complex4 256K, 10, 2
	hg_pass1sclevels9complex3 384K, 10, 4
	hg_pass1sclevels9complex3 384K, 10, 2
	hg_pass1sclevels9complex4 512K, 10, 4
	hg_pass1sclevels9complex4 512K, 10, 2

	hg_pass1sclevels9complex3 768K, 11, 4
	hg_pass1sclevels9complex3 768K, 11, 2
	hg_pass1sclevels9complex4 1024K, 11, 4
	hg_pass1sclevels9complex4 1024K, 11, 2
	hg_pass1sclevels10complex3 1536K, 11, 1

	hg_pass1sclevels8complex3 768K, 12, 2
	hg_pass1sclevels8complex4 1024K, 12, 4
	hg_pass1sclevels8complex4 1024K, 12, 1
	hg_pass1sclevels9complex3 1536K, 12, 4
	hg_pass1sclevels9complex3 1536K, 12, 1
	hg_pass1sclevels9complex4 2048K, 12, 4
	hg_pass1sclevels9complex4 2048K, 12, 1
	hg_pass1sclevels10complex3 3072K, 12, 2
	hg_pass1sclevels10complex3 3072K, 12, 1
	hg_pass1sclevels10complex4 4096K, 12, 2
	hg_pass1sclevels10complex4 4096K, 12, 1
	hg_pass1sclevels11complex3 6M, 12, 1
	hg_pass1sclevels11complex4 8M, 12, 1

;; Some P4's with really small L2 caches perform better
;; on larger FFTs if they do not prefetch

PREFETCHING = 0

	xpass2gen 10
	xpass2gen 11

	hg_pass1sclevels10pfa6 768K, 10, 1
	hg_pass1sclevels9pfa7 896K, 11, 1
	hg_pass1sclevels9pfa8 1024K, 11, 1
	hg_pass1sclevels10pfa5 1280K, 11, 1
	hg_pass1sclevels10pfa6 1536K, 11, 1
	hg_pass1sclevels10pfa7 1792K, 11, 1
	hg_pass1sclevels10pfa8 2048K, 11, 1

	hg_pass1sclevels9complex3 768K, 11, 1
	hg_pass1sclevels9complex4 1024K, 11, 1
	hg_pass1sclevels10complex3 1536K, 11, 1
	hg_pass1sclevels10complex4 2048K, 11, 1

_TEXT	ENDS
END
