; Copyright 2001-2009 Mersenne Research, Inc.  All rights reserved
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

xarch	TEXTEQU	<K8>

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

PREFETCHING = 64
TLB_PRIMING = 0

;; Generate pass 2 routines optmized for this architecture

	xpass2gen 8
	xpass2gen 10
	xpass2gen 11
	xpass2gen 12
	xpass2gen 13

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
	hg_pass1sclevels8pfa6 192K, 10, 2
	hg_pass1sclevels8pfa7 224K, 10, 2
	hg_pass1sclevels8pfa8 256K, 10, 2
	hg_pass1sclevels9pfa5 320K, 10, 1
	hg_pass1sclevels9pfa6 384K, 10, 1
	hg_pass1sclevels9pfa7 448K, 10, 1
	hg_pass1sclevels9pfa8 512K, 10, 1

	hg_pass1sclevels8pfa5 640K, 12, 2
	hg_pass1sclevels8pfa6 768K, 12, 2
	hg_pass1sclevels8pfa7 896K, 12, 2
	hg_pass1sclevels8pfa7 896K, 12, 1
	hg_pass1sclevels8pfa8 1024K, 12, 2
	hg_pass1sclevels8pfa8 1024K, 12, 1
	hg_pass1sclevels9pfa5 1280K, 12, 1
	hg_pass1sclevels9pfa6 1536K, 12, 1
	hg_pass1sclevels9pfa7 1792K, 12, 1
	hg_pass1sclevels9pfa8 2048K, 12, 1
	hg_pass1sclevels10pfa5 2560K, 12, 1
	hg_pass1sclevels10pfa6 3072K, 12, 1
	hg_pass1sclevels10pfa7 3584K, 12, 1
	hg_pass1sclevels10pfa8 4096K, 12, 1
	hg_pass1sclevels11pfa5 5M, 12, 4
	hg_pass1sclevels11pfa5 5M, 12, 1
	hg_pass1sclevels11pfa6 6M, 12, 4
	hg_pass1sclevels11pfa6 6M, 12, 1
	hg_pass1sclevels11pfa7 7M, 12, 4
	hg_pass1sclevels11pfa7 7M, 12, 1
	hg_pass1sclevels11pfa8 8M, 12, 4
	hg_pass1sclevels11pfa8 8M, 12, 1

	hg_pass1sclevels9pfa5 2560K, 13, 1
	hg_pass1sclevels9pfa6 3072K, 13, 1
	hg_pass1sclevels9pfa7 3584K, 13, 1
	hg_pass1sclevels9pfa8 4096K, 13, 1

	hg_pass1levels6complex3	12K, 8, 4
	hg_pass1levels6complex4 16K, 8, 4
	hg_pass1levels7complex3 24K, 8, 4
	hg_pass1levels7complex4 32K, 8, 4
	hg_pass1levels5complex3 48K, 11, 4
	hg_pass1levels5complex4 64K, 11, 4
	hg_pass1levels6complex3 96K, 11, 4
	hg_pass1levels6complex4	128K, 11, 4

	hg_pass1sclevels8complex3 192K, 10, 2
	hg_pass1sclevels8complex4 256K, 10, 1
	hg_pass1sclevels9complex3 384K, 10, 1
	hg_pass1sclevels9complex4 512K, 10, 1

	hg_pass1sclevels8complex3 768K, 12, 2
	hg_pass1sclevels8complex4 1024K, 12, 1
	hg_pass1sclevels9complex3 1536K, 12, 1
	hg_pass1sclevels9complex4 2048K, 12, 1
	hg_pass1sclevels10complex3 3072K, 12, 4
	hg_pass1sclevels10complex3 3072K, 12, 1
	hg_pass1sclevels10complex4 4096K, 12, 1
	hg_pass1sclevels11complex3 6M, 12, 4
	hg_pass1sclevels11complex3 6M, 12, 1
	hg_pass1sclevels11complex4 8M, 12, 4
	hg_pass1sclevels11complex4 8M, 12, 1

	hg_pass1sclevels9complex3 3072K, 13, 1
	hg_pass1sclevels9complex4 4096K, 13, 1

	;; A K8-class CPU is suited for these big FFTs.
	;; We should move these to the K10 file once
	;; we've done some K10-specific optimizations

	hg_pass1sclevels11pfa5 10M, 13, 4
	hg_pass1sclevels11pfa6 12M, 13, 4
	hg_pass1sclevels11pfa7 14M, 13, 4
	hg_pass1sclevels11pfa8 16M, 13, 4
	hg_pass1sclevels12pfa5 20M, 13, 2
	hg_pass1sclevels12pfa6 24M, 13, 2
	hg_pass1sclevels12pfa7 28M, 13, 2
	hg_pass1sclevels12pfa8 32M, 13, 2

	hg_pass1sclevels11complex3 12M, 13, 4
	hg_pass1sclevels11complex4 16M, 13, 4
	hg_pass1sclevels12complex3 24M, 13, 2
	hg_pass1sclevels12complex4 32M, 13, 2

_TEXT	ENDS
END
