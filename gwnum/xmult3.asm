; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes of 40K and above are supported.
; This code does two passes, 11 levels on the second pass.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xfft3.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac
INCLUDE xnormal.mac
INCLUDE xpass1.mac

EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC

EXTRNP	xpass2_11_levels
EXTRNP	xpass2_11_levels_p

_TEXT SEGMENT

;; For small FFTs, use the v22 distance between two pass 2 data blocks.
;; Pass 2 does 11 FFT levels 2 sets of data (2 * 2^11 complex values =
;; 2^13 doubles = 64KB).

PREFETCHING = 1
blkdst = (65536+4096+128)

;; All the FFT routines for each FFT length

PROCP	_xmm_gw_ffts3
	EXPANDING = 2
	xfft	40K
	xfft	48K
	xfft	48Kp
	xfft	56K
	xfft	64K
	xfft	64Kp
	xfft	80K
	xfft	96K
	xfft	96Kp
	xfft	112K
	xfft	128K
	xfft	128Kp
	xfft	160K
	xfft	192K
	xfft	192Kp
	xfft	224K
	xfft	256K
	xfft	256Kp
	xfft	320K
	xfft	384K
	xfft	384Kp
	xfft	448K
	xfft	512K
	xfft	512Kp

INCLUDE xpass1sc.mac

	xfftclm	640K, 4
	xfftclm	640K, 2
	xfftclm	768K, 4
	xfftclm	768K, 2
	xfftclm	768Kp, 4
	xfftclm	768Kp, 2
	xfftclm	896K, 4
	xfftclm	896K, 2
	xfftclm	896K, 1
	xfftclm	1024K, 4
	xfftclm	1024K, 2
	xfftclm	1024K, 1
	xfftclm	1024Kp, 4
	xfftclm	1024Kp, 2
	xfftclm	1024Kp, 1
	xfftclm	1280K, 4
	xfftclm	1280K, 2
	xfftclm	1280K, 1
	xfftclm	1536K, 4
	xfftclm	1536K, 2
	xfftclm	1536K, 1
	xfftclm	1536Kp, 4
	xfftclm	1536Kp, 2
	xfftclm	1536Kp, 1
	xfftclm	1792K, 4
	xfftclm	1792K, 2
	xfftclm	1792K, 1
	xfftclm	1792K, 0
	xfftclm	2048K, 4
	xfftclm	2048K, 2
	xfftclm	2048K, 1
	xfftclm	2048K, 0
	xfftclm	2048Kp, 4
	xfftclm	2048Kp, 2
	xfftclm	2048Kp, 1
	xfftclm	2048Kp, 0
	xfftclm	2560K, 2
	xfftclm	2560K, 1
	xfftclm	2560K, 0
	xfftclm	3072K, 2
	xfftclm	3072K, 1
	xfftclm	3072K, 0
	xfftclm	3072Kp, 2
	xfftclm	3072Kp, 1
	xfftclm	3072Kp, 0
	xfftclm	3584K, 2
	xfftclm	3584K, 1
	xfftclm	3584K, 0
	xfftclm	4096K, 2
	xfftclm	4096K, 1
	xfftclm	4096K, 0
	xfftclm	4096Kp, 2
	xfftclm	4096Kp, 1
	xfftclm	4096Kp, 0

ENDPP	_xmm_gw_ffts3

_TEXT	ENDS
END
