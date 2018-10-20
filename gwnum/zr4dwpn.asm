; Copyright 2011-2018 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the AVX-512 radix-4/8 DJB FFT with delayed sin/cos multiplies and partial normalization.
;

	TITLE   setup

zfft_type TEXTEQU <r4dwpn>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE zarch.mac
INCLUDE zbasics.mac
INCLUDE zmult.mac
INCLUDE zonepass.mac
INCLUDE zr4.mac
INCLUDE zr4dwpnpass1sc.mac
INCLUDE zr4dwpnpass2.mac

_TEXT SEGMENT

;; Generate real shared pass 1 routines optimized for this architecture

buildfor SKX, zpass1gen 2,0,,1		; Build real one-pass FFT wrapper

build421 SKX, SKX, SKX,	zpass1gen 128,0,
build421 SKX, SKX, SKX,	zpass1gen 192,0,
build421 SKX, SKX, SKX,	zpass1gen 640,0,
build421 SKX, SKX, SKX,	zpass1gen 768,0,
build421 SKX, SKX, SKX,	zpass1gen 896,0,
build421 SKX, SKX, SKX,	zpass1gen 960,0,
build421 SKX, SKX, SKX,	zpass1gen 1024,0,
build421 SKX, SKX, SKX,	zpass1gen 1152,0,
build421 SKX, SKX, SKX,	zpass1gen 1280,0,
build421 SKX, SKX, SKX,	zpass1gen 1344,0,
build421 SKX, SKX, SKX,	zpass1gen 1536,0,
build421 SKX, SKX, SKX,	zpass1gen 1920,0,
build421 SKX, SKX, SKX,	zpass1gen 2048,0,
build421 SKX, SKX, SKX,	zpass1gen 2304,0,
build421 SKX, SKX, SKX,	zpass1gen 3072,0,
IFDEF NOT_YET
build421 SKX, SKX, SKX,	zpass1gen 256,0,
build421 SKX, SKX, SKX,	zpass1gen 320,0,
build421 SKX, SKX, SKX,	zpass1gen 384,0,
build421 SKX, SKX, SKX,	zpass1gen 448,0,
build421 SKX, SKX, SKX,	zpass1gen 512,0,
build421 SKX, SKX, SKX,	zpass1gen 1792,0,
ENDIF

;; Generate all-complex shared pass 1 routines optimized for this architecture

buildfor SKX, zpass1gen 2,1,,1		; Build complex one-pass FFT wrapper

;build421 SKX, SKX, SKX,	zpass1gen 80,1,			BUG - pass 1 length 80 not working
build421 SKX, SKX, SKX,	zpass1gen 128,1,
build421 SKX, SKX, SKX,	zpass1gen 192,1,
build421 SKX, SKX, SKX,	zpass1gen 640,1,
build421 SKX, SKX, SKX,	zpass1gen 768,1,
build421 SKX, SKX, SKX,	zpass1gen 896,1,
build421 SKX, SKX, SKX,	zpass1gen 960,1,
build421 SKX, SKX, SKX,	zpass1gen 1024,1,
build421 SKX, SKX, SKX,	zpass1gen 1152,1,
build421 SKX, SKX, SKX,	zpass1gen 1280,1,
build421 SKX, SKX, SKX,	zpass1gen 1344,1,
build421 SKX, SKX, SKX,	zpass1gen 1536,1,
build421 SKX, SKX, SKX,	zpass1gen 1920,1,
build421 SKX, SKX, SKX,	zpass1gen 2048,1,
build421 SKX, SKX, SKX,	zpass1gen 2304,1,
build421 SKX, SKX, SKX,	zpass1gen 3072,1,
IFDEF NOT_YET
build421 SKX, SKX, SKX,	zpass1gen 256,1,
build421 SKX, SKX, SKX,	zpass1gen 384,1,
build421 SKX, SKX, SKX,	zpass1gen 512,1,
ENDIF

;; Generate real and complex pass 2 routines optimized for this architecture

buildfor SKX,	zpass2gen 64
buildfor SKX,	zpass2gen 320
buildfor SKX,	zpass2gen 384
buildfor SKX,	zpass2gen 448
buildfor SKX,	zpass2gen 512
buildfor SKX,	zpass2gen 640
buildfor SKX,	zpass2gen 768
buildfor SKX,	zpass2gen 1024
buildfor SKX,	zpass2gen 1600
buildfor SKX,	zpass2gen 1920
buildfor SKX,	zpass2gen 2240
buildfor SKX,	zpass2gen 2304
buildfor SKX,	zpass2gen 2560
buildfor SKX,	zpass2gen 2688
buildfor SKX,	zpass2gen 3072
buildfor SKX,	zpass2gen 3136
buildfor SKX,	zpass2gen 3200
buildfor SKX,	zpass2gen 3584
buildfor SKX,	zpass2gen 3840
buildfor SKX,	zpass2gen 4096
buildfor SKX,	zpass2gen 4480
buildfor SKX,	zpass2gen 4608
buildfor SKX,	zpass2gen 5120
buildfor SKX,	zpass2gen 5376
buildfor SKX,	zpass2gen 6144
buildfor SKX,	zpass2gen 6400
buildfor SKX,	zpass2gen 7168
buildfor SKX,	zpass2gen 7680
buildfor SKX,	zpass2gen 8192
buildfor SKX,	zpass2gen 9216
buildfor SKX,	zpass2gen 10240
buildfor SKX,	zpass2gen 12288
buildfor SKX,	zpass2gen 12800
buildfor SKX,	zpass2gen 15360
buildfor SKX,	zpass2gen 16384
buildfor SKX,	zpass2gen 17920
buildfor SKX,	zpass2gen 18432
buildfor SKX,	zpass2gen 20480
buildfor SKX,	zpass2gen 21504
buildfor SKX,	zpass2gen 24576
buildfor SKX,	zpass2gen 25088
buildfor SKX,	zpass2gen 28672
buildfor SKX,	zpass2gen 30720
buildfor SKX,	zpass2gen 32768
IFDEF NOT_YET
buildfor SKX,	zpass2gen 48
buildfor SKX,	zpass2gen 80
buildfor SKX,	zpass2gen 192
buildfor SKX,	zpass2gen 8
buildfor SKX,	zpass2gen 1280
buildfor SKX,	zpass2gen 1536
buildfor SKX,	zpass2gen 11
buildfor SKX,	zpass2gen 25600
ENDIF

_TEXT	ENDS
END
