; Copyright 2011-2024 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the small one pass AVX-512 FFTs.
;

	TITLE   setup

zfft_type TEXTEQU <r4>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE zarch.mac
INCLUDE zbasics.mac
INCLUDE zmult.mac
INCLUDE zonepass.mac
INCLUDE zr4.mac
INCLUDE zr4dwpnpass1sc.mac

_TEXT SEGMENT

;; To reduce executable size, we call common code to do the very few real FFT words in the last FFT level

PROCF	zreal_fft_final1
	int_prolog 0,16,0
	zr8_sixteen_reals_eight_complex_fft_final rsi, 8*128, 128, 2*128, 4*128
	int_epilog 0,16,0
zreal_fft_final1 ENDP

PROCF	zreal_square1
	int_prolog 0,16,0
	zr8_sixteen_reals_eight_complex_with_square rsi, 8*128, 128, 2*128, 4*128
	int_epilog 0,16,0
zreal_square1 ENDP

PROCF	zreal_mult1
	int_prolog 0,16,0
	zr8_sixteen_reals_eight_complex_with_mult rsi, 8*128, 128, 2*128, 4*128
	int_epilog 0,16,0
zreal_mult1 ENDP

PROCF	zreal_mulf1
	int_prolog 0,16,0
	zr8_sixteen_reals_eight_complex_with_mulf rsi, 8*128, 128, 2*128, 4*128
	int_epilog 0,16,0
zreal_mulf1 ENDP

EXTRN zcomplex_square_opcode:PROC
EXTRN zcomplex_mult_opcode:PROC
EXTRN zcomplex_mulf_opcode:PROC

;; Implement the small custom-coded one pass FFTs

buildfor SKX,	zonepass 128, 0
buildfor SKX,	zonepass 256, 0
buildfor SKX,	zonepass 384, 0
buildfor SKX,	zonepass 512, 0
buildfor SKX,	zonepass 640, 0
buildfor SKX,	zonepass 768, 0
buildfor SKX,	zonepass 896, 0
buildfor SKX,	zonepass 1K, 0
buildfor SKX,	zonepass 1536, 0
buildfor SKX,	zonepass 2K, 0
buildfor SKX,	zonepass 2560, 0
buildfor SKX,	zonepass 3K, 0
buildfor SKX,	zonepass 4K, 0
buildfor SKX,	zonepass 5K, 0
buildfor SKX,	zonepass 6K, 0
buildfor SKX,	zonepass 7K, 0
buildfor SKX,	zonepass 8K, 0
IFDEF NOT_YET_IMPL
buildfor SKX,	zonepass 1280, 0
buildfor     ,	zonepass 10K, 0
buildfor     ,	zonepass 12K, 0
buildfor     ,	zonepass 16K, 0
buildfor     ,	zonepass 20K, 0
buildfor     ,	zonepass 18K, 0
buildfor     ,	zonepass 24K, 0
buildfor     ,	zonepass 32K, 0
ENDIF

buildfor SKX,	zonepass 128, 1
buildfor SKX,	zonepass 256, 1
buildfor SKX,	zonepass 384, 1
buildfor SKX,	zonepass 512, 1
buildfor SKX,	zonepass 640, 1
buildfor SKX,	zonepass 768, 1
buildfor SKX,	zonepass 896, 1
buildfor SKX,	zonepass 1K, 1
buildfor SKX,	zonepass 1536, 1
buildfor SKX,	zonepass 2K, 1
buildfor SKX,	zonepass 2560, 1
buildfor SKX,	zonepass 3K, 1
buildfor SKX,	zonepass 4K, 1
buildfor SKX,	zonepass 5K, 1
buildfor SKX,	zonepass 6K, 1
buildfor SKX,	zonepass 7K, 1
buildfor SKX,	zonepass 8K, 1
IFDEF NOT_YET_IMPL
buildfor SKX,	zonepass 1280, 1
buildfor SKX,	zonepass 10K, 1
buildfor SKX,	zonepass 12K, 1
buildfor SKX,	zonepass 16K, 1
buildfor     ,	zonepass 20K, 1
buildfor     ,	zonepass 18K, 1
buildfor     ,	zonepass 24K, 1
buildfor     ,	zonepass 32K, 1
ENDIF

_TEXT	ENDS
END
