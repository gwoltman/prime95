; Copyright 1995-2009 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the setup, common routines, and global variables
; for the various discrete-weighted transforms
;

	TITLE   gwdata

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE xmult.mac

VERSION_NUMBER = 2601		;; Version 26.1

;; List of possible suffixes for the FFTs

SFX_486		EQU	1	;; 486/Pentium version - no suffix

SFX_PPRO	EQU	2	;; Pentium Pro / P2 optimized
SFX_PPRO_CLM2	EQU	4

SFX_P3		EQU	8	;; P3 and Athlon optimized
SFX_P3_CLM2	EQU	10h

SFX_P4		EQU	20h	;; P4 optimized (SSE2, no suffix)
SFX_P4_CLM0	EQU	40h
SFX_P4_CLM1	EQU	80h
SFX_P4_CLM2	EQU	100h
SFX_P4_CLM4	EQU	200h
SFX_P4_CLM8	EQU	400h
SFX_P4_CLM10	EQU	SFX_P4_CLM1 + SFX_P4_CLM0
SFX_P4_CLM20	EQU	SFX_P4_CLM2 + SFX_P4_CLM0
SFX_P4_CLM21	EQU	SFX_P4_CLM2 + SFX_P4_CLM1
SFX_P4_CLM210	EQU	SFX_P4_CLM2 + SFX_P4_CLM1 + SFX_P4_CLM0
SFX_P4_CLM41	EQU	SFX_P4_CLM4 + SFX_P4_CLM1
SFX_P4_CLM42	EQU	SFX_P4_CLM4 + SFX_P4_CLM2
SFX_P4_CLM420	EQU	SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM0
SFX_P4_CLM421	EQU	SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM1
SFX_P4_CLM4210	EQU	SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM1 + SFX_P4_CLM0
SFX_P4_CLM84	EQU	SFX_P4_CLM8 + SFX_P4_CLM4
SFX_P4_CLM8421	EQU	SFX_P4_CLM8 + SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM1

SFX_AMD64	EQU	800h	;; AMD64 optimized (SSE2, AMD suffix)
SFX_AMD64_CLM0	EQU	1000h
SFX_AMD64_CLM1	EQU	2000h
SFX_AMD64_CLM2	EQU	4000h
SFX_AMD64_CLM4	EQU	8000h
SFX_AMD64_CLM8	EQU	10000h
SFX_AMD64_CLM10 EQU	SFX_AMD64_CLM1 + SFX_AMD64_CLM0
SFX_AMD64_CLM20 EQU	SFX_AMD64_CLM2 + SFX_AMD64_CLM0
SFX_AMD64_CLM21 EQU	SFX_AMD64_CLM2 + SFX_AMD64_CLM1
SFX_AMD64_CLM210 EQU	SFX_AMD64_CLM2 + SFX_AMD64_CLM1 + SFX_AMD64_CLM0
SFX_AMD64_CLM41	EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM1
SFX_AMD64_CLM42	EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM2
SFX_AMD64_CLM421 EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM2 + SFX_AMD64_CLM1
SFX_AMD64_CLM4210 EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM2 + SFX_AMD64_CLM1 + SFX_AMD64_CLM0
SFX_AMD64_CLM84 EQU	SFX_AMD64_CLM8 + SFX_AMD64_CLM4
SFX_AMD64_CLM8421 EQU	SFX_AMD64_CLM8 + SFX_AMD64_CLM4 + SFX_AMD64_CLM2 + SFX_AMD64_CLM1

exfft	MACRO fft_length, flags, levels
	exfft1	fft_length, flags, 0, 1, 2, 4, 8, levels
	ENDM
exfft1	MACRO fft_length, flags, c0, c1, c2, c4, c8, levels
	IFNDEF X86_64
	IF flags AND SFX_486
	EXTRN	fft&fft_length:PROC
	ENDIF
	IF flags AND SFX_PPRO
	EXTRN	fft&fft_length&PPRO:PROC
	ENDIF
	IF flags AND SFX_PPRO_CLM2
	EXTRN	fft&fft_length&c2&PPRO:PROC
	ENDIF
	IF flags AND SFX_P3
	EXTRN	fft&fft_length&P3:PROC
	ENDIF
	IF flags AND SFX_P3_CLM2
	EXTRN	fft&fft_length&c2&P3:PROC
	ENDIF
	ENDIF
	IF flags AND SFX_P4
	EXTRN	x&fft_length:PROC
	ENDIF
	IF flags AND SFX_P4_CLM0
	EXTRN	x&fft_length&c0&levels:PROC
	ENDIF
	IF flags AND SFX_P4_CLM1
	EXTRN	x&fft_length&c1&levels:PROC
	ENDIF
	IF flags AND SFX_P4_CLM2
	EXTRN	x&fft_length&c2&levels:PROC
	ENDIF
	IF flags AND SFX_P4_CLM4
	EXTRN	x&fft_length&c4&levels:PROC
	ENDIF
	IF flags AND SFX_P4_CLM8
	EXTRN	x&fft_length&c8&levels:PROC
	ENDIF
	IF flags AND SFX_AMD64
	EXTRN	x&fft_length&AMD:PROC
	ENDIF
	IF flags AND SFX_AMD64_CLM0
	EXTRN	x&fft_length&c0&levels&AMD:PROC
	ENDIF
	IF flags AND SFX_AMD64_CLM1
	EXTRN	x&fft_length&c1&levels&AMD:PROC
	ENDIF
	IF flags AND SFX_AMD64_CLM2
	EXTRN	x&fft_length&c2&levels&AMD:PROC
	ENDIF
	IF flags AND SFX_AMD64_CLM4
	EXTRN	x&fft_length&c4&levels&AMD:PROC
	ENDIF
	IF flags AND SFX_AMD64_CLM8
	EXTRN	x&fft_length&c8&levels&AMD:PROC
	ENDIF
	ENDM

	exfft	32, SFX_PPRO + SFX_P4
	exfft	40, SFX_PPRO
	exfft	48, SFX_PPRO + SFX_P4
	exfft	56, SFX_PPRO
	exfft	64, SFX_PPRO + SFX_P4
	exfft	80, SFX_PPRO + SFX_P4
	exfft	96, SFX_PPRO + SFX_P4
	exfft	112, SFX_PPRO + SFX_P4
	exfft	128, SFX_PPRO + SFX_P4
	exfft	160, SFX_PPRO + SFX_P4
	exfft	192, SFX_PPRO + SFX_P4
	exfft	224, SFX_PPRO + SFX_P4
	exfft	256, SFX_PPRO + SFX_P4
	exfft	320, SFX_PPRO + SFX_P4
	exfft	384, SFX_PPRO + SFX_P4
	exfft	448, SFX_PPRO + SFX_P4
	exfft	512, SFX_PPRO + SFX_P4
	exfft	640, SFX_PPRO + SFX_P4
	exfft	768, SFX_PPRO + SFX_P4
	exfft	896, SFX_PPRO + SFX_P4
	exfft	1024, SFX_PPRO + SFX_P4
	exfft	1280, SFX_PPRO + SFX_P4
	exfft	1536, SFX_PPRO + SFX_P4
	exfft	1792, SFX_PPRO + SFX_P4
	exfft	2048, SFX_PPRO + SFX_P4
	exfft	2560, SFX_PPRO + SFX_P4
	exfft	3072, SFX_PPRO + SFX_P4
	exfft	3584, SFX_PPRO + SFX_P4
	exfft	4096, SFX_PPRO + SFX_P4
	exfft	5120, SFX_PPRO + SFX_P4
	exfft	6144, SFX_PPRO + SFX_P4
	exfft	7168, SFX_PPRO + SFX_P4
	exfft	8192, SFX_PPRO + SFX_P4
	exfft	10K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	12K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	14K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	16K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	20K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	24K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	28K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	32K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	40K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	48K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	56K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	64K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	80K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	96K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	112K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	128K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	160K, SFX_PPRO + SFX_P3
allfft	exfft	160K, SFX_P4 + SFX_AMD64
allfft	exfft	160K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	160K, SFX_P4_CLM2 + SFX_AMD64_CLM2, 10
	exfft	192K, SFX_PPRO + SFX_P3
allfft	exfft	192K, SFX_P4 + SFX_AMD64
allfft	exfft	192K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	192K, SFX_P4_CLM41 + SFX_AMD64_CLM2, 10
	exfft	224K, SFX_PPRO + SFX_P3
allfft	exfft	224K, SFX_P4 + SFX_AMD64
allfft	exfft	224K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	224K, SFX_P4_CLM42 + SFX_AMD64_CLM2, 10
	exfft	256K, SFX_PPRO + SFX_P3
allfft	exfft	256K, SFX_P4 + SFX_AMD64
allfft	exfft	256K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	256K, SFX_P4_CLM42 + SFX_AMD64_CLM2, 10
	exfft	320K, SFX_PPRO + SFX_P3
allfft	exfft	320K, SFX_P4 + SFX_AMD64
allfft	exfft	320K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
allfft	exfft	320K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	320K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 10
	exfft	384K, SFX_PPRO + SFX_P3
allfft	exfft	384K, SFX_P4 + SFX_AMD64
allfft	exfft	384K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
allfft	exfft	384K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	384K, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	448K, SFX_PPRO + SFX_P3
allfft	exfft	448K, SFX_P4 + SFX_AMD64
allfft	exfft	448K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
allfft	exfft	448K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	448K, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	512K, SFX_PPRO + SFX_P3
allfft	exfft	512K, SFX_P4 + SFX_AMD64
allfft	exfft	512K, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 10
allfft	exfft	512K, SFX_P4_CLM421 + SFX_AMD64_CLM4, 11
	exfft	512K, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	640K, SFX_PPRO + SFX_P3
allfft	exfft	640K, SFX_P4_CLM1, 8
allfft	exfft	640K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 10
allfft	exfft	640K, SFX_P4_CLM8421 + SFX_AMD64_CLM84, 11
allfft	exfft	640K, SFX_P4_CLM8421 + SFX_AMD64_CLM84, 11
allfft	exfft	640K9, SFX_P4_CLM1, 11
allfft	exfft	640K, SFX_P4_CLM42 + SFX_AMD64_CLM8421, 12
	exfft	640K, SFX_P4_CLM41, 11
	exfft	640K, SFX_P4_CLM2 + SFX_AMD64_CLM2, 12
	exfft	768K, SFX_PPRO + SFX_P3
allfft	exfft	768K, SFX_P4_CLM1, 8
allfft	exfft	768K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 10
allfft	exfft	768K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	768K9, SFX_P4_CLM10, 11
allfft	exfft	768K, SFX_P4_CLM42 + SFX_AMD64_CLM421, 12
	exfft	768K, SFX_P4_CLM0, 10
	exfft	768K, SFX_P4_CLM42, 11
	exfft	768K, SFX_P4_CLM2 + SFX_AMD64_CLM2, 12
	exfft	896K, SFX_PPRO + SFX_P3
allfft	exfft	896K, SFX_P4_CLM1, 8
allfft	exfft	896K, SFX_P4_CLM210 + SFX_AMD64_CLM2, 10
allfft	exfft	896K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	896K9, SFX_P4_CLM10, 11
allfft	exfft	896K, SFX_P4_CLM42 + SFX_AMD64_CLM421, 12
	exfft	896K, SFX_P4_CLM420, 11
	exfft	896K, SFX_P4_CLM2 + SFX_AMD64_CLM21, 12
	exfft	1024K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	1024K, SFX_P4_CLM1, 8
allfft	exfft	1024K, SFX_P4_CLM210 + SFX_AMD64_CLM2, 10
allfft	exfft	1024K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	1024K9, SFX_P4_CLM10, 11
allfft	exfft	1024K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
	exfft	1024K, SFX_P4_CLM420, 11
	exfft	1024K, SFX_P4_CLM1 + SFX_AMD64_CLM21, 12
	exfft	1280K, SFX_PPRO + SFX_P3
allfft	exfft	1280K, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	1280K, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 11
allfft	exfft	1280K9, SFX_P4_CLM10, 11
allfft	exfft	1280K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	1280K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	1280K, SFX_P4_CLM20, 11
	exfft	1280K, SFX_P4_CLM4 + SFX_AMD64_CLM1, 12
	exfft	1536K, SFX_PPRO + SFX_P3
allfft	exfft	1536K, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	1536K, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 11
allfft	exfft	1536K9, SFX_P4_CLM10, 11
allfft	exfft	1536K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	1536K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	1536K, SFX_P4_CLM10, 11
	exfft	1536K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	1792K, SFX_PPRO + SFX_P3
allfft	exfft	1792K, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	1792K, SFX_P4_CLM4210 + SFX_AMD64_CLM4210, 11
allfft	exfft	1792K9, SFX_P4_CLM10, 11
allfft	exfft	1792K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	1792K9, SFX_P4_CLM10, 12
allfft	exfft	1792K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	1792K, SFX_P4_CLM10, 11
	exfft	1792K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	2048K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	2048K, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	2048K, SFX_P4_CLM4210 + SFX_AMD64_CLM4210, 11
allfft	exfft	2048K9, SFX_P4_CLM10, 11
allfft	exfft	2048K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	2048K9, SFX_P4_CLM10, 12
allfft	exfft	2048K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	2048K, SFX_P4_CLM0, 11
	exfft	2048K, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	2560K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	2560K, SFX_P4_CLM0, 10
allfft	exfft	2560K, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	2560K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	2560K, SFX_P4_CLM4 + SFX_AMD64_CLM421, 13
	exfft	2560K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	2560K, SFX_AMD64_CLM1, 13
	exfft	3072K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	3072K, SFX_P4_CLM0, 10
allfft	exfft	3072K, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	3072K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	3072K, SFX_P4_CLM4 + SFX_AMD64_CLM421, 13
	exfft	3072K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	3072K, SFX_AMD64_CLM1, 13
	exfft	3584K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	3584K, SFX_P4_CLM0, 10
allfft	exfft	3584K, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	3584K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	3584K, SFX_P4_CLM4 + SFX_AMD64_CLM421, 13
	exfft	3584K, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	3584K, SFX_AMD64_CLM1, 13
	exfft	4096K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	4096K, SFX_P4_CLM0, 10
allfft	exfft	4096K, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	4096K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	4096K, SFX_P4_CLM42 + SFX_AMD64_CLM421, 13
	exfft	4096K, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	4096K, SFX_AMD64_CLM1, 13
allfft	exfft	5M, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	5M, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	5M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	5M, SFX_P4_CLM41 + SFX_AMD64_CLM41, 12
allfft	exfft	6M, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	6M, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	6M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	6M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	7M, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	7M, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	7M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	7M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	8M, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	8M, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	8M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	8M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	10M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	10M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	10M, SFX_P4_CLM41 + SFX_AMD64_CLM41, 13
allfft	exfft	12M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	12M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	12M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 13
allfft	exfft	14M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	14M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	14M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 13
allfft	exfft	16M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	16M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	16M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 13
allfft	exfft	20M, SFX_P4_CLM210 + SFX_AMD64_CLM4210, 13
	exfft	20M, SFX_P4_CLM21 + SFX_AMD64_CLM20, 13
allfft	exfft	24M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	24M, SFX_P4_CLM210 + SFX_AMD64_CLM20, 13
allfft	exfft	28M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	28M, SFX_P4_CLM210 + SFX_AMD64_CLM20, 13
allfft	exfft	32M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	32M, SFX_P4_CLM210 + SFX_AMD64_CLM10, 13

	exfft	32p, SFX_PPRO + SFX_P4
	exfft	48p, SFX_PPRO + SFX_P4
	exfft	64p, SFX_PPRO + SFX_P4
	exfft	96p, SFX_PPRO + SFX_P4
	exfft	128p, SFX_PPRO + SFX_P4
	exfft	192p, SFX_PPRO + SFX_P4
	exfft	256p, SFX_PPRO + SFX_P4
	exfft	384p, SFX_PPRO + SFX_P4
	exfft	512p, SFX_PPRO + SFX_P4
	exfft	768p, SFX_PPRO + SFX_P4
	exfft	1024p, SFX_PPRO + SFX_P4
	exfft	1536p, SFX_PPRO + SFX_P4
	exfft	2048p, SFX_PPRO + SFX_P4
	exfft	3072p, SFX_PPRO + SFX_P4
	exfft	4096p, SFX_PPRO + SFX_P4
	exfft	6144p, SFX_PPRO + SFX_P4
	exfft	8192p, SFX_PPRO + SFX_P4
	exfft	12Kp, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	16Kp, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	24Kp, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	32Kp, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	48Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	64Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	96Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	128Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	192Kp, SFX_PPRO + SFX_P3
allfft	exfft	192Kp, SFX_P4 + SFX_AMD64
allfft	exfft	192Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	192Kp, SFX_P4_CLM1 + SFX_AMD64_CLM2, 10
	exfft	256Kp, SFX_PPRO + SFX_P3
allfft	exfft	256Kp, SFX_P4 + SFX_AMD64
allfft	exfft	256Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	256Kp, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	384Kp, SFX_PPRO + SFX_P3
allfft	exfft	384Kp, SFX_P4 + SFX_AMD64
allfft	exfft	384Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
allfft	exfft	384Kp, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	384Kp, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	512Kp, SFX_PPRO + SFX_P3
allfft	exfft	512Kp, SFX_P4 + SFX_AMD64
allfft	exfft	512Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 10
allfft	exfft	512Kp, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	512Kp, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	768Kp, SFX_PPRO + SFX_P3
allfft	exfft	768Kp, SFX_P4_CLM1, 8
allfft	exfft	768Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 10
allfft	exfft	768Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	768Kp9, SFX_P4_CLM10, 11
allfft	exfft	768Kp, SFX_P4_CLM42 + SFX_AMD64_CLM421, 12
	exfft	768Kp, SFX_P4_CLM420, 11
	exfft	768Kp, SFX_P4_CLM2 + SFX_AMD64_CLM2, 12
	exfft	1024Kp, SFX_PPRO + SFX_P3
allfft	exfft	1024Kp, SFX_P4_CLM1, 8
allfft	exfft	1024Kp, SFX_P4_CLM210 + SFX_AMD64_CLM2, 10
allfft	exfft	1024Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	1024Kp9, SFX_P4_CLM10, 11
allfft	exfft	1024Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
	exfft	1024Kp, SFX_P4_CLM420, 11
	exfft	1024Kp, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	1536Kp, SFX_PPRO + SFX_P3
allfft	exfft	1536Kp, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	1536Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 11
allfft	exfft	1536Kp9, SFX_P4_CLM10, 11
allfft	exfft	1536Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	1536Kp, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	1536Kp, SFX_P4_CLM10, 11
	exfft	1536Kp, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	2048Kp, SFX_PPRO + SFX_P3
allfft	exfft	2048Kp, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	2048Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM4210, 11
allfft	exfft	2048Kp9, SFX_P4_CLM10, 11
allfft	exfft	2048Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	2048Kp9, SFX_P4_CLM10, 12
allfft	exfft	2048Kp, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	2048Kp, SFX_P4_CLM0, 11
	exfft	2048Kp, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	3072Kp, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	3072Kp, SFX_P4_CLM0, 10
allfft	exfft	3072Kp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	3072Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	3072Kp, SFX_P4_CLM4 + SFX_AMD64_CLM421, 13
	exfft	3072Kp, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	3072Kp, SFX_AMD64_CLM1, 13
	exfft	4096Kp, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	4096Kp, SFX_P4_CLM0, 10
allfft	exfft	4096Kp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	4096Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	4096Kp, SFX_P4_CLM42 + SFX_AMD64_CLM421, 13
	exfft	4096Kp, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	4096Kp, SFX_AMD64_CLM1, 13

	exfft	v26_4096Kp, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12

allfft	exfft	6Mp, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	6Mp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	6Mp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	6Mp, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	8Mp, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	8Mp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	8Mp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	8Mp, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	12Mp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	12Mp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	12Mp, SFX_P4_CLM21 + SFX_AMD64_CLM41, 13
allfft	exfft	16Mp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	16Mp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 13
	exfft	16Mp, SFX_P4_CLM21 + SFX_AMD64_CLM41, 13
allfft	exfft	24Mp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	24Mp, SFX_P4_CLM210 + SFX_AMD64_CLM20, 13
allfft	exfft	32Mp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	32Mp, SFX_P4_CLM210 + SFX_AMD64_CLM10, 13

;
; Global variables needed by FFT setup code
;

IFNDEF X86_64
_GWDATA SEGMENT PAGE PUBLIC 'DATA'
ELSE
_GWDATA SEGMENT PAGE
ENDIF

; Jmptable definitions

RPF	EQU	80000000h	; FFT requires prefetch capability
RPFW	EQU	40000000h	; FFT requires prefetchw (3DNow!) capability
V26	EQU	20000000h	; FFT is experimental

IFNDEF	X86_64
jmptable DD	755,	32,	0.0000036,	672
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft32PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			3, 1, 1, 1, 0
	DD	939,	40,	0.0000057,	948
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft40PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			4, 1, 1, 1, 0
	DD	1113,	48,	0.0000065,	1128
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			5, 1, 1, 1, 0
	DD	1303,	56,	0.0000084,	1356
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft56PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			6, 1, 1, 1, 0
	DD	1499,	64,	0.0000083,	1392
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			7, 1, 1, 1, 0
	DD	1857,	80,	0.0000121,	1848
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft80PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			9, 2, 1, 1, 0
	DD	2211,	96,	0.0000141,	2208
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			11, 2, 1, 1, 0
	DD	2585,	112,	0.0000179,	2616
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft112PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			13, 3, 1, 1, 0
	DD	2953,	128,	0.0000178,	2832
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			15, 3, 1, 1, 0
	DD	3663,	160,	0.0000296,	3840
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft160PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			19, 4, 2, 1, 0
	DD	4359,	192,	0.000035,	4608
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft192PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			23, 5, 2, 1, 0
	DD	5093,	224,	0.000045,	5424
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft224PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			27, 6, 3, 1, 0
	DD	5833,	256,	0.000045,	5712
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft256PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			31, 7, 3, 1, 0
	DD	7243,	320,	0.000062,	7680
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft320PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			39, 9, 2, 1, 0
	DD	8639,	384,	0.000075,	9216
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft384PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			47, 11, 2, 1, 0
	DD	10085,	448,	0.000093,	10800
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft448PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			55, 13, 3, 1, 0
	DD	11537,	512,	0.000097,	11472
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft512PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			63, 15, 3, 1, 0
	DD	14301,	640,	0.000140,	15552
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft640PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			79, 19, 4, 2, 0
	DD	17047,	768,	0.000167,	18672
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft768PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			95, 23, 5, 2, 0
	DD	19881,	896,	0.000210,	21840
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft896PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			111, 27, 6, 3, 0
	DD	22799,	1024,	0.000218,	22992
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1024PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			127, 31, 7, 3, 0
	DD	28295,	1280,	0.000302,	31152
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1280PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			159, 39, 9, 2, 0
	DD	33761,	1536,	0.000365,	37392
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1536PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			191, 47, 11, 2, 0
	DD	39411,	1792,	0.000456,	43680
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1792PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			223, 55, 13, 3, 0
	DD	45061,	2048,	0.000490,	46032
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft2048PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			255, 63, 15, 3, 0
	DD	55825,	2560,	0.000708,	62544
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft2560PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			319, 79, 19, 4, 2, 0
	DD	66519,	3072,	0.000851,	75072
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft3072PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			383, 95, 23, 5, 2, 0
	DD	77599,	3584,	0.00107,	87648
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft3584PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			447, 111, 27, 6, 3, 0
	DD	89047,	4096,	0.00113,	92112
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft4096PPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			511, 127, 31, 7, 3, 0
	DD	110400,	5120,	0.00152,	24848
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft5120PPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	131100,	6144,	0.00191,	28016
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft6144PPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	152800,	7168,	0.00226,	31232
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft7168PPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	175300,	8192,	0.00242,	34400
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft8192PPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	217700,	10240,	0.00333,	40880
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft10KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	258200,	12288,	0.00397,	47264
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft12KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	301400,	14336,	0.00488,	53696
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft14KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	346100,	16384,	0.00522,	59936
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft16KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	430300,	20480,	0.00692,	72800
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft20KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	511600,	24576,	0.00826,	85568
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft24KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	596100,	28672,	0.0101,		98384
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft28KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	683700,	32768,	0.0109,		111008
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft32KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	848800,	40960,	0.0151,		136832
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft40KP3
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	848800,	40960,	0.0151,		136832
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft40KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1009000, 49152,	0.0184,		162416
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48KP3
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1009000, 49152,	0.0184,		162416
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1177000, 57344,	0.0227,		188048
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft56KP3
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1177000, 57344,	0.0227,		188048
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft56KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1350000, 65536,	0.0252,		213152
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64KP3
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1350000, 65536,	0.0252,		213152
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1678000, 81920, 0.0360,		264752
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft80KP3
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1678000, 81920, 0.0360,		264752
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft80KPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1994000, 98304, 0.0445,		297536
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	1994000, 98304, 0.0445,		297536
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	2324000, 114688, 0.0548,	341072
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft112KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	2324000, 114688, 0.0548,	341072
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft112KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	2664000, 131072, 0.0604,	384416
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	2664000, 131072, 0.0604,	384416
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	3310000, 163840, 0.0830,	471680
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			5120		;; scratch area size
	DD			OFFSET fft160KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	3310000, 163840, 0.0830,	471680
	DD			2		;; Flags, min_l2_cache, clm
	DD			5120		;; scratch area size
	DD			OFFSET fft160KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	3933000, 196608, 0.0982,	558704
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			OFFSET fft192KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	3933000, 196608, 0.0982,	558704
	DD			2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			OFFSET fft192KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	4593000, 229376, 0.1193,	645776
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			7168		;; scratch area size
	DD			OFFSET fft224KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	4593000, 229376, 0.1193,	645776
	DD			2		;; Flags, min_l2_cache, clm
	DD			7168		;; scratch area size
	DD			OFFSET fft224KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	5264000, 262144, 0.1316,	732320
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			OFFSET fft256KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	5264000, 262144, 0.1316,	732320
	DD			2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			OFFSET fft256KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	6545000, 327680, 0.1726,	906800
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			OFFSET fft320KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	6545000, 327680, 0.1726,	906800
	DD			2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			OFFSET fft320KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	7772000, 393216, 0.2107,	1080848
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft384KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	7772000, 393216, 0.2107,	1080848
	DD			2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft384KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	9071000, 458752, 0.2520,	1254944
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			OFFSET fft448KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	9071000, 458752, 0.2520,	1254944
	DD			2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			OFFSET fft448KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	10380000, 524288, 0.2808,	1428128
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft512KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	10380000, 524288, 0.2808,	1428128
	DD			2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft512KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	12890000, 655360, 0.372,	1777232
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			OFFSET fft640KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	12890000, 655360, 0.372,	1777232
	DD			2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			OFFSET fft640KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	15310000, 786432, 0.453,	2125376
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft768KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	15310000, 786432, 0.453,	2125376
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft768KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	17890000, 917504, 0.536,	2473568
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			OFFSET fft896KP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	17890000, 917504, 0.536,	2473568
	DD			2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			OFFSET fft896KPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	20460000, 1048576, 0.600,	2819744
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft1024K2P3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	20460000, 1048576, 0.600,	2819744
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft1024K2PPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	25390000, 1310720, 0.776,	3474992
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			OFFSET fft1280KP3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	25390000, 1310720, 0.776,	3474992
	DD			2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			OFFSET fft1280KPPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	30190000, 1572864, 0.934,	4140560
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft1536KP3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	30190000, 1572864, 0.934,	4140560
	DD			2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft1536KPPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	35200000, 1835008, 1.113,	4806176
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			OFFSET fft1792KP3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	35200000, 1835008, 1.113,	4806176
	DD			2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			OFFSET fft1792KPPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	40300000, 2097152, 1.226,	5470880
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft2048K2P3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	40300000, 2097152, 1.226,	5470880
	DD			2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft2048K2PPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	50020000, 2621440, 1.636,	6803024
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			OFFSET fft2560K2P3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	50020000, 2621440, 1.636,	6803024
	DD			2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			OFFSET fft2560K2PPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	59510000, 3145728, 1.990,	8134208
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft3072K2P3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	59510000, 3145728, 1.990,	8134208
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft3072K2PPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	69360000, 3670016, 2.380,	9465440
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			OFFSET fft3584K2P3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	69360000, 3670016, 2.380,	9465440
	DD			2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			OFFSET fft3584K2PPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	79370000, 4194304, 2.604,	10794656
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft4096K2P3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	79370000, 4194304, 2.604,	10794656
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft4096K2PPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	0
jmptablep DD	755,	32,	0.000004,	976
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft32pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			4, 1, 1, 1, 0
	DD	1111,	48,	0.000007,	1608
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			6, 3, 1, 1, 0
	DD	1485,	64,	0.000010,	1952
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			8, 4, 1, 1, 0
	DD	2199,	96,	0.0000141,	3072
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			12, 3, 1, 1, 0
	DD	2947,	128,	0.000021,	3904
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			16, 4, 1, 1, 1, 0
	DD	4345,	192,	0.000035,	6288
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft192pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			24, 6, 3, 1, 0
	DD	5817,	256,	0.000051,	7808
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft256pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			32, 8, 4, 1, 1, 0
	DD	8607,	384,	0.000075,	12432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft384pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			48, 12, 3, 1, 0
	DD	11515,	512,	0.000106,	15616
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft512pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			64, 16, 4, 1, 1, 0
	DD	17001,	768,	0.000167,	25008
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft768pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			96, 24, 6, 3, 1, 0
	DD	22701,	1024,	0.000249,	31232
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1024pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			128, 32, 8, 4, 1, 0
	DD	33569,	1536,	0.000365,	49872
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1536pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			192, 48, 12, 3, 0
	DD	44951,	2048,	0.000582,	62464
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft2048pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			256, 64, 16, 4, 1, 0
	DD	66319,	3072,	0.000851,	99888
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft3072pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			384, 96, 24, 6, 3, 0
	DD	88747,	4096,	0.00135,	124928
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft4096pPPRO
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			512, 128, 32, 8, 4, 0
	DD	130600,	6144,	0.00191,	26512
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft6144pPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	174000,	8192,	0.00284,	32960
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft8192pPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	257700,	12288,	0.00397,	46000
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft12KpPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	344700,	16384,	0.00588,	58752
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft16KpPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	508600,	24576,	0.00826,	84688
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft24KpPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	679400,	32768,	0.01299,	110336
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft32KpPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1006000, 49152,	0.0184,		162352
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48KpP3
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1006000, 49152,	0.0184,		162352
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48KpPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1345000, 65536,	0.03283,	213504
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64KpP3
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1345000, 65536,	0.03283,	213504
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64KpPPRO
	DD			8		;; FFT levels done in pass2
	DD			0
	DD	1983000, 98304, 0.0445,		290512
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96KpP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	1983000, 98304, 0.0445,		290512
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96KpPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	2652000, 131072, 0.0719,	377600
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128KpP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	2652000, 131072, 0.0719,	377600
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128KpPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	3924000, 196608, 0.0982,	552496
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			OFFSET fft192KpP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	3924000, 196608, 0.0982,	552496
	DD			2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			OFFSET fft192KpPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	5242000, 262144, 0.155,		726528
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			OFFSET fft256KpP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	5242000, 262144, 0.155,		726528
	DD			2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			OFFSET fft256KpPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	7733000, 393216, 0.2107,	1076176
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft384KpP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	7733000, 393216, 0.2107,	1076176
	DD			2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft384KpPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	10320000, 524288, 0.322,	1424384
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft512KpP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	10320000, 524288, 0.322,	1424384
	DD			2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft512KpPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	15260000, 786432, 0.453,	2123824
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft768KpP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	15260000, 786432, 0.453,	2123824
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft768KpPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	20360000, 1048576, 0.681,	2820096
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft1024KpP3
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	20360000, 1048576, 0.681,	2820096
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft1024KpPPRO
	DD			10		;; FFT levels done in pass2
	DD			0
	DD	30070000, 1572864, 0.934,	4111312
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+5*64	;; scratch area size
	DD			OFFSET fft1536KpP3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	30070000, 1572864, 0.934,	4111312
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+5*64	;; scratch area size
	DD			OFFSET fft1536KpPPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	40110000, 2097152, 1.380,	5442560
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+7*64	;; scratch area size
	DD			OFFSET fft2048KpP3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	40110000, 2097152, 1.380,	5442560
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+7*64	;; scratch area size
	DD			OFFSET fft2048KpPPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	59360000, 3145728, 1.990,	8108080
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft3072Kp2P3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	59360000, 3145728, 1.990,	8108080
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft3072Kp2PPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	79100000, 4194304, 2.919,	10770432
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft4096Kp2P3
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	79100000, 4194304, 2.919,	10770432
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft4096Kp2PPRO
	DD			12		;; FFT levels done in pass2
	DD			0
	DD	0
ENDIF

;; Jump tables for the Pentium 4 SSE2 optimized code

CELE_D	EQU	257		;; 256K L2 cache of a Celeron D
WILLI	EQU	256		;; 256K L2 cache of a Willamette P4

xjmptable DD	743,	32,	0.00000111,	896
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x32
	DD			0		;; FFT levels done in pass2
	DD			4, 4
	DD			1, 1, 1, 1, 1, 0
	DD	1099,	48,	0.00000144,	1408
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x48
	DD			0		;; FFT levels done in pass2
	DD			6, 6
	DD			2, 1, 1, 1, 1, 0
	DD	1469,	64,	0.00000178,	1920
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x64
	DD			0		;; FFT levels done in pass2
	DD			8, 8
	DD			3, 1, 1, 1, 1, 0
	DD	1827,	80,	0.00000222,	2176
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x80
	DD			0		;; FFT levels done in pass2
	DD			10, 8*2048+2
	DD			4, 2, 1, 1, 1, 0
	DD	2179,	96,	0.00000259,	2432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x96
	DD			0		;; FFT levels done in pass2
	DD			12, 12
	DD			5, 2, 1, 1, 1, 0
	DD	2539,	112,	0.00000311,	2944
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x112
	DD			0		;; FFT levels done in pass2
	DD			14, (8*2048+4)*2048+2
	DD			6, 3, 1, 1, 1, 0
	DD	2905,	128,	0.00000319,	3328
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x128
	DD			0		;; FFT levels done in pass2
	DD			16, 16
	DD			7, 3, 1, 1, 1, 0
	DD	3613,	160,	0.00000450,	4736
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x160
	DD			0		;; FFT levels done in pass2
	DD			20, 16*2048+4
	DD			9, 9, 2, 1, 4, 0
	DD	4311,	192,	0.00000542,	5632
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x192
	DD			0		;; FFT levels done in pass2
	DD			24, 24
	DD			11, 11, 2, 1, 5, 0
	DD	5029,	224,	0.00000663,	6656
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x224
	DD			0		;; FFT levels done in pass2
	DD			28, (16*2048+8)*2048+4
	DD			13, 13, 3, 1, 6, 0
	DD	5755,	256,	0.00000691,	7296
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x256
	DD			0		;; FFT levels done in pass2
	DD			32, 32
	DD			15, 15, 3, 1, 7, 0
	DD	7149,	320,	0.00000928,	8448
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x320
	DD			0		;; FFT levels done in pass2
	DD			40, 32*2048+8
	DD			19, 9, 2, 4, 1, 0
	DD	8527,	384,	0.0000111,	9984
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x384
	DD			0		;; FFT levels done in pass2
	DD			48, 48
	DD			23, 11, 2, 5, 1, 0
	DD	9933,	448,	0.0000133,	11648
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x448
	DD			0		;; FFT levels done in pass2
	DD			56, (32*2048+16)*2048+8
	DD			27, 13, 3, 6, 1, 0
	DD	11359,	512,	0.0000143,	13056
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x512
	DD			0		;; FFT levels done in pass2
	DD			64, 64
	DD			31, 15, 3, 7, 1, 0
	DD	14119,	640,	0.0000215,	17408
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x640
	DD			0		;; FFT levels done in pass2
	DD			80, 64*2048+16
	DD			39, 19, 9, 9, 4*256+2, 0
	DD	16839,	768,	0.0000260,	20736
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x768
	DD			0		;; FFT levels done in pass2
	DD			96, 96
	DD			47, 23, 11, 11, 5*256+2, 0
	DD	19639,	896,	0.0000321,	24448
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x896
	DD			0		;; FFT levels done in pass2
	DD			112, (64*2048+32)*2048+16
	DD			55, 27, 13, 13, 6*256+3, 0
	DD	22477,	1024,	0.0000349,	26112
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x1024
	DD			0		;; FFT levels done in pass2
	DD			128, 128
	DD			63, 31, 15, 15, 7*256+3, 0
	DD	27899,	1280,	0.0000494,	33664
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x1280
	DD			0		;; FFT levels done in pass2
	DD			160, 128*2048+32
	DD			79, 39, 9, 4*256+19, 2, 0
	DD	33289,	1536,	0.0000601,	40320
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x1536
	DD			0		;; FFT levels done in pass2
	DD			192, 192
	DD			95, 47, 11, 5*256+23, 2, 0
	DD	38799,	1792,	0.0000719,	47232
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x1792
	DD			0		;; FFT levels done in pass2
	DD			224, (128*2048+64)*2048+32
	DD			111, 55, 13, 6*256+27, 3, 0
	DD	44339,	2048,	0.0000773,	52224
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x2048
	DD			0		;; FFT levels done in pass2
	DD			256, 256
	DD			127, 63, 15, 7*256+31, 3, 0
	DD	55099,	2560,	0.000111,	68096
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x2560
	DD			0		;; FFT levels done in pass2
	DD			320, 256*2048+64
	DD			159, 79, 9*256+19, 9*256+39, 4*256+2, 0
	DD	65729,	3072,	0.000131,	81792
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x3072
	DD			0		;; FFT levels done in pass2
	DD			384, 384
	DD			191, 95, 11*256+23, 11*256+47, 5*256+2, 0
	DD	76559,	3584,	0.000165,	95488
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x3584
	DD			0		;; FFT levels done in pass2
	DD			448, (256*2048+128)*2048+64
	DD			223, 111, 13*256+27, 13*256+55, 6*256+3, 0
	DD	87549,	4096,	0.000163,	104448
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x4096
	DD			0		;; FFT levels done in pass2
	DD			512, 512
	DD			255, 127, 15*256+31, 15*256+63, 7*256+3, 0
	DD	108800,	5120,	0.000215,	135296
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x5120
	DD			0		;; FFT levels done in pass2
	DD			640, 512*2048+128
	DD			319, 159, 9*256+39, 19*256+79, 4*256+2, 0
	DD	129900,	6144,	0.000276,	162432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x6144
	DD			0		;; FFT levels done in pass2
	DD			768, 768
	DD			383, 191, 11*256+47, 23*256+95, 5*256+2, 0
	DD	151300,	7168,	0.000374,	189568
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x7168
	DD			0		;; FFT levels done in pass2
	DD			896, (512*2048+256)*2048+128
	DD			447, 223, 13*256+55, 27*256+111, 6*256+3, 0
	DD	172700,	8192,	0.000398,	208896
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x8192
	DD			0		;; FFT levels done in pass2
	DD			1024, 1024
	DD			511, 255, 15*256+63, 31*256+127, 7*256+3, 0
	DD	214400,	10240,	0.000470,	59904
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x10KAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	214400,	10240,	0.000470,	59904
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x10K
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	255300,	12288,	0.000590,	69632
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x12KAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	255300,	12288,	0.000590,	69632
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x12K
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	297300,	14336,	0.000716,	79488
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x14KAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	297300,	14336,	0.000716,	79488
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x14K
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	340400,	16384,	0.000787,	89088
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x16KAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	340400,	16384,	0.000787,	89088
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x16K
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	423300,	20480,	0.00103,	107904
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x20KAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	423300,	20480,	0.00103,	107904
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x20K
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	504600,	24576,	0.00132,	127232
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x24KAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	504600,	24576,	0.00132,	127232
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x24K
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	587500,	28672,	0.00156,	146688
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x28KAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	587500,	28672,	0.00156,	146688
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x28K
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	671400,	32768,	0.00175,	165888
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x32KAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	671400,	32768,	0.00175,	165888
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x32K
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	835200,	40960,	0.00225,	189696
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x40KAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	835200,	40960,	0.00225,	189696
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x40K
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	835200,	40960,	0.00225,	156672
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x40KAMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	835200,	40960,	0.00225,	156672
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x40K
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	995500, 49152,	0.00279,	206208
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x48KAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	995500, 49152,	0.00279,	206208
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x48K
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	995500, 49152,	0.00279,	177664
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x48KAMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	995500, 49152,	0.00279,	177664
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x48K
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	1158000, 57344,	0.00327,	222976
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x56KAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	1158000, 57344,	0.00327,	222976
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x56K
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	1158000, 57344,	0.00327,	198784
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x56KAMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	1158000, 57344,	0.00327,	198784
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x56K
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	1325000, 65536,	0.00367,	239488
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x64KAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	1325000, 65536,	0.00367,	239488
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x64K
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	1325000, 65536,	0.00367,	219648
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x64KAMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	1325000, 65536,	0.00367,	219648
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x64K
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	1648000, 81920, 0.00474,	273408
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x80KAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	1648000, 81920, 0.00474,	273408
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x80K
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	1648000, 81920, 0.00474,	260992
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x80KAMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	1648000, 81920, 0.00474,	260992
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x80K
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	1966000, 98304, 0.00584,	306688
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x96KAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	1966000, 98304, 0.00584,	306688
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x96K
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	1966000, 98304, 0.00584,	302848
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x96KAMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	1966000, 98304, 0.00584,	302848
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x96K
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	2287000, 114688, 0.00693,	340096
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x112KAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	2287000, 114688, 0.00693,	340096
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x112K
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	2287000, 114688, 0.00693,	344832
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x112KAMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	2287000, 114688, 0.00693,	344832
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x112K
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	2614000, 131072, 0.00779,	373248
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x128KAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	2614000, 131072, 0.00779,	373248
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x128K
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	2614000, 131072, 0.00779,	386560
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x128KAMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	2614000, 131072, 0.00779,	386560
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x128K
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	3251000, 163840, 0.00914,	439168
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x160KAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	3251000, 163840, 0.00914,	471424
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET x160K410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	3251000, 163840, 0.00914,	471424
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			10240+4*128	;; scratch area size
	DP			OFFSET x160K210AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	3251000, 163840, 0.00914,	471424
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			5120+4*128	;; scratch area size
allfft	DP			OFFSET x160K110AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	3251000, 163840, 0.00914,	439168
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x160K
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	3251000, 163840, 0.00914,	471424
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET x160K410
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	3251000, 163840, 0.00914,	471424
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			10240+4*128	;; scratch area size
	DP			OFFSET x160K210
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	3251000, 163840, 0.00914,	471424
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			5120+4*128	;; scratch area size
allfft	DP			OFFSET x160K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	3875000, 196608, 0.0114,	505600
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x192KAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	3875000, 196608, 0.0114,	555392
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x192K410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	3875000, 196608, 0.0114,	555392
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET x192K210AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	3875000, 196608, 0.0114,	555392
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			6144+5*128	;; scratch area size
allfft	DP			OFFSET x192K110AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	3875000, 196608, 0.0114,	505600
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x192K
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	3875000, 196608, 0.0114,	555392
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			24576+5*128	;; scratch area size
	DP			OFFSET x192K410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
allfft	DD	3875000, 196608, 0.0114,	555392
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			12288+5*128	;; scratch area size
allfft	DP			OFFSET x192K210
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	3875000, 196608, 0.0114,	555392
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			6144+5*128	;; scratch area size
	DP			OFFSET x192K110
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	4512000, 229376, 0.0134,	572160
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x224KAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	4512000, 229376, 0.0134,	639616
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET x224K410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	4512000, 229376, 0.0134,	639616
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			14336+6*128	;; scratch area size
	DP			OFFSET x224K210AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	4512000, 229376, 0.0134,	639616
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			7168+6*128	;; scratch area size
allfft	DP			OFFSET x224K110AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	4512000, 229376, 0.0134,	572160
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x224K
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	4512000, 229376, 0.0134,	639616
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			28672+6*128	;; scratch area size
	DP			OFFSET x224K410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
	DD	4512000, 229376, 0.0134,	639616
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			14336+6*128	;; scratch area size
	DP			OFFSET x224K210
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	4512000, 229376, 0.0134,	639616
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			7168+6*128	;; scratch area size
allfft	DP			OFFSET x224K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	5158000, 262144, 0.0150,	638464
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x256KAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	5158000, 262144, 0.0150,	721920
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x256K410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	5158000, 262144, 0.0150,	721920
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			16384+7*128	;; scratch area size
	DP			OFFSET x256K210AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	5158000, 262144, 0.0150,	721920
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			8192+7*128	;; scratch area size
allfft	DP			OFFSET x256K110AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	5158000, 262144, 0.0150,	638464
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x256K
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	5158000, 262144, 0.0150,	721920
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			32768+7*128	;; scratch area size
	DP			OFFSET x256K410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
	DD	5158000, 262144, 0.0150,	721920
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			16384+7*128	;; scratch area size
	DP			OFFSET x256K210
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	5158000, 262144, 0.0150,	721920
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			8192+7*128	;; scratch area size
allfft	DP			OFFSET x256K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	6421000, 327680, 0.0192,	772480
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET x320K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	6421000, 327680, 0.0192,	772480
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x320KAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	6421000, 327680, 0.0192,	890624
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET x320K410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	6421000, 327680, 0.0192,	890624
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET x320K210AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	6421000, 327680, 0.0192,	890624
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET x320K110AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	6421000, 327680, 0.0192,	772480
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET x320K411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	6421000, 327680, 0.0192,	772480
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x320K
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	6421000, 327680, 0.0192,	890624
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			40960+9*128	;; scratch area size
	DP			OFFSET x320K410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
allfft	DD	6421000, 327680, 0.0192,	890624
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET x320K210
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	6421000, 327680, 0.0192,	890624
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET x320K110
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	7651000, 393216, 0.0238,	905600
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x384K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	7651000, 393216, 0.0238,	905600
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x384KAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	7651000, 393216, 0.0238,	1058432
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x384K410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	7651000, 393216, 0.0238,	1058432
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET x384K210AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	7651000, 393216, 0.0238,	1058432
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET x384K110AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	7651000, 393216, 0.0238,	905600
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x384K411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	7651000, 393216, 0.0238,	905600
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x384K
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	7651000, 393216, 0.0238,	1058432
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET x384K410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
	DD	7651000, 393216, 0.0238,	1058432
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			24576+11*128	;; scratch area size
	DP			OFFSET x384K210
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	7651000, 393216, 0.0238,	1058432
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET x384K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	8908000, 458752, 0.0283,	1038976
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET x448K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	8908000, 458752, 0.0283,	1038976
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x448KAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	8908000, 458752, 0.0283,	1226496
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET x448K410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	8908000, 458752, 0.0283,	1226496
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET x448K210AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	8908000, 458752, 0.0283,	1226496
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET x448K110AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	8908000, 458752, 0.0283,	1038976
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET x448K411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	8908000, 458752, 0.0283,	1038976
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x448K
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	8908000, 458752, 0.0283,	1226496
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			57344+13*128	;; scratch area size
	DP			OFFSET x448K410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
	DD	8908000, 458752, 0.0283,	1226496
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			28672+13*128	;; scratch area size
	DP			OFFSET x448K210
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	8908000, 458752, 0.0283,	1226496
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET x448K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x512K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x512KAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	10180000, 524288, 0.0319,	1392640
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET x512K410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	10180000, 524288, 0.0319,	1392640
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x512K210AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	10180000, 524288, 0.0319,	1392640
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x512K110AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x512K411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET x512K211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			8192+7*128	;; scratch area size
allfft	DP			OFFSET x512K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x512K
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	10180000, 524288, 0.0319,	1392640
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET x512K410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
	DD	10180000, 524288, 0.0319,	1392640
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET x512K210
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	10180000, 524288, 0.0319,	1392640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x512K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	10180000, 524288, 0.0319,	1392640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x512K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	12650000, 655360, 0.0410,	1214848
allfft	DD			RPFW+512*65536+2*8;; Flags, min_l2_cache, clm
allfft	DD			40960+4*128	;; scratch area size
allfft	DP			OFFSET x640K812AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_8+1, 0
allfft	DD	12650000, 655360, 0.0410,	1214848
allfft	DD			RPFW+512*65536+2*4;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET x640K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
	DD	12650000, 655360, 0.0410,	1214848
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			10240+4*128	;; scratch area size
	DP			OFFSET x640K212AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
allfft	DD	12650000, 655360, 0.0410,	1214848
allfft	DD			RPFW+512*65536+2*1;; Flags, min_l2_cache, clm
allfft	DD			5120+4*128	;; scratch area size
allfft	DP			OFFSET x640K112AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			RPFW+2*8	;; Flags, min_l2_cache, clm
allfft	DD			81920+9*128	;; scratch area size
allfft	DP			OFFSET x640K811AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_8+1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET x640K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	12650000, 655360, 0.0410,	1214848
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET x640K412
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
	DD	12650000, 655360, 0.0410,	1214848
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			10240+4*128	;; scratch area size
	DP			OFFSET x640K212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			512*65536+2*8	;; Flags, min_l2_cache, clm
allfft	DD			81920+9*128	;; scratch area size
allfft	DP			OFFSET x640K811
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_8+1, 0
	DD	12650000, 655360, 0.0410,	1437440
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			40960+9*128	;; scratch area size
	DP			OFFSET x640K411
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET x640K211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
	DD	12650000, 655360, 0.0410,	1437440
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET x640K111
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_1+1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			10240+9*128	;; scratch area size
allfft	DP			OFFSET x640K9111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	12650000, 655360, 0.0410,	1730816
allfft	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET x640K410
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	12650000, 655360, 0.0410,	1730816
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET x640K210
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
allfft	DD	12650000, 655360, 0.0410,	1730816
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET x640K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	12650000, 655360, 0.0410,	3123840
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x640K18
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			GAP2_8_1+1, 0
allfft	DD	15070000, 786432, 0.0507,	1413504
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x768K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
	DD	15070000, 786432, 0.0507,	1413504
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET x768K212AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
allfft	DD	15070000, 786432, 0.0507,	1413504
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			6144+5*128	;; scratch area size
allfft	DP			OFFSET x768K112AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x768K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x768K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	15070000, 786432, 0.0507,	1413504
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x768K412
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
	DD	15070000, 786432, 0.0507,	1413504
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET x768K212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	15070000, 786432, 0.0507,	1703552
	DD			CELE_D*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET x768K411
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	15070000, 786432, 0.0507,	1703552
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			24576+11*128	;; scratch area size
	DP			OFFSET x768K211
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_2+1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET x768K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET x768K9111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET x768K011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET x768K9011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	15070000, 786432, 0.0507,	2066816
allfft	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x768K410
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	15070000, 786432, 0.0507,	2066816
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x768K210
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
allfft	DD	15070000, 786432, 0.0507,	2066816
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET x768K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
	DD	15070000, 786432, 0.0507,	2066816
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET x768K010
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_0+1, 0
allfft	DD	15070000, 786432, 0.0507,	3746560
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x768K18
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			GAP2_8_1+1, 0
allfft	DD	17550000, 917504, 0.0607,	1612416
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET x896K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
	DD	17550000, 917504, 0.0607,	1612416
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			14336+6*128	;; scratch area size
	DP			OFFSET x896K212AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	17550000, 917504, 0.0607,	1612416
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			7168+6*128	;; scratch area size
	DP			OFFSET x896K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET x896K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET x896K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	17550000, 917504, 0.0607,	1612416
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET x896K412
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
	DD	17550000, 917504, 0.0607,	1612416
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			14336+6*128	;; scratch area size
	DP			OFFSET x896K212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	17550000, 917504, 0.0607,	1969920
	DD			CELE_D*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			57344+13*128	;; scratch area size
	DP			OFFSET x896K411
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	17550000, 917504, 0.0607,	1969920
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			28672+13*128	;; scratch area size
	DP			OFFSET x896K211
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_2+1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET x896K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET x896K9111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
	DD	17550000, 917504, 0.0607,	1969920
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET x896K011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET x896K9011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	17550000, 917504, 0.0607,	2402816
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET x896K210
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
allfft	DD	17550000, 917504, 0.0607,	2402816
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET x896K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	17550000, 917504, 0.0607,	2402816
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET x896K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	17550000, 917504, 0.0607,	4369280
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x896K18
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			GAP2_8_1+1, 0
allfft	DD	20050000, 1048576, 0.0676,	1809408
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x1024K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
	DD	20050000, 1048576, 0.0676,	1809408
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			16384+7*128	;; scratch area size
	DP			OFFSET x1024K212AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	20050000, 1048576, 0.0676,	1809408
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET x1024K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET x1024K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x1024K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	20050000, 1048576, 0.0676,	1809408
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x1024K412
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	20050000, 1048576, 0.0676,	1809408
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET x1024K212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	20050000, 1048576, 0.0676,	1809408
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET x1024K112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
	DD	20050000, 1048576, 0.0676,	2234368
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET x1024K411
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	20050000, 1048576, 0.0676,	2234368
	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET x1024K211
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_2+1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			705*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x1024K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x1024K9111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
	DD	20050000, 1048576, 0.0676,	2234368
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x1024K011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			600*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x1024K9011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	20050000, 1048576, 0.0676,	2734080
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET x1024K210
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
allfft	DD	20050000, 1048576, 0.0676,	2734080
allfft	DD			200*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x1024K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	20050000, 1048576, 0.0676,	2734080
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x1024K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	20050000, 1048576, 0.0676,	4975104
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x1024K18
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			GAP2_8_1+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2099584
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET x1280K413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2207488
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET x1280K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2207488
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET x1280K212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	24930000, 1310720, 0.0892,	2207488
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET x1280K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET x1280K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			RPFW+256*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET x1280K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET x1280K111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2099584
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET x1280K413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	24930000, 1310720, 0.0892,	2207488
	DD			CELE_D*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			40960+9*128	;; scratch area size
	DP			OFFSET x1280K412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2207488
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET x1280K212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2207488
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			10240+9*128	;; scratch area size
allfft	DP			OFFSET x1280K112
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET x1280K411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	24930000, 1310720, 0.0892,	2769152
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			40960+19*128	;; scratch area size
	DP			OFFSET x1280K211
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_2+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET x1280K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET x1280K9111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
	DD	24930000, 1310720, 0.0892,	2769152
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			20480+19*128	;; scratch area size
	DP			OFFSET x1280K011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET x1280K9011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	24930000, 1310720, 0.0892,	3409536
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET x1280K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	24930000, 1310720, 0.0892,	3409536
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET x1280K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	29690000, 1572864, 0.113,	2429312
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x1536K413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	29690000, 1572864, 0.113,	2604672
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x1536K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	29690000, 1572864, 0.113,	2604672
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET x1536K212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	29690000, 1572864, 0.113,	2604672
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET x1536K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x1536K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x1536K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET x1536K111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	29690000, 1572864, 0.113,	2429312
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x1536K413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	29690000, 1572864, 0.113,	2604672
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET x1536K412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	29690000, 1572864, 0.113,	2604672
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET x1536K212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	29690000, 1572864, 0.113,	2604672
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET x1536K112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x1536K411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x1536K211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
	DD	29690000, 1572864, 0.113,	3301760
	DD			WILLI*65536+2*1	;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET x1536K111
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_1+1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET x1536K9111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
	DD	29690000, 1572864, 0.113,	3301760
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET x1536K011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET x1536K9011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	29690000, 1572864, 0.113,	4081280
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x1536K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	29690000, 1572864, 0.113,	4081280
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x1536K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	34560000, 1835008, 0.135,	2759296
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET x1792K413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET x1792K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET x1792K212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	34560000, 1835008, 0.135,	3002112
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET x1792K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET x1792K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET x1792K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET x1792K111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET x1792K011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	34560000, 1835008, 0.135,	2759296
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET x1792K413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	34560000, 1835008, 0.135,	3002112
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			57344+13*128	;; scratch area size
	DP			OFFSET x1792K412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET x1792K212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	34560000, 1835008, 0.135,	3002112
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET x1792K112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET x1792K9112
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET x1792K9012
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET x1792K411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET x1792K211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
	DD	34560000, 1835008, 0.135,	3834368
	DD			WILLI*65536+2*1	;; Flags, min_l2_cache, clm
	DD			28672+27*128	;; scratch area size
	DP			OFFSET x1792K111
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_1+1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET x1792K9111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
	DD	34560000, 1835008, 0.135,	3834368
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			28672+27*128	;; scratch area size
	DP			OFFSET x1792K011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET x1792K9011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	34560000, 1835008, 0.135,	4753024
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET x1792K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	34560000, 1835008, 0.135,	4753024
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET x1792K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	39500000, 2097152, 0.155,	3087360
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x2048K413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	39500000, 2097152, 0.155,	3397632
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET x2048K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	39500000, 2097152, 0.155,	3397632
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x2048K212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	39500000, 2097152, 0.155,	3397632
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x2048K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET x2048K411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET x2048K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x2048K111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x2048K011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	39500000, 2097152, 0.155,	3087360
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x2048K413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	39500000, 2097152, 0.155,	3397632
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET x2048K412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	39500000, 2097152, 0.155,	3397632
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x2048K112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
	DD	39500000, 2097152, 0.155,	3397632
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET x2048K212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
allfft	DD	39500000, 2097152, 0.155,	3397632
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x2048K9112
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	39500000, 2097152, 0.155,	3397632
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x2048K9012
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET x2048K411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET x2048K211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x2048K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x2048K9111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
	DD	39500000, 2097152, 0.155,	4362240
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET x2048K011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x2048K9011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	39500000, 2097152, 0.155,	5416960
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x2048K110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	39500000, 2097152, 0.155,	5416960
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x2048K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	49100000, 2621440, 0.204,	3747584
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET x2560K413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	49100000, 2621440, 0.204,	3747584
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET x2560K213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	49100000, 2621440, 0.204,	3747584
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET x2560K113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	49100000, 2621440, 0.204,	4194560
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET x2560K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	49100000, 2621440, 0.204,	4194560
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET x2560K212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	49100000, 2621440, 0.204,	4194560
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			20480+19*128	;; scratch area size
	DP			OFFSET x2560K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET x2560K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET x2560K111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET x2560K011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	49100000, 2621440, 0.204,	3747584
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET x2560K413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	49100000, 2621440, 0.204,	4194560
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			81920+19*128	;; scratch area size
	DP			OFFSET x2560K412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	49100000, 2621440, 0.204,	4194560
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET x2560K212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	49100000, 2621440, 0.204,	4194560
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			20480+19*128	;; scratch area size
	DP			OFFSET x2560K112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET x2560K211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET x2560K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET x2560K011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	49100000, 2621440, 0.204,	6768768
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x2560K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	58520000, 3145728, 0.259,	4406912
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x3072K413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	58520000, 3145728, 0.259,	4406912
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET x3072K213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	58520000, 3145728, 0.259,	4406912
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET x3072K113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	58520000, 3145728, 0.259,	4989312
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x3072K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	58520000, 3145728, 0.259,	4989312
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x3072K212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	58520000, 3145728, 0.259,	4989312
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET x3072K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET x3072K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x3072K111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x3072K011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	58520000, 3145728, 0.259,	4406912
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x3072K413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	58520000, 3145728, 0.259,	4989312
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			98304+23*128	;; scratch area size
	DP			OFFSET x3072K412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	58520000, 3145728, 0.259,	4989312
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x3072K212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	58520000, 3145728, 0.259,	4989312
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET x3072K112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET x3072K211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x3072K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x3072K011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	58520000, 3145728, 0.259,	8112384
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x3072K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	68130000, 3670016, 0.323,	5066496
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET x3584K413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	68130000, 3670016, 0.323,	5066496
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET x3584K213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	68130000, 3670016, 0.323,	5066496
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET x3584K113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	68130000, 3670016, 0.323,	5784064
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET x3584K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	68130000, 3670016, 0.323,	5784064
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET x3584K212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	68130000, 3670016, 0.323,	5784064
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			28672+27*128	;; scratch area size
	DP			OFFSET x3584K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			114688+55*128	;; scratch area size
allfft	DP			OFFSET x3584K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET x3584K111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET x3584K011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	68130000, 3670016, 0.323,	5066496
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET x3584K413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	68130000, 3670016, 0.323,	5784064
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			114688+27*128	;; scratch area size
	DP			OFFSET x3584K412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	68130000, 3670016, 0.323,	5784064
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			57344+27*128	;; scratch area size
	DP			OFFSET x3584K212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	68130000, 3670016, 0.323,	5784064
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			28672+27*128	;; scratch area size
	DP			OFFSET x3584K112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			114688+55*128	;; scratch area size
allfft	DP			OFFSET x3584K211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET x3584K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET x3584K011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	68130000, 3670016, 0.323,	9456000
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x3584K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	77910000, 4194304, 0.382,	5724160
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET x4096K413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	77910000, 4194304, 0.382,	5724160
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x4096K213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	77910000, 4194304, 0.382,	5724160
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x4096K113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	77910000, 4194304, 0.382,	6574080
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET x4096K412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	77910000, 4194304, 0.382,	6574080
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET x4096K212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	77910000, 4194304, 0.382,	6574080
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET x4096K112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET x4096K211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x4096K111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x4096K011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	77910000, 4194304, 0.382,	5724160
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET x4096K413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	77910000, 4194304, 0.382,	5724160
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x4096K213
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	77910000, 4194304, 0.382,	6574080
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			131072+31*128	;; scratch area size
	DP			OFFSET x4096K412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	77910000, 4194304, 0.382,	6574080
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			65536+31*128	;; scratch area size
	DP			OFFSET x4096K212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	77910000, 4194304, 0.382,	6574080
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET x4096K112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET x4096K211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x4096K111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x4096K011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	77910000, 4194304, 0.382,	10782720
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x4096K010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET x5M413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET x5M213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET x5M113AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_1+1, 0
	DD	96830000, 5242880, 0.485,	8167040
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			163840+39*128	;; scratch area size
	DP			OFFSET x5M412AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	96830000, 5242880, 0.485,	8167040
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET x5M212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	96830000, 5242880, 0.485,	8167040
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			40960+39*128	;; scratch area size
	DP			OFFSET x5M112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	96830000, 5242880, 0.485,	10756224
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x5M111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	96830000, 5242880, 0.485,	10756224
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x5M011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET x5M413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET x5M213
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	96830000, 5242880, 0.485,	8167040
	DD			1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			163840+39*128	;; scratch area size
	DP			OFFSET x5M412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	96830000, 5242880, 0.485,	8167040
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET x5M212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	96830000, 5242880, 0.485,	8167040
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			40960+39*128	;; scratch area size
	DP			OFFSET x5M112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	96830000, 5242880, 0.485,	8167040
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET x5M012
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
allfft	DD	96830000, 5242880, 0.485,	10756224
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x5M111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	96830000, 5242880, 0.485,	10756224
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x5M011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x6M413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x6M213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET x6M113AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_1+1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET x6M412AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	115300000, 6291456, 0.668,	9756288
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET x6M212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET x6M112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	115300000, 6291456, 0.668,	12886272
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x6M111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	115300000, 6291456, 0.668,	12886272
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x6M011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x6M413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x6M213
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET x6M412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			98304+47*128	;; scratch area size
	DP			OFFSET x6M212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET x6M112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	115300000, 6291456, 0.668,	12886272
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x6M111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	115300000, 6291456, 0.668,	12886272
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x6M011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET x7M413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET x7M213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET x7M113AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_1+1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			229376+55*128	;; scratch area size
	DP			OFFSET x7M412AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	134200000, 7340032, 0.886,	11345536
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			114688+55*128	;; scratch area size
allfft	DP			OFFSET x7M212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			57344+55*128	;; scratch area size
	DP			OFFSET x7M112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	134200000, 7340032, 0.886,	15016320
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x7M111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	134200000, 7340032, 0.886,	15016320
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x7M011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET x7M413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET x7M213
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			229376+55*128	;; scratch area size
	DP			OFFSET x7M412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			114688+55*128	;; scratch area size
	DP			OFFSET x7M212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			57344+55*128	;; scratch area size
	DP			OFFSET x7M112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	134200000, 7340032, 0.886,	15016320
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x7M111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	134200000, 7340032, 0.886,	15016320
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x7M011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET x8M413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET x8M213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x8M113AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_1+1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET x8M412AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	153400000, 8388608, 1.042,	12926976
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET x8M212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET x8M112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	153400000, 8388608, 1.042,	17129472
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x8M111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	153400000, 8388608, 1.042,	17129472
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x8M011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET x8M413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET x8M213
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x8M113
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_1+1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET x8M412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			131072+63*128	;; scratch area size
	DP			OFFSET x8M212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET x8M112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	153400000, 8388608, 1.042,	17129472
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x8M111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	153400000, 8388608, 1.042,	17129472
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x8M011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
	DD	190700000, 10485760, 1.100,	13639296
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			163840+39*128	;; scratch area size
	DP			OFFSET x10M413AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
allfft	DD	190700000, 10485760, 1.100,	13639296
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET x10M213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	190700000, 10485760, 1.100,	13639296
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			40960+39*128	;; scratch area size
	DP			OFFSET x10M113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			163840+79*128	;; scratch area size
allfft	DP			OFFSET x10M212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x10M112AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x10M012AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
	DD	190700000, 10485760, 1.100,	13639296
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			163840+39*128	;; scratch area size
	DP			OFFSET x10M413
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
allfft	DD	190700000, 10485760, 1.100,	13639296
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET x10M213
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	190700000, 10485760, 1.100,	13639296
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			40960+39*128	;; scratch area size
	DP			OFFSET x10M113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			163840+79*128	;; scratch area size
allfft	DP			OFFSET x10M212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x10M112
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x10M012
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET x12M413AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
allfft	DD	227300000, 12582912, 1.400,	16277120
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET x12M213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET x12M113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			196608+95*128	;; scratch area size
allfft	DP			OFFSET x12M212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x12M112AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x12M012AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET x12M413
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			98304+47*128	;; scratch area size
	DP			OFFSET x12M213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET x12M113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			196608+95*128	;; scratch area size
allfft	DP			OFFSET x12M212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x12M112
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x12M012
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			229376+55*128	;; scratch area size
	DP			OFFSET x14M413AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
allfft	DD	264600000, 14680064, 1.700,	18914944
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			114688+55*128	;; scratch area size
allfft	DP			OFFSET x14M213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			57344+55*128	;; scratch area size
	DP			OFFSET x14M113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			229376+111*128	;; scratch area size
allfft	DP			OFFSET x14M212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x14M112AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x14M012AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			229376+55*128	;; scratch area size
	DP			OFFSET x14M413
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			114688+55*128	;; scratch area size
	DP			OFFSET x14M213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			57344+55*128	;; scratch area size
	DP			OFFSET x14M113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			229376+111*128	;; scratch area size
allfft	DP			OFFSET x14M212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x14M112
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x14M012
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET x16M413AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
allfft	DD	302600000, 16777216, 2.100,	21544960
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET x16M213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET x16M113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			262144+127*128	;; scratch area size
allfft	DP			OFFSET x16M212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x16M112AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x16M012AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET x16M413
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			131072+63*128	;; scratch area size
	DP			OFFSET x16M213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET x16M113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			262144+127*128	;; scratch area size
allfft	DP			OFFSET x16M212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x16M112
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x16M012
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_0+1, 0
allfft	DD	376100000, 20971520, 2.700,	26828928
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			327680+79*128	;; scratch area size
allfft	DP			OFFSET x20M413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	376100000, 20971520, 2.700,	26828928
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			163840+79*128	;; scratch area size
	DP			OFFSET x20M213AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
allfft	DD	376100000, 20971520, 2.700,	26828928
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x20M113AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_1+1, 0
	DD	376100000, 20971520, 2.700,	26828928
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			81920+79*128	;; scratch area size
	DP			OFFSET x20M013AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	376100000, 20971520, 2.700,	26828928
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			163840+79*128	;; scratch area size
	DP			OFFSET x20M213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	376100000, 20971520, 2.700,	26828928
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			81920+79*128	;; scratch area size
	DP			OFFSET x20M113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	376100000, 20971520, 2.700,	26828928
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET x20M013
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_0+1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			196608+95*128	;; scratch area size
	DP			OFFSET x24M213AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
allfft	DD	448000000, 25165824, 3.300,	32104704
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x24M113AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_1+1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET x24M013AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			196608+95*128	;; scratch area size
	DP			OFFSET x24M213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET x24M113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET x24M013
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			229376+111*128	;; scratch area size
	DP			OFFSET x28M213AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
allfft	DD	521500000, 29360128, 3.900,	37380480
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET x28M113AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_1+1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			114688+111*128	;; scratch area size
	DP			OFFSET x28M013AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			229376+111*128	;; scratch area size
	DP			OFFSET x28M213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			114688+111*128	;; scratch area size
	DP			OFFSET x28M113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			114688+111*128	;; scratch area size
	DP			OFFSET x28M013
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
allfft	DD	596000000, 33554432, 4.500,	42639360
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			262144+127*128	;; scratch area size
allfft	DP			OFFSET x32M213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET x32M113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET x32M013AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			262144+127*128	;; scratch area size
	DP			OFFSET x32M213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET x32M113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET x32M013
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	0
xjmptablep DD	739,	32,	0.00000111,	1152
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x32p
	DD			0		;; FFT levels done in pass2
	DD			4, 4
	DD			2, 1, 1, 1, 1, 0
	DD	1095,	48,	0.00000144,	1920
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x48p
	DD			0		;; FFT levels done in pass2
	DD			6, 6
	DD			3, 1, 1, 1, 1, 0
	DD	1465,	64,	0.00000178,	2432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x64p
	DD			0		;; FFT levels done in pass2
	DD			8, 8
	DD			4, 1, 1, 1, 1, 0
	DD	2173,	96,	0.00000259,	3328
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x96p
	DD			0		;; FFT levels done in pass2
	DD			12, 12
	DD			6, 3, 1, 1, 1, 0
	DD	2897,	128,	0.00000319,	4352
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x128p
	DD			0		;; FFT levels done in pass2
	DD			16, 16
	DD			8, 4, 1, 2, 1, 0
	DD	4295,	192,	0.00000542,	6784
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x192p
	DD			0		;; FFT levels done in pass2
	DD			24, 24
	DD			12, 6, 3, 3, 1, 0
	DD	5729,	256,	0.00000691,	8576
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x256p
	DD			0		;; FFT levels done in pass2
	DD			32, 32
	DD			16, 8, 4, 4, 1, 0
	DD	8493,	384,	0.0000111,	13312
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x384p
	DD			0		;; FFT levels done in pass2
	DD			48, 48
	DD			24, 12, 3, 6, 1, 0
	DD	11319,	512,	0.0000143,	17152
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x512p
	DD			0		;; FFT levels done in pass2
	DD			64, 64
	DD			32, 16, 4, 2*256+8, 1, 0
	DD	16779,	768,	0.0000260,	26624
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x768p
	DD			0		;; FFT levels done in pass2
	DD			96, 96
	DD			48, 24, 3*256+6, 3*256+12, 1, 0
	DD	22381,	1024,	0.0000349,	34304
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x1024p
	DD			0		;; FFT levels done in pass2
	DD			128, 128
	DD			64, 32, 4*256+8, 4*256+16, 1, 0
	DD	33189,	1536,	0.0000601,	52992
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x1536p
	DD			0		;; FFT levels done in pass2
	DD			192, 192
	DD			96, 48, 12, 6*256+24, 3, 0
	DD	44221,	2048,	0.0000773,	68608
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x2048p
	DD			0		;; FFT levels done in pass2
	DD			256, 256
	DD			128, 64, 16, 8*256+32, 2*256+4, 0
	DD	65519,	3072,	0.000131,	106112
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x3072p
	DD			0		;; FFT levels done in pass2
	DD			384, 384
	DD			192, 96, 6*256+24, 12*256+48, 3*256+3, 0
	DD	87271,	4096,	0.000172,	137216
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x4096p
	DD			0		;; FFT levels done in pass2
	DD			512, 512
	DD			256, 128, 8*256+32, 16*256+64, 4*256+4, 0
	DD	129600,	6144,	0.000291,	211968
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x6144p
	DD			0		;; FFT levels done in pass2
	DD			768, 768
	DD			384, 192, 12*256+48, 24*256+96, 6*256+3, 0
	DD	172400,	8192,	0.000395,	274432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x8192p
	DD			0		;; FFT levels done in pass2
	DD			1024, 1024
	DD			512, 256, 16*256+64, 32*256+128, 8*256+4, 0
	DD	254800,	12288,	0.000626,	64896
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x12KpAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	254800,	12288,	0.000626,	64896
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x12Kp
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	339100,	16384,	0.000857,	84224
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x16KpAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	339100,	16384,	0.000857,	84224
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x16Kp
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	503400,	24576,	0.00135,	123904
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x24KpAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	503400,	24576,	0.00135,	123904
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x24Kp
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	669600,	32768,	0.00179,	162816
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x32KpAMD
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	669600,	32768,	0.00179,	162816
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x32Kp
	DD			8		;; FFT levels done in pass2
	DD			GAP2_8_4+1, 0
	DD	992800, 49152,	0.00303,	145280
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x48KpAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	992800, 49152,	0.00303,	145280
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x48Kp
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	1320000, 65536,	0.00404,	178560
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x64KpAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	1320000, 65536,	0.00404,	178560
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x64Kp
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	1962000, 98304, 0.00644,	245632
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x96KpAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	1962000, 98304, 0.00644,	245632
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x96Kp
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	2610000, 131072, 0.00871,	312064
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x128KpAMD
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	2610000, 131072, 0.00871,	312064
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET x128Kp
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
allfft	DD	3867000, 196608, 0.0136,	445952
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x192KpAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	3867000, 196608, 0.0136,	529536
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x192Kp410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
	DD	3867000, 196608, 0.0136,	529536
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET x192Kp210AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	3867000, 196608, 0.0136,	529536
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			6144+5*128	;; scratch area size
allfft	DP			OFFSET x192Kp110AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	3867000, 196608, 0.0136,	445952
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x192Kp
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	3867000, 196608, 0.0136,	529536
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x192Kp410
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	3867000, 196608, 0.0136,	529536
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			12288+5*128	;; scratch area size
allfft	DP			OFFSET x192Kp210
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	3867000, 196608, 0.0136,	529536
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			6144+5*128	;; scratch area size
	DP			OFFSET x192Kp110
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	5146000, 262144, 0.0179,	579072
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x256KpAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	5146000, 262144, 0.0179,	697344
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x256Kp410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	5146000, 262144, 0.0179,	697344
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET x256Kp210AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	5146000, 262144, 0.0179,	697344
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET x256Kp110AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	5146000, 262144, 0.0179,	579072
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x256Kp
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	5146000, 262144, 0.0179,	697344
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			32768+7*128	;; scratch area size
	DP			OFFSET x256Kp410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
	DD	5146000, 262144, 0.0179,	697344
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			16384+7*128	;; scratch area size
	DP			OFFSET x256Kp210
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	5146000, 262144, 0.0179,	697344
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			8192+7*128	;; scratch area size
allfft	DP			OFFSET x256Kp110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	7635000, 393216, 0.0273,	846976
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x384Kp411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	7635000, 393216, 0.0273,	846976
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x384KpAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	7635000, 393216, 0.0273,	1036288
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x384Kp410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	7635000, 393216, 0.0273,	1036288
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET x384Kp210AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	7635000, 393216, 0.0273,	1036288
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET x384Kp110AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	7635000, 393216, 0.0273,	846976
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x384Kp411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	7635000, 393216, 0.0273,	846976
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x384Kp
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	7635000, 393216, 0.0273,	1036288
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET x384Kp410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
	DD	7635000, 393216, 0.0273,	1036288
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			24576+11*128	;; scratch area size
	DP			OFFSET x384Kp210
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	7635000, 393216, 0.0273,	1036288
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET x384Kp110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	10150000, 524288, 0.0368,	1113088
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x512Kp411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	10150000, 524288, 0.0368,	1113088
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x512KpAMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	10150000, 524288, 0.0368,	1372160
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET x512Kp410AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	10150000, 524288, 0.0368,	1372160
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x512Kp210AMD
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
	DD	10150000, 524288, 0.0368,	1372160
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x512Kp110AMD
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_1+1, 0
allfft	DD	10150000, 524288, 0.0368,	1113088
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x512Kp411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	10150000, 524288, 0.0368,	1113088
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET x512Kp
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	10150000, 524288, 0.0368,	1372160
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET x512Kp410
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_4+1, 0
	DD	10150000, 524288, 0.0368,	1372160
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET x512Kp210
	DD			10		;; FFT levels done in pass2
	DD			GAP2_10_2+1, 0
allfft	DD	10150000, 524288, 0.0368,	1372160
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x512Kp110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	15040000, 786432, 0.0566,	1285248
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x768Kp412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
	DD	15040000, 786432, 0.0566,	1285248
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET x768Kp212AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
allfft	DD	15040000, 786432, 0.0566,	1285248
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			6144+5*128	;; scratch area size
allfft	DP			OFFSET x768Kp112AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_1+1, 0
allfft	DD	15040000, 786432, 0.0566,	1648640
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x768Kp411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	15040000, 786432, 0.0566,	1285248
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET x768Kp412
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
	DD	15040000, 786432, 0.0566,	1285248
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET x768Kp212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	15040000, 786432, 0.0566,	1648640
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET x768Kp411
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	15040000, 786432, 0.0566,	1648640
	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
	DD			24576+11*128	;; scratch area size
	DP			OFFSET x768Kp211
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_2+1, 0
allfft	DD	15040000, 786432, 0.0566,	1648640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET x768Kp111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	15040000, 786432, 0.0566,	1648640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET x768Kp9111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
	DD	15040000, 786432, 0.0566,	1648640
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET x768Kp011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	15040000, 786432, 0.0566,	1648640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET x768Kp9011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	15040000, 786432, 0.0566,	2050176
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x768Kp410
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_4+1, 0
allfft	DD	15040000, 786432, 0.0566,	2050176
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x768Kp210
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_2+1, 0
allfft	DD	15040000, 786432, 0.0566,	2050176
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET x768Kp110
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_1+1, 0
allfft	DD	15040000, 786432, 0.0566,	2050176
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET x768Kp010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	15040000, 786432, 0.0566,	3790464
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x768Kp18
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			GAP2_8_1+1, 0
allfft	DD	20000000, 1048576, 0.0761,	1682432
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET x1024Kp412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	20000000, 1048576, 0.0761,	1682432
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET x1024Kp212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	20000000, 1048576, 0.0761,	1682432
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET x1024Kp112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET x1024Kp411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			RPFW+256*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x1024Kp211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
	DD	20000000, 1048576, 0.0761,	1682432
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			32768+7*128	;; scratch area size
	DP			OFFSET x1024Kp412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	20000000, 1048576, 0.0761,	1682432
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET x1024Kp212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	20000000, 1048576, 0.0761,	2181120
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET x1024Kp411
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_4+1, 0
	DD	20000000, 1048576, 0.0761,	1682432
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET x1024Kp112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
	DD	20000000, 1048576, 0.0761,	2181120
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET x1024Kp211
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_2+1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x1024Kp111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x1024Kp9111_1
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
	DD	20000000, 1048576, 0.0761,	2181120
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x1024Kp011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET x1024Kp9011_1
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	29640000, 1572864, 0.125,	2480128
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x1536Kp412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	29640000, 1572864, 0.125,	2480128
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET x1536Kp212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	29640000, 1572864, 0.125,	2480128
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET x1536Kp112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x1536Kp411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x1536Kp211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET x1536Kp111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x1536Kp411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	29640000, 1572864, 0.125,	2480128
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+23*128	;; scratch area size
	DP			OFFSET x1536Kp412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	29640000, 1572864, 0.125,	2480128
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET x1536Kp212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	29640000, 1572864, 0.125,	2480128
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			12288+23*128	;; scratch area size
	DP			OFFSET x1536Kp112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x1536Kp211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
	DD	29640000, 1572864, 0.125,	3252352
	DD			WILLI*65536+2*1	;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET x1536Kp111
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_1+1, 0
	DD	29640000, 1572864, 0.125,	3252352
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET x1536Kp011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	39390000, 2097152, 0.169,	3274752
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET x2048Kp412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	39390000, 2097152, 0.169,	3274752
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x2048Kp212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	39390000, 2097152, 0.169,	3274752
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x2048Kp112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET x2048Kp411AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET x2048Kp211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x2048Kp111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x2048Kp011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET x2048Kp411
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_4+1, 0
	DD	39390000, 2097152, 0.169,	3274752
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET x2048Kp412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	39390000, 2097152, 0.169,	3274752
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x2048Kp212
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	39390000, 2097152, 0.169,	3274752
	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x2048Kp112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET x2048Kp211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET x2048Kp111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
	DD	39390000, 2097152, 0.169,	4317184
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET x2048Kp011
	DD			11		;; FFT levels done in pass2
	DD			GAP2_11_0+1, 0
allfft	DD	58410000, 3145728, 0.289,	4143104
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x3072Kp413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	58410000, 3145728, 0.289,	4143104
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET x3072Kp213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	58410000, 3145728, 0.289,	4143104
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET x3072Kp113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	58410000, 3145728, 0.289,	4870272
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET x3072Kp412AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_4+1, 0
allfft	DD	58410000, 3145728, 0.289,	4870272
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET x3072Kp212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	58410000, 3145728, 0.289,	4870272
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET x3072Kp112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET x3072Kp211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x3072Kp111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x3072Kp011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	58410000, 3145728, 0.289,	4143104
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET x3072Kp413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	58410000, 3145728, 0.289,	4870272
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			98304+23*128	;; scratch area size
	DP			OFFSET x3072Kp412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	58410000, 3145728, 0.289,	4870272
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			49152+23*128	;; scratch area size
	DP			OFFSET x3072Kp212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	58410000, 3145728, 0.289,	4870272
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET x3072Kp112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET x3072Kp211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x3072Kp111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET x3072Kp011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	58410000, 3145728, 0.289,	8132736
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x3072Kp010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
allfft	DD	77700000, 4194304, 0.425,	5462016
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET x4096Kp413AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
allfft	DD	77700000, 4194304, 0.425,	5462016
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET x4096Kp213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	77700000, 4194304, 0.425,	5462016
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET x4096Kp113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	77700000, 4194304, 0.425,	6459392
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET x4096Kp212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	77700000, 4194304, 0.425,	6459392
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET x4096Kp112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0

	DD	77700000, 4194304, 0.425,	16459392
	DD			RPFW+V26+2*1	;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET xv26_4096Kp112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0

allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET x4096Kp211AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x4096Kp111AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x4096Kp011AMD
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0

	DD	77700000, 4194304, 0.425,	66459392
	DD			V26+1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			131072+31*128	;; scratch area size
	DP			OFFSET xv26_4096Kp412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	77700000, 4194304, 0.425,	66459392
	DD			V26+512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			65536+31*128	;; scratch area size
	DP			OFFSET xv26_4096Kp212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	77700000, 4194304, 0.425,	66459392
	DD			V26+2*1		;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET xv26_4096Kp112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0

	DD	77700000, 4194304, 0.425,	6459392
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			131072+31*128	;; scratch area size
	DP			OFFSET x4096Kp412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	77700000, 4194304, 0.425,	6459392
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			65536+31*128	;; scratch area size
	DP			OFFSET x4096Kp212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	77700000, 4194304, 0.425,	6459392
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET x4096Kp112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0

allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET x4096Kp211
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_2+1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x4096Kp111
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_1+1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x4096Kp011
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			GAP2_11_0+1, 0
allfft	DD	77700000, 4194304, 0.425,	4*6459392
allfft	DD			2*1	;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET x4096Kp010
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			GAP2_10_0+1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET x6Mp412AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	115000000, 6291456, 0.668,	9650176
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET x6Mp212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET x6Mp112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET x6Mp412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			98304+47*128	;; scratch area size
	DP			OFFSET x6Mp212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET x6Mp112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET x8Mp412AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
allfft	DD	153100000, 8388608, 1.042,	12828672
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET x8Mp212AMD
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			GAP2_12_2+1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET x8Mp112AMD
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
allfft	DD	153100000, 8388608, 1.042,	3*12828672
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET x8Mp413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET x8Mp412
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_4+1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			131072+63*128	;; scratch area size
	DP			OFFSET x8Mp212
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_2+1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET x8Mp112
	DD			12		;; FFT levels done in pass2
	DD			GAP2_12_1+1, 0
	DD	226800000, 12582912, 1.400,	16031744
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET x12Mp413AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
allfft	DD	226800000, 12582912, 1.400,	16031744
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET x12Mp213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	226800000, 12582912, 1.400,	16031744
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET x12Mp113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	226800000, 12582912, 1.400,	16031744
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			196608+47*128	;; scratch area size
allfft	DP			OFFSET x12Mp413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	226800000, 12582912, 1.400,	16031744
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			98304+47*128	;; scratch area size
	DP			OFFSET x12Mp213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	226800000, 12582912, 1.400,	16031744
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET x12Mp113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
	DD	301900000, 16777216, 2.100,	21307392
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET x16Mp413AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_4+1, 0
allfft	DD	301900000, 16777216, 2.100,	21307392
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET x16Mp213AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_2+1, 0
	DD	301900000, 16777216, 2.100,	21307392
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET x16Mp113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	301900000, 16777216, 2.100,	21307392
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			262144+63*128	;; scratch area size
allfft	DP			OFFSET x16Mp413
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_4+1, 0
	DD	301900000, 16777216, 2.100,	21307392
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			131072+63*128	;; scratch area size
	DP			OFFSET x16Mp213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	301900000, 16777216, 2.100,	21307392
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET x16Mp113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
allfft	DD	301900000, 16777216, 2.100,	21307392
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET x16Mp013
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_0+1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			196608+95*128	;; scratch area size
	DP			OFFSET x24Mp213AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
allfft	DD	446900000, 25165824, 3.300,	31883392
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET x24Mp113AMD
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			GAP2_13_1+1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET x24Mp013AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			196608+95*128	;; scratch area size
	DP			OFFSET x24Mp213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET x24Mp113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET x24Mp013
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET x32Mp113AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET x32Mp013AMD
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			262144+127*128	;; scratch area size
	DP			OFFSET x32Mp213
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_2+1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET x32Mp113
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_1+1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET x32Mp013
	DD			13		;; FFT levels done in pass2
	DD			GAP2_13_0+1, 0
	DD	0

	;; Align so that other GWDATA areas are also aligned on a cache line
	align 128
_GWDATA ENDS

;; FFT setup routines

_TEXT SEGMENT

; gwinfo1 (resptr)
;	Return address of jmp tables for C code to examine
; Windows 32-bit (_gwinfo1)
; Linux 32-bit (gwinfo1)
;	Parameter resptr = [esp+4]
; Windows 64-bit (gwinfo1) - leaf routine, no unwind info necessary
;	Parameter resptr = rcx
; Linux 64-bit (gwinfo1)
;	Parameter resptr = rdi

PROCL	gwinfo1
	IFNDEF X86_64
	mov	ecx, [esp+4]		; Address of data struct to return info
	ENDIF
	IFDEF LINUX64
	mov	rcx, rdi		; Address of data struct to return info
	ENDIF
	mov	rax, OFFSET xjmptable	; P4 mersenne mod FFTs
	mov	[rcx+0*SZPTR], rax
	mov	rax, OFFSET xjmptablep	; P4 2^N+1 mod FFTs
	mov	[rcx+1*SZPTR], rax
	IFNDEF X86_64
	mov	rax, OFFSET jmptable	; x86 mersenne mod FFTs
	mov	[rcx+2*SZPTR], rax
	mov	rax, OFFSET jmptablep	; x86 2^N+1 mod FFTs
	mov	[rcx+3*SZPTR], rax
	ENDIF
	mov	eax, VERSION_NUMBER
	mov	[rcx+4*SZPTR], eax
	ret
gwinfo1 ENDP

_TEXT	ENDS
END
