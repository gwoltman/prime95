; Copyright 1995-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the setup, common routines, and global variables
; for the various discrete-weighted transforms
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

VERSION_NUMBER = 2411		;; Version 24.11

INCLUDE	unravel.mac

EXTRN	_CPU_FLAGS:DWORD
EXTRN	_KARG:QWORD
EXTRN	_PARG:DWORD
EXTRN	_CARG:DWORD
EXTRN	_FFTLEN:DWORD
EXTRN	_RATIONAL_FFT:DWORD
EXTRN	_BITS_PER_WORD:DWORD
EXTRN	_ZERO_PADDED_FFT:DWORD
EXTRN	_INFT:EXTPTR
EXTRN	_GWVERNUM:DWORD
EXTRN	_SRCARG:EXTPTR
EXTRN	_SRC2ARG:EXTPTR
EXTRN	_DESTARG:EXTPTR
EXTRN	_DEST2ARG:EXTPTR
EXTRN	_NUMARG:DWORD
EXTRN	_NUM2ARG:DWORD
EXTRN	_GWPROCPTRS:DWORD

IFNDEF X86_64
EXTRN	gwadd1:PROC
EXTRN	gwaddq1:PROC
EXTRN	gwsub1:PROC
EXTRN	gwsubq1:PROC
EXTRN	gwaddsub1:PROC
EXTRN	gwaddsubq1:PROC
EXTRN	gwcopyzero1:PROC
EXTRN	gwadd2:PROC
EXTRN	gwaddq2:PROC
EXTRN	gwsub2:PROC
EXTRN	gwsubq2:PROC
EXTRN	gwaddsub2:PROC
EXTRN	gwaddsubq2:PROC
EXTRN	gwcopyzero2:PROC
ENDIF

EXTRN	gwxadd1:PROC
EXTRN	gwxaddq1:PROC
EXTRN	gwxsub1:PROC
EXTRN	gwxsubq1:PROC
EXTRN	gwxaddsub1:PROC
EXTRN	gwxaddsubq1:PROC
EXTRN	gwxcopyzero1:PROC
EXTRN	gwxaddf1:PROC
EXTRN	gwxsubf1:PROC
EXTRN	gwxaddsubf1:PROC
EXTRN	gwxadd2:PROC
EXTRN	gwxaddq2:PROC
EXTRN	gwxsub2:PROC
EXTRN	gwxsubq2:PROC
EXTRN	gwxaddsub2:PROC
EXTRN	gwxaddsubq2:PROC
EXTRN	gwxcopyzero2:PROC
EXTRN	gwxaddf2:PROC
EXTRN	gwxsubf2:PROC
EXTRN	gwxaddsubf2:PROC
EXTRN	gwxadd3:PROC
EXTRN	gwxaddq3:PROC
EXTRN	gwxsub3:PROC
EXTRN	gwxsubq3:PROC
EXTRN	gwxaddsub3:PROC
EXTRN	gwxaddsubq3:PROC
EXTRN	gwxcopyzero3:PROC
EXTRN	gwxaddf3:PROC
EXTRN	gwxsubf3:PROC
EXTRN	gwxaddsubf3:PROC

;; List of possible suffixes for the normalization routines

exnorm	MACRO num, prefix, suffix
	exnorm1	&prefix, &num&suffix
	exnorm1	&prefix, &num&e&suffix
	exnorm1	&prefix, &num&c&suffix
	exnorm1	&prefix, &num&ec&suffix
	exnorm1	&prefix, &num&z&suffix
	exnorm1	&prefix, &num&ze&suffix
	exnorm1	&prefix, &num&zp&suffix
	exnorm1	&prefix, &num&zpe&suffix
	exnorm1	&prefix, &num&zpc&suffix
	exnorm1	&prefix, &num&zpec&suffix
	ENDM
exnorm1	MACRO prefix, suffix
	EXTRN	&prefix&r&suffix:PROC
	EXTRN	&prefix&i&suffix:PROC
	ENDM

IFNDEF X86_64
exnorm 1
exnorm 2
ENDIF
exnorm 1, x
exnorm 2, x
exnorm 2, x, AMD
exnorm 3, x
exnorm 3, x, AMD

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
SFX_P4_CLM210	EQU	SFX_P4_CLM2 + SFX_P4_CLM1 + SFX_P4_CLM0
SFX_P4_CLM42	EQU	SFX_P4_CLM4 + SFX_P4_CLM2
SFX_P4_CLM421	EQU	SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM1
SFX_P4_CLM4210	EQU	SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM1 + SFX_P4_CLM0

SFX_AMD64	EQU	400h	;; AMD64 optimized (SSE2, AMD suffix)
SFX_AMD64_CLM0	EQU	800h
SFX_AMD64_CLM1	EQU	1000h
SFX_AMD64_CLM2	EQU	2000h
SFX_AMD64_CLM4	EQU	4000h
SFX_AMD64_CLM210 EQU	SFX_AMD64_CLM2 + SFX_AMD64_CLM1 + SFX_AMD64_CLM0
SFX_AMD64_CLM42	EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM2
SFX_AMD64_CLM421 EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM2 + SFX_AMD64_CLM1
SFX_AMD64_CLM4210 EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM2 + SFX_AMD64_CLM1 + SFX_AMD64_CLM0

exfft	MACRO fft_length, x
	exfft1	fft_length, _1, x, 0, 1, 2, 4
	exfft1	fft_length, _2, x, 0, 1, 2, 4
	exfft1	fft_length, _3, x, 0, 1, 2, 4
	exfft1	fft_length, _4, x, 0, 1, 2, 4
	ENDM
exfft1	MACRO fft_length, suffix, x, c0, c1, c2, c4
	IFNDEF X86_64
	IF x AND SFX_486
	EXTRN	fft&fft_length&suffix:PROC
	ENDIF
	IF x AND SFX_PPRO
	EXTRN	fft&fft_length&suffix&PPRO:PROC
	ENDIF
	IF x AND SFX_PPRO_CLM2
	EXTRN	fft&fft_length&c2&suffix&PPRO:PROC
	ENDIF
	IF x AND SFX_P3
	EXTRN	fft&fft_length&suffix&P3:PROC
	ENDIF
	IF x AND SFX_P3_CLM2
	EXTRN	fft&fft_length&c2&suffix&P3:PROC
	ENDIF
	ENDIF
	IF x AND SFX_P4
	EXTRN	xfft&fft_length&suffix:PROC
	ENDIF
	IF x AND SFX_P4_CLM0
	EXTRN	xfft&fft_length&c0&suffix:PROC
	ENDIF
	IF x AND SFX_P4_CLM1
	EXTRN	xfft&fft_length&c1&suffix:PROC
	ENDIF
	IF x AND SFX_P4_CLM2
	EXTRN	xfft&fft_length&c2&suffix:PROC
	ENDIF
	IF x AND SFX_P4_CLM4
	EXTRN	xfft&fft_length&c4&suffix:PROC
	ENDIF
	IF x AND SFX_AMD64
	EXTRN	xfft&fft_length&suffix&AMD:PROC
	ENDIF
	IF x AND SFX_AMD64_CLM0
	EXTRN	xfft&fft_length&c0&suffix&AMD:PROC
	ENDIF
	IF x AND SFX_AMD64_CLM1
	EXTRN	xfft&fft_length&c1&suffix&AMD:PROC
	ENDIF
	IF x AND SFX_AMD64_CLM2
	EXTRN	xfft&fft_length&c2&suffix&AMD:PROC
	ENDIF
	IF x AND SFX_AMD64_CLM4
	EXTRN	xfft&fft_length&c4&suffix&AMD:PROC
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
exfft	160K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	192K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	224K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	256K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	320K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	384K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	448K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	512K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	640K, SFX_PPRO + SFX_P3 + SFX_P4_CLM42 + SFX_AMD64_CLM42
exfft	768K, SFX_PPRO + SFX_P3 + SFX_P4_CLM42 + SFX_AMD64_CLM42
exfft	896K, SFX_PPRO + SFX_P3 + SFX_P4_CLM421 + SFX_AMD64_CLM421
exfft	1024K, SFX_PPRO_CLM2 + SFX_P3_CLM2 + SFX_P4_CLM421 + SFX_AMD64_CLM421
exfft	1280K, SFX_PPRO + SFX_P3 + SFX_P4_CLM421 + SFX_AMD64_CLM421
exfft	1536K, SFX_PPRO + SFX_P3 + SFX_P4_CLM421 + SFX_AMD64_CLM421
exfft	1792K, SFX_PPRO + SFX_P3 + SFX_P4_CLM4210 + SFX_AMD64_CLM4210
exfft	2048K, SFX_PPRO_CLM2 + SFX_P3_CLM2 + SFX_P4_CLM4210 + SFX_AMD64_CLM4210
exfft	2560K, SFX_PPRO_CLM2 + SFX_P3_CLM2 + SFX_P4_CLM210 + SFX_AMD64_CLM210
exfft	3072K, SFX_PPRO_CLM2 + SFX_P3_CLM2 + SFX_P4_CLM210 + SFX_AMD64_CLM210
exfft	3584K, SFX_PPRO_CLM2 + SFX_P3_CLM2 + SFX_P4_CLM210 + SFX_AMD64_CLM210
exfft	4096K, SFX_PPRO_CLM2 + SFX_P3_CLM2 + SFX_P4_CLM210 + SFX_AMD64_CLM210

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
exfft	192Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	256Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	384Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	512Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
exfft	768Kp, SFX_PPRO + SFX_P3 + SFX_P4_CLM42 + SFX_AMD64_CLM42
exfft	1024Kp, SFX_PPRO + SFX_P3 + SFX_P4_CLM421 + SFX_AMD64_CLM421
exfft	1536Kp, SFX_PPRO + SFX_P3 + SFX_P4_CLM421 + SFX_AMD64_CLM421
exfft	2048Kp, SFX_PPRO + SFX_P3 + SFX_P4_CLM4210 + SFX_AMD64_CLM4210
exfft	3072Kp, SFX_PPRO_CLM2 + SFX_P3_CLM2 + SFX_P4_CLM210 + SFX_AMD64_CLM210
exfft	4096Kp, SFX_PPRO_CLM2 + SFX_P3_CLM2 + SFX_P4_CLM210 + SFX_AMD64_CLM210

PUBLIC	MINUS_CARG
PUBLIC	SQRTHALF
PUBLIC	HALF
PUBLIC	BIGVAL
PUBLIC	BIGBIGVAL
PUBLIC	ttmp_ff_inv
PUBLIC	P309
PUBLIC	M809
PUBLIC	M262
PUBLIC	M382
PUBLIC	P951
PUBLIC	P588
PUBLIC	M162
PUBLIC	P618
PUBLIC	P623
PUBLIC	M358
PUBLIC	P404
PUBLIC	P975
PUBLIC	P445
PUBLIC	P180
PUBLIC	M223
PUBLIC	M901
PUBLIC	M691
PUBLIC	P866
PUBLIC	P433
PUBLIC	P577
PUBLIC	P25
PUBLIC	P75
PUBLIC	P3
PUBLIC	pass2_premults
PUBLIC	plus1_premults
PUBLIC	addcount1
PUBLIC	normcount1
PUBLIC	count1
PUBLIC	count2
PUBLIC	count3
PUBLIC	count4
PUBLIC	count5
PUBLIC	xsincos_complex
PUBLIC	sincos1
PUBLIC	sincos2
PUBLIC	sincos3
PUBLIC	sincos4
PUBLIC	sincos5
PUBLIC	sincos6
PUBLIC	sincos7
PUBLIC	sincos8
PUBLIC	sincos9
PUBLIC	sincos10
PUBLIC	norm_grp_mults
PUBLIC	norm_col_mults
PUBLIC	norm_biglit_array
PUBLIC	carries
PUBLIC	scratch_area
PUBLIC	zero_fft
PUBLIC	const_fft
PUBLIC	zpad_addr
PUBLIC	ffttype
PUBLIC	cache_line_multiplier
PUBLIC	pass1blkdst
PUBLIC	normblkdst
PUBLIC	normblkdst8
PUBLIC	loopcount1
PUBLIC	loopcount2
PUBLIC	loopcount3
PUBLIC	loopcount4
PUBLIC	norm_ptr1
PUBLIC	norm_ptr2
PUBLIC	norm_ptr3
PUBLIC	normval1
PUBLIC	normval2
PUBLIC	normval3
PUBLIC	normval4

PUBLIC	XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
PUBLIC	XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8
PUBLIC	XMM_TWO
PUBLIC	XMM_HALF
PUBLIC	XMM_SQRTHALF
PUBLIC	XMM_BIGVAL
PUBLIC	XMM_BIGBIGVAL
PUBLIC	XMM_MINUS_CARG
PUBLIC	XMM_MINUS_CARG_TIMES_MULCONST
PUBLIC	XMM_K_LO
PUBLIC	XMM_K_TIMES_MULCONST_LO
PUBLIC	XMM_K_HI
PUBLIC	XMM_K_TIMES_MULCONST_HI
PUBLIC	XMM_K_HI_1
PUBLIC	XMM_K_TIMES_MULCONST_HI_1
PUBLIC	XMM_K_HI_2
PUBLIC	XMM_K_TIMES_MULCONST_HI_2
PUBLIC	XMM_ZPAD0, XMM_ZPAD1, XMM_ZPAD2, XMM_ZPAD3
PUBLIC	XMM_ZPAD4, XMM_ZPAD5, XMM_ZPAD6
PUBLIC	XMM_ABSVAL
PUBLIC	XMM_LIMIT_BIGMAX
PUBLIC	XMM_LIMIT_BIGMAX_NEG
PUBLIC	XMM_LIMIT_INVERSE
PUBLIC	XMM_TTP_FUDGE
PUBLIC	XMM_TTMP_FUDGE
PUBLIC	XMM_MAXERR
PUBLIC	XMM_SUMOUT
PUBLIC	XMM_NORM012_FF
PUBLIC	XMM_MULCONST

PUBLIC	XMM_P309
PUBLIC	XMM_M809
PUBLIC	XMM_M262
PUBLIC	XMM_M382
PUBLIC	XMM_P951
PUBLIC	XMM_P588
PUBLIC	XMM_M162
PUBLIC	XMM_P618

PUBLIC	XMM_P25
PUBLIC	XMM_P75
PUBLIC	XMM_P3
PUBLIC	XMM_P433
PUBLIC	XMM_P577
PUBLIC	XMM_P866

PUBLIC	XMM_P623
PUBLIC	XMM_M358
PUBLIC	XMM_P404
PUBLIC	XMM_P975
PUBLIC	XMM_P445
PUBLIC	XMM_P180

;
; Global variables needed in multiplication routines
;

IFNDEF X86_64
_GWDATA SEGMENT PAGE PUBLIC 'DATA'
ELSE
_GWDATA SEGMENT PAGE
ENDIF

; These are the SSE2 floating point constants that are computed in C code
	
XMM_LIMIT_BIGMAX DQ	32 DUP (0.0)	; Normalization constants
XMM_LIMIT_BIGMAX_NEG DQ	32 DUP (0.0)
XMM_LIMIT_INVERSE DQ	32 DUP (0.0)
XMM_SQRTHALF	DQ	0.5, 0.5
XMM_BIGVAL	DQ	0.0, 0.0	; Used to round double to integer
XMM_BIGBIGVAL	DQ	0.0, 0.0	; Used to round to big word multiple
XMM_NORM012_FF	DQ	0.0, 0.0	; Used in xnorm012 macros (FFTLEN/2)
XMM_MINUS_CARG	DQ	0.0, 0.0	; CARG negated and stored as double
XMM_K_LO	DQ	0.0, 0.0	; KARG, bottom big word bits
XMM_K_HI	DQ	0.0, 0.0	; KARG, top bits
XMM_K_HI_1	DQ	0.0, 0.0	; XMM_K_HI, bottom big word bits
XMM_K_HI_2	DQ	0.0, 0.0	; XMM_K_HI, top bits
XMM_MINUS_CARG_TIMES_MULCONST DQ 0.0, 0.0
XMM_K_TIMES_MULCONST_LO DQ 0.0, 0.0	; KARG*mulconst, bottom big word bits
XMM_K_TIMES_MULCONST_HI DQ 0.0, 0.0	; KARG*mulconst, top bits
XMM_K_TIMES_MULCONST_HI_1 DQ 0.0, 0.0	; XMM_K_TIMES_MULCONST_HI, low bits
XMM_K_TIMES_MULCONST_HI_2 DQ 0.0, 0.0	; XMM_K_TIMES_MULCONST_HI, top bits
XMM_P309	DQ	0.309, 0.309	; Used in five reals pfa
XMM_M809	DQ	-0.809, -0.809
XMM_M262	DQ	-0.262, -0.262
XMM_M382	DQ	-0.382, -0.382
XMM_P951	DQ	0.951, 0.951
XMM_P588	DQ	0.588, 0.588
XMM_M162	DQ	-0.162, -0.162
XMM_P618	DQ	0.618, 0.618
XMM_P25		DQ	0.25, 0.25	; Used in six reals pfa
XMM_P75		DQ	0.75, 0.75
XMM_P3		DQ	3.0, 3.0
XMM_P433	DQ	0.433, 0.433
XMM_P577	DQ	0.577, 0.577
XMM_P866	DQ	0.866, 0.866
XMM_P623	DQ	0.623, 0.623	; Used in seven reals pfa
XMM_M358	DQ	-0.358, -0.358
XMM_P404	DQ	4.04, 4.04
XMM_P975	DQ	0.975, 0.975
XMM_P445	DQ	0.445, 0.445
XMM_P180	DQ	1.80, 1.80

; SSE2 temporaries and constants

XMM_ZPAD0	DQ	0.0		; Zero padding temporaries
XMM_ZPAD1	DQ	0.0
XMM_ZPAD2	DQ	0.0
XMM_ZPAD3	DQ	0.0
XMM_ZPAD4	DQ	0.0
XMM_ZPAD5	DQ	0.0
XMM_ZPAD6	DQ	0.0
XMM_ZPAD7	DQ	0.0
XMM_TMP1	DQ	0.0, 0.0
XMM_TMP2	DQ	0.0, 0.0
XMM_TMP3	DQ	0.0, 0.0
XMM_TMP4	DQ	0.0, 0.0
XMM_TMP5	DQ	0.0, 0.0
XMM_TMP6	DQ	0.0, 0.0
XMM_TMP7	DQ	0.0, 0.0
XMM_TMP8	DQ	0.0, 0.0
XMM_ABSVAL	DD	7FFFFFFFh,0FFFFFFFFh,7FFFFFFFh,0FFFFFFFFh
					; Used to compute absolute values
XMM_TWO		DQ	2.0, 2.0
XMM_HALF	DQ	0.5, 0.5
XMM_SUMOUT	DQ	0.0, 0.0	; Used in normalization macros
XMM_MAXERR	DQ	0.0, 0.0	; Used in normalization macros
XMM_MULCONST	DQ	0.0, 0.0
XMM_TTP_FUDGE	DQ	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
		DQ	0.5, 1.0, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0
		DQ	1.0, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0, 0.5
		DQ	0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5
XMM_TTMP_FUDGE	DQ	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
		DQ	2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0
		DQ	1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0
		DQ	2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0

; These are the x87 floating point constants that are computed in C code
	
MINUS_CARG	DQ	0.0		; CARG negated and stored as double
SQRTHALF	DQ	0.0		; Used in all ffts
ttmp_ff_inv	DQ	0.0		; Inverse FFT adjust (2/FFTLEN)
BIGVAL		DD	0.0		; Used to round to an integer
BIGBIGVAL	DD	0.0		; Used to round to multiple of big word 
P309		DQ	0.309		; Used in five_reals_fft/unfft
M809		DQ	-0.809		; Used in five_reals_fft/unfft
M262		DQ	-2.618		; Used in five_reals_fft/unfft
M382		DQ	-0.382		; Used in five_reals_fft/unfft
P951		DQ	0.951		; Used in five_reals_fft/unfft
P588		DQ	0.588		; Used in five_reals_fft/unfft
M162		DQ	-1.617		; Used in five_reals_fft/unfft
P618		DQ	0.618		; Used in five_reals_fft/unfft
P623		DQ	0.623		; Used in seven_reals_fft/unfft
M358		DQ	-0.358		; Used in seven_reals_fft/unfft
P404		DQ	4.040		; Used in seven_reals_fft/unfft
P975		DQ	0.975		; Used in seven_reals_fft/unfft
P445		DQ	0.445		; Used in seven_reals_fft/unfft
P180		DQ	1.802		; Used in seven_reals_fft/unfft
M223		DQ	-0.223		; Used in seven_reals_fft/unfft
M901		DQ	-0.901		; Used in seven_reals_fft/unfft
M691		DQ	-0.691		; Used in seven_reals_fft/unfft
P866		DQ	0.866		; Used in six_reals_fft/unfft
P433		DQ	0.433		; Used in six_reals_fft/unfft
P577		DQ	0.577		; Used in six_reals_fft/unfft

; Floating point constants in x87 code not initialized in C code

P25		DD	0.25		; Used in six_reals_fft/unfft
P75		DD	0.75		; Used in six_reals_fft/unfft
P3		DD	3.0		; Used in six_reals_fft/unfft
HALF		DD	0.5		; Used in all ffts

; Counters and variables

pass2_premults	DPTR	0		; Address of pass 2 premultiplier data
xsincos_complex	DPTR	0		; Addr of pass2 complex sin/cos data
sincos1		DPTR	0
sincos2		DPTR	0
sincos3		DPTR	0
sincos4		DPTR	0
sincos5		DPTR	0
sincos6		DPTR	0
sincos7		DPTR	0
sincos8		DPTR	0
sincos9		DPTR	0
sincos10	DPTR	0
norm_grp_mults	DPTR	0		; Ptr to array of normalize multipliers
norm_col_mults	DPTR	0
norm_biglit_array DPTR	0		; Ptr to byte array of big/lit flags
carries		DPTR	0		; Ptr to array of carries (2 pass FFT)
scratch_area	DPTR	0		; Scratch area for pass 1 of SSE2 FFTs
plus1_premults	DPTR	0		; Address of 2^N+1 premultiplier data
zpad_addr	DPTR	0		; Address of XMM_ZPAD0
pass1blkdst	DPTR	0		; Dist between blocks in pass 1
normblkdst	DPTR	0		; Dist between blocks in normalization
normblkdst8	DPTR	0		; Dist between 8 blocks in normalize
norm_ptr1	DPTR	0
norm_ptr2	DPTR	0
norm_ptr3	DPTR	0
normval2	DPTR	0
normval3	DPTR	0

addcount1	DD	0		; Loop counters used in adding
normcount1	DD	0		; Loop counters used in normalizing
count1		DD	0		; Counters used in common fft code
count2		DD	0
count3		DD	0
count4		DD	0
count5		DD	0
zero_fft	DD	0		; TRUE if zero upper half normalize
const_fft	DD	0		; TRUE if mul-by-const in normalize
ffttype		DD	0		; Type of fft (1, 2, 3, or 4)
cache_line_multiplier DD 0		; Cache line multiplier from jmptable
loopcount1	DD	0
loopcount2	DD	0
loopcount3	DD	0
loopcount4	DD	0
normval1	DD	0
normval4	DD	0

RPF	EQU	80000000h	; FFT requires prefetch capability
RPFW	EQU	40000000h	; FFT requires prefetchw (3DNow!) capability

IFNDEF	X86_64
jmptable DD	755,	32,	0.0000036,	672
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft32_1PPRO, OFFSET fft32_2PPRO
	DD			OFFSET fft32_3PPRO, OFFSET fft32_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			3, 1, 1, 1, 0
	DD	939,	40,	0.0000057,	948
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft40_1PPRO, OFFSET fft40_2PPRO
	DD			OFFSET fft40_3PPRO, OFFSET fft40_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			4, 1, 1, 1, 0
	DD	1113,	48,	0.0000065,	1128
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft48_1PPRO, OFFSET fft48_2PPRO
	DD			OFFSET fft48_3PPRO, OFFSET fft48_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			5, 1, 1, 1, 0
	DD	1303,	56,	0.0000084,	1356
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft56_1PPRO, OFFSET fft56_2PPRO
	DD			OFFSET fft56_3PPRO, OFFSET fft56_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			6, 1, 1, 1, 0
	DD	1499,	64,	0.0000083,	1392
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft64_1PPRO, OFFSET fft64_2PPRO
	DD			OFFSET fft64_3PPRO, OFFSET fft64_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			7, 1, 1, 1, 0
	DD	1857,	80,	0.0000121,	1848
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft80_1PPRO, OFFSET fft80_2PPRO
	DD			OFFSET fft80_3PPRO, OFFSET fft80_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			9, 2, 1, 1, 0
	DD	2211,	96,	0.0000141,	2208
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft96_1PPRO, OFFSET fft96_2PPRO
	DD			OFFSET fft96_3PPRO, OFFSET fft96_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			11, 2, 1, 1, 0
	DD	2585,	112,	0.0000179,	2616
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft112_1PPRO, OFFSET fft112_2PPRO
	DD			OFFSET fft112_3PPRO, OFFSET fft112_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			13, 3, 1, 1, 0
	DD	2953,	128,	0.0000178,	2832
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft128_1PPRO, OFFSET fft128_2PPRO
	DD			OFFSET fft128_3PPRO, OFFSET fft128_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			15, 3, 1, 1, 0
	DD	3663,	160,	0.0000296,	3840
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft160_1PPRO, OFFSET fft160_2PPRO
	DD			OFFSET fft160_3PPRO, OFFSET fft160_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			19, 4, 2, 1, 0
	DD	4359,	192,	0.000035,	4608
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft192_1PPRO, OFFSET fft192_2PPRO
	DD			OFFSET fft192_3PPRO, OFFSET fft192_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			23, 5, 2, 1, 0
	DD	5093,	224,	0.000045,	5424
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft224_1PPRO, OFFSET fft224_2PPRO
	DD			OFFSET fft224_3PPRO, OFFSET fft224_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			27, 6, 3, 1, 0
	DD	5833,	256,	0.000045,	5712
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft256_1PPRO, OFFSET fft256_2PPRO
	DD			OFFSET fft256_3PPRO, OFFSET fft256_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			31, 7, 3, 1, 0
	DD	7243,	320,	0.000062,	7680
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft320_1PPRO, OFFSET fft320_2PPRO
	DD			OFFSET fft320_3PPRO, OFFSET fft320_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			39, 9, 2, 1, 0
	DD	8639,	384,	0.000075,	9216
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft384_1PPRO, OFFSET fft384_2PPRO
	DD			OFFSET fft384_3PPRO, OFFSET fft384_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			47, 11, 2, 1, 0
	DD	10085,	448,	0.000093,	10800
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft448_1PPRO, OFFSET fft448_2PPRO
	DD			OFFSET fft448_3PPRO, OFFSET fft448_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			55, 13, 3, 1, 0
	DD	11537,	512,	0.000097,	11472
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft512_1PPRO, OFFSET fft512_2PPRO
	DD			OFFSET fft512_3PPRO, OFFSET fft512_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			63, 15, 3, 1, 0
	DD	14301,	640,	0.000140,	15552
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft640_1PPRO, OFFSET fft640_2PPRO
	DD			OFFSET fft640_3PPRO, OFFSET fft640_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			79, 19, 4, 2, 0
	DD	17047,	768,	0.000167,	18672
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft768_1PPRO, OFFSET fft768_2PPRO
	DD			OFFSET fft768_3PPRO, OFFSET fft768_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			95, 23, 5, 2, 0
	DD	19881,	896,	0.000210,	21840
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft896_1PPRO, OFFSET fft896_2PPRO
	DD			OFFSET fft896_3PPRO, OFFSET fft896_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			111, 27, 6, 3, 0
	DD	22799,	1024,	0.000218,	22992
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft1024_1PPRO, OFFSET fft1024_2PPRO
	DD			OFFSET fft1024_3PPRO, OFFSET fft1024_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			127, 31, 7, 3, 0
	DD	28295,	1280,	0.000302,	31152
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft1280_1PPRO, OFFSET fft1280_2PPRO
	DD			OFFSET fft1280_3PPRO, OFFSET fft1280_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			159, 39, 9, 2, 0
	DD	33761,	1536,	0.000365,	37392
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft1536_1PPRO, OFFSET fft1536_2PPRO
	DD			OFFSET fft1536_3PPRO, OFFSET fft1536_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			191, 47, 11, 2, 0
	DD	39411,	1792,	0.000456,	43680
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft1792_1PPRO, OFFSET fft1792_2PPRO
	DD			OFFSET fft1792_3PPRO, OFFSET fft1792_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			223, 55, 13, 3, 0
	DD	45061,	2048,	0.000490,	46032
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft2048_1PPRO, OFFSET fft2048_2PPRO
	DD			OFFSET fft2048_3PPRO, OFFSET fft2048_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			255, 63, 15, 3, 0
	DD	55825,	2560,	0.000708,	62544
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft2560_1PPRO, OFFSET fft2560_2PPRO
	DD			OFFSET fft2560_3PPRO, OFFSET fft2560_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			319, 79, 19, 4, 2, 0
	DD	66519,	3072,	0.000851,	75072
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft3072_1PPRO, OFFSET fft3072_2PPRO
	DD			OFFSET fft3072_3PPRO, OFFSET fft3072_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			383, 95, 23, 5, 2, 0
	DD	77599,	3584,	0.00107,	87648
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft3584_1PPRO, OFFSET fft3584_2PPRO
	DD			OFFSET fft3584_3PPRO, OFFSET fft3584_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			447, 111, 27, 6, 3, 0
	DD	89047,	4096,	0.00113,	92112
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft4096_1PPRO, OFFSET fft4096_2PPRO
	DD			OFFSET fft4096_3PPRO, OFFSET fft4096_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			511, 127, 31, 7, 3, 0
	DD	110400,	5120,	0.00152,	24848
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft5120_1PPRO, OFFSET fft5120_2PPRO
	DD			OFFSET fft5120_3PPRO, OFFSET fft5120_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10, 1
	DD			9, 1, 1, 1, 0
	DD	131100,	6144,	0.00191,	28016
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft6144_1PPRO, OFFSET fft6144_2PPRO
	DD			OFFSET fft6144_3PPRO, OFFSET fft6144_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12, 1
	DD			11, 1, 1, 1, 0
	DD	152800,	7168,	0.00226,	31232
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft7168_1PPRO, OFFSET fft7168_2PPRO
	DD			OFFSET fft7168_3PPRO, OFFSET fft7168_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			14, 1
	DD			13, 1, 1, 1, 0
	DD	175300,	8192,	0.00242,	34400
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft8192_1PPRO, OFFSET fft8192_2PPRO
	DD			OFFSET fft8192_3PPRO, OFFSET fft8192_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			16, 1
	DD			15, 1, 1, 1, 0
	DD	217700,	10240,	0.00333,	40880
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft10K_1PPRO, OFFSET fft10K_2PPRO
	DD			OFFSET fft10K_3PPRO, OFFSET fft10K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			20, 1
	DD			19, 1, 1, 1, 0
	DD	258200,	12288,	0.00397,	47264
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft12K_1PPRO, OFFSET fft12K_2PPRO
	DD			OFFSET fft12K_3PPRO, OFFSET fft12K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			24, 1
	DD			23, 1, 1, 1, 0
	DD	301400,	14336,	0.00488,	53696
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft14K_1PPRO, OFFSET fft14K_2PPRO
	DD			OFFSET fft14K_3PPRO, OFFSET fft14K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			28, 1
	DD			27, 1, 1, 1, 0
	DD	346100,	16384,	0.00522,	59936
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft16K_1PPRO, OFFSET fft16K_2PPRO
	DD			OFFSET fft16K_3PPRO, OFFSET fft16K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			32, 1
	DD			31, 1, 1, 1, 0
	DD	430300,	20480,	0.00692,	72800
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft20K_1PPRO, OFFSET fft20K_2PPRO
	DD			OFFSET fft20K_3PPRO, OFFSET fft20K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			40, 1
	DD			39, 1, 1, 1, 0
	DD	511600,	24576,	0.00826,	85568
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft24K_1PPRO, OFFSET fft24K_2PPRO
	DD			OFFSET fft24K_3PPRO, OFFSET fft24K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			48, 1
	DD			47, 1, 1, 1, 0
	DD	596100,	28672,	0.0101,		98384
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft28K_1PPRO, OFFSET fft28K_2PPRO
	DD			OFFSET fft28K_3PPRO, OFFSET fft28K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			56, 1
	DD			55, 1, 1, 1, 0
	DD	683700,	32768,	0.0109,		111008
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft32K_1PPRO, OFFSET fft32K_2PPRO
	DD			OFFSET fft32K_3PPRO, OFFSET fft32K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			64, 1
	DD			63, 1, 1, 1, 0
	DD	848800,	40960,	0.0151,		136832
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft40K_1P3, OFFSET fft40K_2P3
	DD			OFFSET fft40K_3P3, OFFSET fft40K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			80, 1
	DD			79, 1, 1, 1, 0
	DD	848800,	40960,	0.0151,		136832
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft40K_1PPRO, OFFSET fft40K_2PPRO
	DD			OFFSET fft40K_3PPRO, OFFSET fft40K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			80, 1
	DD			79, 1, 1, 1, 0
	DD	1009000, 49152,	0.0184,		162416
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft48K_1P3, OFFSET fft48K_2P3
	DD			OFFSET fft48K_3P3, OFFSET fft48K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			96, 1
	DD			95, 1, 1, 1, 0
	DD	1009000, 49152,	0.0184,		162416
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft48K_1PPRO, OFFSET fft48K_2PPRO
	DD			OFFSET fft48K_3PPRO, OFFSET fft48K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			96, 1
	DD			95, 1, 1, 1, 0
	DD	1177000, 57344,	0.0227,		188048
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft56K_1P3, OFFSET fft56K_2P3
	DD			OFFSET fft56K_3P3, OFFSET fft56K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, 1, 1, 0
	DD	1177000, 57344,	0.0227,		188048
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft56K_1PPRO, OFFSET fft56K_2PPRO
	DD			OFFSET fft56K_3PPRO, OFFSET fft56K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, 1, 1, 0
	DD	1350000, 65536,	0.0252,		213152
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft64K_1P3, OFFSET fft64K_2P3
	DD			OFFSET fft64K_3P3, OFFSET fft64K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 1, 1, 0
	DD	1350000, 65536,	0.0252,		213152
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft64K_1PPRO, OFFSET fft64K_2PPRO
	DD			OFFSET fft64K_3PPRO, OFFSET fft64K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 1, 1, 0
	DD	1678000, 81920, 0.0360,		264752
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft80K_1P3, OFFSET fft80K_2P3
	DD			OFFSET fft80K_3P3, OFFSET fft80K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	1678000, 81920, 0.0360,		264752
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft80K_1PPRO, OFFSET fft80K_2PPRO
	DD			OFFSET fft80K_3PPRO, OFFSET fft80K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	1994000, 98304, 0.0445,		297536
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft96K_1P3, OFFSET fft96K_2P3
	DD			OFFSET fft96K_3P3, OFFSET fft96K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			48, 1
	DD			47, 1, 1, 1, 0
	DD	1994000, 98304, 0.0445,		297536
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft96K_1PPRO, OFFSET fft96K_2PPRO
	DD			OFFSET fft96K_3PPRO, OFFSET fft96K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			48, 1
	DD			47, 1, 1, 1, 0
	DD	2324000, 114688, 0.0548,	341072
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft112K_1P3, OFFSET fft112K_2P3
	DD			OFFSET fft112K_3P3, OFFSET fft112K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			56, 1
	DD			55, 1, 1, 1, 0
	DD	2324000, 114688, 0.0548,	341072
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft112K_1PPRO, OFFSET fft112K_2PPRO
	DD			OFFSET fft112K_3PPRO, OFFSET fft112K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			56, 1
	DD			55, 1, 1, 1, 0
	DD	2664000, 131072, 0.0604,	384416
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft128K_1P3, OFFSET fft128K_2P3
	DD			OFFSET fft128K_3P3, OFFSET fft128K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			64, 1
	DD			63, 1, 1, 1, 0
	DD	2664000, 131072, 0.0604,	384416
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft128K_1PPRO, OFFSET fft128K_2PPRO
	DD			OFFSET fft128K_3PPRO, OFFSET fft128K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			64, 1
	DD			63, 1, 1, 1, 0
	DD	3310000, 163840, 0.0830,	471680+5120
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			5120		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft160K_1P3, OFFSET fft160K_2P3
	DD			OFFSET fft160K_3P3, OFFSET fft160K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			80, 1
	DD			79, 1, 1, 1, 0
	DD	3310000, 163840, 0.0830,	471680+5120
	DD			2		;; Flags, min_l2_cache, clm
	DD			5120		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft160K_1PPRO, OFFSET fft160K_2PPRO
	DD			OFFSET fft160K_3PPRO, OFFSET fft160K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			80, 1
	DD			79, 1, 1, 1, 0
	DD	3933000, 196608, 0.0982,	558704+6144
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft192K_1P3, OFFSET fft192K_2P3
	DD			OFFSET fft192K_3P3, OFFSET fft192K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			96, 1
	DD			95, 1, 1, 1, 0
	DD	3933000, 196608, 0.0982,	558704+6144
	DD			2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft192K_1PPRO, OFFSET fft192K_2PPRO
	DD			OFFSET fft192K_3PPRO, OFFSET fft192K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			96, 1
	DD			95, 1, 1, 1, 0
	DD	4593000, 229376, 0.1193,	645776+7168
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			7168		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft224K_1P3, OFFSET fft224K_2P3
	DD			OFFSET fft224K_3P3, OFFSET fft224K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, 1, 1, 0
	DD	4593000, 229376, 0.1193,	645776+7168
	DD			2		;; Flags, min_l2_cache, clm
	DD			7168		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft224K_1PPRO, OFFSET fft224K_2PPRO
	DD			OFFSET fft224K_3PPRO, OFFSET fft224K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, 1, 1, 0
	DD	5264000, 262144, 0.1316,	732320+8192+3*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft256K_1P3, OFFSET fft256K_2P3
	DD			OFFSET fft256K_3P3, OFFSET fft256K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 1, 1, 0
	DD	5264000, 262144, 0.1316,	732320+8192+3*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft256K_1PPRO, OFFSET fft256K_2PPRO
	DD			OFFSET fft256K_3PPRO, OFFSET fft256K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 1, 1, 0
	DD	6545000, 327680, 0.1726,	906800+10240+4*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft320K_1P3, OFFSET fft320K_2P3
	DD			OFFSET fft320K_3P3, OFFSET fft320K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	6545000, 327680, 0.1726,	906800+10240+4*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft320K_1PPRO, OFFSET fft320K_2PPRO
	DD			OFFSET fft320K_3PPRO, OFFSET fft320K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	7772000, 393216, 0.2107,	1080848+12288+5*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft384K_1P3, OFFSET fft384K_2P3
	DD			OFFSET fft384K_3P3, OFFSET fft384K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 1, 1, 0
	DD	7772000, 393216, 0.2107,	1080848+12288+5*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft384K_1PPRO, OFFSET fft384K_2PPRO
	DD			OFFSET fft384K_3PPRO, OFFSET fft384K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 1, 1, 0
	DD	9071000, 458752, 0.2520,	1254944+14336+6*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft448K_1P3, OFFSET fft448K_2P3
	DD			OFFSET fft448K_3P3, OFFSET fft448K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, 1, 1, 0
	DD	9071000, 458752, 0.2520,	1254944+14336+6*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft448K_1PPRO, OFFSET fft448K_2PPRO
	DD			OFFSET fft448K_3PPRO, OFFSET fft448K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, 1, 1, 0
	DD	10380000, 524288, 0.2808,	1428128+16384+7*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft512K_1P3, OFFSET fft512K_2P3
	DD			OFFSET fft512K_3P3, OFFSET fft512K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 1, 1, 0
	DD	10380000, 524288, 0.2808,	1428128+16384+7*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft512K_1PPRO, OFFSET fft512K_2PPRO
	DD			OFFSET fft512K_3PPRO, OFFSET fft512K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 1, 1, 0
	DD	12890000, 655360, 0.372,	1777232+20480+9*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft640K_1P3, OFFSET fft640K_2P3
	DD			OFFSET fft640K_3P3, OFFSET fft640K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 1, 1, 0
	DD	12890000, 655360, 0.372,	1777232+20480+9*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft640K_1PPRO, OFFSET fft640K_2PPRO
	DD			OFFSET fft640K_3PPRO, OFFSET fft640K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 1, 1, 0
	DD	15310000, 786432, 0.453,	2125376+24576+11*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft768K_1P3, OFFSET fft768K_2P3
	DD			OFFSET fft768K_3P3, OFFSET fft768K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 1, 1, 0
	DD	15310000, 786432, 0.453,	2125376+24576+11*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft768K_1PPRO, OFFSET fft768K_2PPRO
	DD			OFFSET fft768K_3PPRO, OFFSET fft768K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 1, 1, 0
	DD	17890000, 917504, 0.536,	2473568+28672+13*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft896K_1P3, OFFSET fft896K_2P3
	DD			OFFSET fft896K_3P3, OFFSET fft896K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, 1, 1, 0
	DD	17890000, 917504, 0.536,	2473568+28672+13*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft896K_1PPRO, OFFSET fft896K_2PPRO
	DD			OFFSET fft896K_3PPRO, OFFSET fft896K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, 1, 1, 0
	DD	20460000, 1048576, 0.600,	2819744+32768+15*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft1024K2_1P3, OFFSET fft1024K2_2P3
	DD			OFFSET fft1024K2_3P3, OFFSET fft1024K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 1, 1, 0
	DD	20460000, 1048576, 0.600,	2819744+32768+15*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft1024K2_1PPRO, OFFSET fft1024K2_2PPRO
	DD			OFFSET fft1024K2_3PPRO, OFFSET fft1024K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 1, 1, 0
	DD	25390000, 1310720, 0.776,	3474992+10240+4*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft1280K_1P3, OFFSET fft1280K_2P3
	DD			OFFSET fft1280K_3P3, OFFSET fft1280K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	25390000, 1310720, 0.776,	3474992+10240+4*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft1280K_1PPRO, OFFSET fft1280K_2PPRO
	DD			OFFSET fft1280K_3PPRO, OFFSET fft1280K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	30190000, 1572864, 0.934,	4140560+12288+5*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft1536K_1P3, OFFSET fft1536K_2P3
	DD			OFFSET fft1536K_3P3, OFFSET fft1536K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 1, 1, 0
	DD	30190000, 1572864, 0.934,	4140560+12288+5*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft1536K_1PPRO, OFFSET fft1536K_2PPRO
	DD			OFFSET fft1536K_3PPRO, OFFSET fft1536K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 1, 1, 0
	DD	35200000, 1835008, 1.113,	4806176+14336+6*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft1792K_1P3, OFFSET fft1792K_2P3
	DD			OFFSET fft1792K_3P3, OFFSET fft1792K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, 1, 1, 0
	DD	35200000, 1835008, 1.113,	4806176+14336+6*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft1792K_1PPRO, OFFSET fft1792K_2PPRO
	DD			OFFSET fft1792K_3PPRO, OFFSET fft1792K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, 1, 1, 0
	DD	40300000, 2097152, 1.226,	5470880+16384+7*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft2048K2_1P3, OFFSET fft2048K2_2P3
	DD			OFFSET fft2048K2_3P3, OFFSET fft2048K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 1, 1, 0
	DD	40300000, 2097152, 1.226,	5470880+16384+7*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft2048K2_1PPRO, OFFSET fft2048K2_2PPRO
	DD			OFFSET fft2048K2_3PPRO, OFFSET fft2048K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 1, 1, 0
	DD	50020000, 2621440, 1.636,	6803024+20480+9*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft2560K2_1P3, OFFSET fft2560K2_2P3
	DD			OFFSET fft2560K2_3P3, OFFSET fft2560K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 1, 1, 0
	DD	50020000, 2621440, 1.636,	6803024+20480+9*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft2560K2_1PPRO, OFFSET fft2560K2_2PPRO
	DD			OFFSET fft2560K2_3PPRO, OFFSET fft2560K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 1, 1, 0
	DD	59510000, 3145728, 1.990,	8134208+24576+11*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft3072K2_1P3, OFFSET fft3072K2_2P3
	DD			OFFSET fft3072K2_3P3, OFFSET fft3072K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 1, 1, 0
	DD	59510000, 3145728, 1.990,	8134208+24576+11*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft3072K2_1PPRO, OFFSET fft3072K2_2PPRO
	DD			OFFSET fft3072K2_3PPRO, OFFSET fft3072K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 1, 1, 0
	DD	69360000, 3670016, 2.380,	9465440+28672+13*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft3584K2_1P3, OFFSET fft3584K2_2P3
	DD			OFFSET fft3584K2_3P3, OFFSET fft3584K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, 1, 1, 0
	DD	69360000, 3670016, 2.380,	9465440+28672+13*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft3584K2_1PPRO, OFFSET fft3584K2_2PPRO
	DD			OFFSET fft3584K2_3PPRO, OFFSET fft3584K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, 1, 1, 0
	DD	79370000, 4194304, 2.604,	10794656+32768+15*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft4096K2_1P3, OFFSET fft4096K2_2P3
	DD			OFFSET fft4096K2_3P3, OFFSET fft4096K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 1, 1, 0
	DD	79370000, 4194304, 2.604,	10794656+32768+15*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft4096K2_1PPRO, OFFSET fft4096K2_2PPRO
	DD			OFFSET fft4096K2_3PPRO, OFFSET fft4096K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 1, 1, 0
	DD	0
jmptablep DD	755,	32,	0.000004,	976
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft32p_1PPRO, OFFSET fft32p_2PPRO
	DD			OFFSET fft32p_3PPRO, OFFSET fft32p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			4, 1, 1, 1, 0
	DD	1111,	48,	0.000007,	1608
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft48p_1PPRO, OFFSET fft48p_2PPRO
	DD			OFFSET fft48p_3PPRO, OFFSET fft48p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			6, 3, 1, 1, 0
	DD	1485,	64,	0.000010,	1952
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft64p_1PPRO, OFFSET fft64p_2PPRO
	DD			OFFSET fft64p_3PPRO, OFFSET fft64p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			8, 4, 1, 1, 0
	DD	2199,	96,	0.0000141,	3072
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft96p_1PPRO, OFFSET fft96p_2PPRO
	DD			OFFSET fft96p_3PPRO, OFFSET fft96p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			12, 3, 1, 1, 0
	DD	2947,	128,	0.000021,	3904
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft128p_1PPRO, OFFSET fft128p_2PPRO
	DD			OFFSET fft128p_3PPRO, OFFSET fft128p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			16, 4, 1, 1, 1, 0
	DD	4345,	192,	0.000035,	6288
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft192p_1PPRO, OFFSET fft192p_2PPRO
	DD			OFFSET fft192p_3PPRO, OFFSET fft192p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			24, 6, 3, 1, 0
	DD	5817,	256,	0.000051,	7808
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft256p_1PPRO, OFFSET fft256p_2PPRO
	DD			OFFSET fft256p_3PPRO, OFFSET fft256p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			32, 8, 4, 1, 1, 0
	DD	8607,	384,	0.000075,	12432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft384p_1PPRO, OFFSET fft384p_2PPRO
	DD			OFFSET fft384p_3PPRO, OFFSET fft384p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			48, 12, 3, 1, 0
	DD	11515,	512,	0.000106,	15616
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft512p_1PPRO, OFFSET fft512p_2PPRO
	DD			OFFSET fft512p_3PPRO, OFFSET fft512p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			64, 16, 4, 1, 1, 0
	DD	17001,	768,	0.000167,	25008
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft768p_1PPRO, OFFSET fft768p_2PPRO
	DD			OFFSET fft768p_3PPRO, OFFSET fft768p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			96, 24, 6, 3, 1, 0
	DD	22701,	1024,	0.000249,	31232
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft1024p_1PPRO, OFFSET fft1024p_2PPRO
	DD			OFFSET fft1024p_3PPRO, OFFSET fft1024p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			128, 32, 8, 4, 1, 0
	DD	33569,	1536,	0.000365,	49872
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft1536p_1PPRO, OFFSET fft1536p_2PPRO
	DD			OFFSET fft1536p_3PPRO, OFFSET fft1536p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			192, 48, 12, 3, 0
	DD	44951,	2048,	0.000582,	62464
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft2048p_1PPRO, OFFSET fft2048p_2PPRO
	DD			OFFSET fft2048p_3PPRO, OFFSET fft2048p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			256, 64, 16, 4, 1, 0
	DD	66319,	3072,	0.000851,	99888
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft3072p_1PPRO, OFFSET fft3072p_2PPRO
	DD			OFFSET fft3072p_3PPRO, OFFSET fft3072p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			384, 96, 24, 6, 3, 0
	DD	88747,	4096,	0.00135,	124928
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DD			OFFSET fft4096p_1PPRO, OFFSET fft4096p_2PPRO
	DD			OFFSET fft4096p_3PPRO, OFFSET fft4096p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			1, 1
	DD			512, 128, 32, 8, 4, 0
	DD	130600,	6144,	0.00191,	26512
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft6144p_1PPRO, OFFSET fft6144p_2PPRO
	DD			OFFSET fft6144p_3PPRO, OFFSET fft6144p_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12, 1
	DD			12, 1, 1, 1, 0
	DD	174000,	8192,	0.00284,	32960
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft8192p_1PPRO, OFFSET fft8192p_2PPRO
	DD			OFFSET fft8192p_3PPRO, OFFSET fft8192p_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			16, 1
	DD			16, 1, 1, 1, 0
	DD	257700,	12288,	0.00397,	46000
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft12Kp_1PPRO, OFFSET fft12Kp_2PPRO
	DD			OFFSET fft12Kp_3PPRO, OFFSET fft12Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			24, 1
	DD			24, 1, 1, 1, 0
	DD	344700,	16384,	0.00588,	58752
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft16Kp_1PPRO, OFFSET fft16Kp_2PPRO
	DD			OFFSET fft16Kp_3PPRO, OFFSET fft16Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			32, 1
	DD			32, 1, 1, 1, 0
	DD	508600,	24576,	0.00826,	84688
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft24Kp_1PPRO, OFFSET fft24Kp_2PPRO
	DD			OFFSET fft24Kp_3PPRO, OFFSET fft24Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			48, 1
	DD			48, 1, 1, 1, 0
	DD	679400,	32768,	0.01299,	110336
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft32Kp_1PPRO, OFFSET fft32Kp_2PPRO
	DD			OFFSET fft32Kp_3PPRO, OFFSET fft32Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			64, 1
	DD			64, 1, 1, 1, 0
	DD	1006000, 49152,	0.0184,		162352
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft48Kp_1P3, OFFSET fft48Kp_2P3
	DD			OFFSET fft48Kp_3P3, OFFSET fft48Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			96, 1
	DD			96, 1, 1, 1, 0
	DD	1006000, 49152,	0.0184,		162352
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft48Kp_1PPRO, OFFSET fft48Kp_2PPRO
	DD			OFFSET fft48Kp_3PPRO, OFFSET fft48Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			96, 1
	DD			96, 1, 1, 1, 0
	DD	1345000, 65536,	0.03283,	213504
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft64Kp_1P3, OFFSET fft64Kp_2P3
	DD			OFFSET fft64Kp_3P3, OFFSET fft64Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 1, 1, 0
	DD	1345000, 65536,	0.03283,	213504
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DD			OFFSET fft64Kp_1PPRO, OFFSET fft64Kp_2PPRO
	DD			OFFSET fft64Kp_3PPRO, OFFSET fft64Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 1, 1, 0
	DD	1983000, 98304, 0.0445,		290512
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft96Kp_1P3, OFFSET fft96Kp_2P3
	DD			OFFSET fft96Kp_3P3, OFFSET fft96Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			48, 1
	DD			48, 1, 1, 1, 0
	DD	1983000, 98304, 0.0445,		290512
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft96Kp_1PPRO, OFFSET fft96Kp_2PPRO
	DD			OFFSET fft96Kp_3PPRO, OFFSET fft96Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			48, 1
	DD			48, 1, 1, 1, 0
	DD	2652000, 131072, 0.0719,	377600
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft128Kp_1P3, OFFSET fft128Kp_2P3
	DD			OFFSET fft128Kp_3P3, OFFSET fft128Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			64, 1
	DD			64, 1, 1, 1, 0
	DD	2652000, 131072, 0.0719,	377600
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft128Kp_1PPRO, OFFSET fft128Kp_2PPRO
	DD			OFFSET fft128Kp_3PPRO, OFFSET fft128Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			64, 1
	DD			64, 1, 1, 1, 0
	DD	3924000, 196608, 0.0982,	552496+6144
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft192Kp_1P3, OFFSET fft192Kp_2P3
	DD			OFFSET fft192Kp_3P3, OFFSET fft192Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			96, 1
	DD			96, 1, 1, 1, 0
	DD	3924000, 196608, 0.0982,	552496+6144
	DD			2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft192Kp_1PPRO, OFFSET fft192Kp_2PPRO
	DD			OFFSET fft192Kp_3PPRO, OFFSET fft192Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			96, 1
	DD			96, 1, 1, 1, 0
	DD	5242000, 262144, 0.155,		726528+8192+3*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft256Kp_1P3, OFFSET fft256Kp_2P3
	DD			OFFSET fft256Kp_3P3, OFFSET fft256Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 1, 1, 0
	DD	5242000, 262144, 0.155,		726528+8192+3*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft256Kp_1PPRO, OFFSET fft256Kp_2PPRO
	DD			OFFSET fft256Kp_3PPRO, OFFSET fft256Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 1, 1, 0
	DD	7733000, 393216, 0.2107,	1076176+12288+5*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft384Kp_1P3, OFFSET fft384Kp_2P3
	DD			OFFSET fft384Kp_3P3, OFFSET fft384Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 1, 1, 0
	DD	7733000, 393216, 0.2107,	1076176+12288+5*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft384Kp_1PPRO, OFFSET fft384Kp_2PPRO
	DD			OFFSET fft384Kp_3PPRO, OFFSET fft384Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 1, 1, 0
	DD	10320000, 524288, 0.322,	1424384+16384+7*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft512Kp_1P3, OFFSET fft512Kp_2P3
	DD			OFFSET fft512Kp_3P3, OFFSET fft512Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 1, 1, 0
	DD	10320000, 524288, 0.322,	1424384+16384+7*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft512Kp_1PPRO, OFFSET fft512Kp_2PPRO
	DD			OFFSET fft512Kp_3PPRO, OFFSET fft512Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 1, 1, 0
	DD	15260000, 786432, 0.453,	2123824+24576+11*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft768Kp_1P3, OFFSET fft768Kp_2P3
	DD			OFFSET fft768Kp_3P3, OFFSET fft768Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 1, 1, 0
	DD	15260000, 786432, 0.453,	2123824+24576+11*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft768Kp_1PPRO, OFFSET fft768Kp_2PPRO
	DD			OFFSET fft768Kp_3PPRO, OFFSET fft768Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 1, 1, 0
	DD	20360000, 1048576, 0.681,	2820096+32768+15*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft1024Kp_1P3, OFFSET fft1024Kp_2P3
	DD			OFFSET fft1024Kp_3P3, OFFSET fft1024Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			512, 1
	DD			512, 1, 1, 1, 0
	DD	20360000, 1048576, 0.681,	2820096+32768+15*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			10		;; FFT levels done in pass2
	DD			OFFSET fft1024Kp_1PPRO, OFFSET fft1024Kp_2PPRO
	DD			OFFSET fft1024Kp_3PPRO, OFFSET fft1024Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 1, 1, 0
	DD	30070000, 1572864, 0.934,	4111312+24576+5*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+5*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft1536Kp_1P3, OFFSET fft1536Kp_2P3
	DD			OFFSET fft1536Kp_3P3, OFFSET fft1536Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 1, 1, 0
	DD	30070000, 1572864, 0.934,	4111312+24576+5*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+5*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft1536Kp_1PPRO, OFFSET fft1536Kp_2PPRO
	DD			OFFSET fft1536Kp_3PPRO, OFFSET fft1536Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 1, 1, 0
	DD	40110000, 2097152, 1.380,	5442560+32768+7*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+7*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft2048Kp_1P3, OFFSET fft2048Kp_2P3
	DD			OFFSET fft2048Kp_3P3, OFFSET fft2048Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 1, 1, 0
	DD	40110000, 2097152, 1.380,	5442560+32768+7*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+7*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft2048Kp_1PPRO, OFFSET fft2048Kp_2PPRO
	DD			OFFSET fft2048Kp_3PPRO, OFFSET fft2048Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 1, 1, 0
	DD	59360000, 3145728, 1.990,	8108080+24576+11*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft3072Kp2_1P3, OFFSET fft3072Kp2_2P3
	DD			OFFSET fft3072Kp2_3P3, OFFSET fft3072Kp2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 1, 1, 0
	DD	59360000, 3145728, 1.990,	8108080+24576+11*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft3072Kp2_1PPRO, OFFSET fft3072Kp2_2PPRO
	DD			OFFSET fft3072Kp2_3PPRO, OFFSET fft3072Kp2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 1, 1, 0
	DD	79100000, 4194304, 2.919,	10770432+32768+15*64
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft4096Kp2_1P3, OFFSET fft4096Kp2_2P3
	DD			OFFSET fft4096Kp2_3P3, OFFSET fft4096Kp2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			512, 1
	DD			512, 1, 1, 1, 0
	DD	79100000, 4194304, 2.919,	10770432+32768+15*64
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			12		;; FFT levels done in pass2
	DD			OFFSET fft4096Kp2_1PPRO, OFFSET fft4096Kp2_2PPRO
	DD			OFFSET fft4096Kp2_3PPRO, OFFSET fft4096Kp2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			512, 1
	DD			512, 1, 1, 1, 0
	DD	0
ENDIF

IFNDEF X86_64
prctab1	DD	OFFSET gwadd1, OFFSET gwaddq1, OFFSET gwsub1, OFFSET gwsubq1
	DD	OFFSET gwaddsub1, OFFSET gwaddsubq1, OFFSET gwcopyzero1
	DD	OFFSET gwaddq1, OFFSET gwsubq1, OFFSET gwaddsubq1
	DD	OFFSET r1, OFFSET r1e, OFFSET r1c, OFFSET r1ec
	DD	OFFSET r1z, OFFSET r1ze, 0, 0
	DD	OFFSET r1zp, OFFSET r1zpe, OFFSET r1zpc, OFFSET r1zpec
	DD	OFFSET i1, OFFSET i1e, OFFSET i1c, OFFSET i1ec
	DD	OFFSET i1z, OFFSET i1ze, 0, 0
	DD	OFFSET i1zp, OFFSET i1zpe, OFFSET i1zpc, OFFSET i1zpec
prctab2	DD	OFFSET gwadd2, OFFSET gwaddq2, OFFSET gwsub2, OFFSET gwsubq2
	DD	OFFSET gwaddsub2, OFFSET gwaddsubq2, OFFSET gwcopyzero2
	DD	OFFSET gwaddq2, OFFSET gwsubq2, OFFSET gwaddsubq2
	DD	OFFSET r2, OFFSET r2e, OFFSET r2c, OFFSET r2ec
	DD	OFFSET r2z, OFFSET r2ze, 0, 0
	DD	OFFSET r2zp, OFFSET r2zpe, OFFSET r2zpc, OFFSET r2zpec
	DD	OFFSET i2, OFFSET i2e, OFFSET i2c, OFFSET i2ec
	DD	OFFSET i2z, OFFSET i2ze, 0, 0
	DD	OFFSET i2zp, OFFSET i2zpe, OFFSET i2zpc, OFFSET i2zpec
ENDIF

xprctab1 DPTR	OFFSET gwxadd1, OFFSET gwxaddq1, OFFSET gwxsub1
	DPTR	OFFSET gwxsubq1, OFFSET gwxaddsub1, OFFSET gwxaddsubq1
	DPTR	OFFSET gwxcopyzero1
	DPTR	OFFSET gwxaddf1, OFFSET gwxsubf1, OFFSET gwxaddsubf1
	DPTR	OFFSET xr1, OFFSET xr1e, OFFSET xr1c, OFFSET xr1ec
	DPTR	OFFSET xr1z, OFFSET xr1ze, 0, 0
	DPTR	OFFSET xr1zp, OFFSET xr1zpe, OFFSET xr1zpc, OFFSET xr1zpec
	DPTR	OFFSET xi1, OFFSET xi1e, OFFSET xi1c, OFFSET xi1ec
	DPTR	OFFSET xi1z, OFFSET xi1ze, 0, 0
	DPTR	OFFSET xi1zp, OFFSET xi1zpe, OFFSET xi1zpc, OFFSET xi1zpec
xprctab2 DPTR	OFFSET gwxadd2, OFFSET gwxaddq2, OFFSET gwxsub2
	DPTR	OFFSET gwxsubq2, OFFSET gwxaddsub2, OFFSET gwxaddsubq2
	DPTR	OFFSET gwxcopyzero2
	DPTR	OFFSET gwxaddf2, OFFSET gwxsubf2, OFFSET gwxaddsubf2
	DPTR	OFFSET xr2, OFFSET xr2e, OFFSET xr2c, OFFSET xr2ec
	DPTR	OFFSET xr2z, OFFSET xr2ze, 0, 0
	DPTR	OFFSET xr2zp, OFFSET xr2zpe, OFFSET xr2zpc, OFFSET xr2zpec
	DPTR	OFFSET xi2, OFFSET xi2e, OFFSET xi2c, OFFSET xi2ec
	DPTR	OFFSET xi2z, OFFSET xi2ze, 0, 0
	DPTR	OFFSET xi2zp, OFFSET xi2zpe, OFFSET xi2zpc, OFFSET xi2zpec
xprctab2a DPTR	OFFSET gwxadd2, OFFSET gwxaddq2, OFFSET gwxsub2
	DPTR	OFFSET gwxsubq2, OFFSET gwxaddsub2, OFFSET gwxaddsubq2
	DPTR	OFFSET gwxcopyzero2
	DPTR	OFFSET gwxaddf2, OFFSET gwxsubf2, OFFSET gwxaddsubf2
	DPTR	OFFSET xr2AMD, OFFSET xr2eAMD, OFFSET xr2cAMD, OFFSET xr2ecAMD
	DPTR	OFFSET xr2zAMD, OFFSET xr2zeAMD, 0, 0
	DPTR	OFFSET xr2zpAMD, OFFSET xr2zpeAMD
	DPTR	OFFSET xr2zpcAMD, OFFSET xr2zpecAMD
	DPTR	OFFSET xi2AMD, OFFSET xi2eAMD, OFFSET xi2cAMD, OFFSET xi2ecAMD
	DPTR	OFFSET xi2zAMD, OFFSET xi2zeAMD, 0, 0
	DPTR	OFFSET xi2zpAMD, OFFSET xi2zpeAMD
	DPTR	OFFSET xi2zpcAMD, OFFSET xi2zpecAMD
xprctab3 DPTR	OFFSET gwxadd3, OFFSET gwxaddq3, OFFSET gwxsub3
	DPTR	OFFSET gwxsubq3, OFFSET gwxaddsub3, OFFSET gwxaddsubq3
	DPTR	OFFSET gwxcopyzero3
	DPTR	OFFSET gwxaddf3, OFFSET gwxsubf3, OFFSET gwxaddsubf3
	DPTR	OFFSET xr3, OFFSET xr3e, OFFSET xr3c, OFFSET xr3ec
	DPTR	OFFSET xr3z, OFFSET xr3ze, 0, 0
	DPTR	OFFSET xr3zp, OFFSET xr3zpe, OFFSET xr3zpc, OFFSET xr3zpec
	DPTR	OFFSET xi3, OFFSET xi3e, OFFSET xi3c, OFFSET xi3ec
	DPTR	OFFSET xi3z, OFFSET xi3ze, 0, 0
	DPTR	OFFSET xi3zp, OFFSET xi3zpe, OFFSET xi3zpc, OFFSET xi3zpec
xprctab3a DPTR	OFFSET gwxadd3, OFFSET gwxaddq3, OFFSET gwxsub3
	DPTR	OFFSET gwxsubq3, OFFSET gwxaddsub3, OFFSET gwxaddsubq3
	DPTR	OFFSET gwxcopyzero3
	DPTR	OFFSET gwxaddf3, OFFSET gwxsubf3, OFFSET gwxaddsubf3
	DPTR	OFFSET xr3AMD, OFFSET xr3eAMD, OFFSET xr3cAMD, OFFSET xr3ecAMD
	DPTR	OFFSET xr3zAMD, OFFSET xr3zeAMD, 0, 0
	DPTR	OFFSET xr3zpAMD, OFFSET xr3zpeAMD
	DPTR	OFFSET xr3zpcAMD, OFFSET xr3zpecAMD
	DPTR	OFFSET xi3AMD, OFFSET xi3eAMD, OFFSET xi3cAMD, OFFSET xi3ecAMD
	DPTR	OFFSET xi3zAMD, OFFSET xi3zeAMD, 0, 0
	DPTR	OFFSET xi3zpAMD, OFFSET xi3zpeAMD
	DPTR	OFFSET xi3zpcAMD, OFFSET xi3zpecAMD

;; Jump tables for the Pentium 4 SSE2 optimized code

xjmptable DD	743,	32,	0.00000111,	720
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft32_1, OFFSET xfft32_2
	DPTR			OFFSET xfft32_3, OFFSET xfft32_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			4, 4
	DD			1, 1, 1, 1, 1, 0
	DD	1099,	48,	0.00000144,	1272
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft48_1, OFFSET xfft48_2
	DPTR			OFFSET xfft48_3, OFFSET xfft48_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			6, 6
	DD			2, 1, 1, 1, 1, 0
	DD	1469,	64,	0.00000178,	1728
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft64_1, OFFSET xfft64_2
	DPTR			OFFSET xfft64_3, OFFSET xfft64_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			8, 8
	DD			3, 1, 1, 1, 1, 0
	DD	1827,	80,	0.00000222,	1992
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft80_1, OFFSET xfft80_2
	DPTR			OFFSET xfft80_3, OFFSET xfft80_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			10, 8*2048+2
	DD			4, 2, 1, 1, 1, 0
	DD	2179,	96,	0.00000259,	2352
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft96_1, OFFSET xfft96_2
	DPTR			OFFSET xfft96_3, OFFSET xfft96_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			12, 12
	DD			5, 2, 1, 1, 1, 0
	DD	2539,	112,	0.00000311,	2808
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft112_1, OFFSET xfft112_2
	DPTR			OFFSET xfft112_3, OFFSET xfft112_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			14, (8*2048+4)*2048+2
	DD			6, 3, 1, 1, 1, 0
	DD	2905,	128,	0.00000319,	3168
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft128_1, OFFSET xfft128_2
	DPTR			OFFSET xfft128_3, OFFSET xfft128_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			16, 16
	DD			7, 3, 1, 1, 1, 0
	DD	3613,	160,	0.00000450,	4656
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft160_1, OFFSET xfft160_2
	DPTR			OFFSET xfft160_3, OFFSET xfft160_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			20, 16*2048+4
	DD			9, 9, 2, 1, 4, 0
	DD	4311,	192,	0.00000542,	5568
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft192_1, OFFSET xfft192_2
	DPTR			OFFSET xfft192_3, OFFSET xfft192_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			24, 24
	DD			11, 11, 2, 1, 5, 0
	DD	5029,	224,	0.00000663,	6576
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft224_1, OFFSET xfft224_2
	DPTR			OFFSET xfft224_3, OFFSET xfft224_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			28, (16*2048+8)*2048+4
	DD			13, 13, 3, 1, 6, 0
	DD	5755,	256,	0.00000691,	7200
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft256_1, OFFSET xfft256_2
	DPTR			OFFSET xfft256_3, OFFSET xfft256_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			32, 32
	DD			15, 15, 3, 1, 7, 0
	DD	7149,	320,	0.00000928,	8256
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft320_1, OFFSET xfft320_2
	DPTR			OFFSET xfft320_3, OFFSET xfft320_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			40, 32*2048+8
	DD			19, 9, 2, 4, 1, 0
	DD	8527,	384,	0.0000111,	9888
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft384_1, OFFSET xfft384_2
	DPTR			OFFSET xfft384_3, OFFSET xfft384_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			48, 48
	DD			23, 11, 2, 5, 1, 0
	DD	9933,	448,	0.0000133,	11616
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft448_1, OFFSET xfft448_2
	DPTR			OFFSET xfft448_3, OFFSET xfft448_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			56, (32*2048+16)*2048+8
	DD			27, 13, 3, 6, 1, 0
	DD	11359,	512,	0.0000143,	12960
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft512_1, OFFSET xfft512_2
	DPTR			OFFSET xfft512_3, OFFSET xfft512_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			64, 64
	DD			31, 15, 3, 7, 1, 0
	DD	14119,	640,	0.0000215,	17280
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft640_1, OFFSET xfft640_2
	DPTR			OFFSET xfft640_3, OFFSET xfft640_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			80, 64*2048+16
	DD			39, 19, 9, 9, 4*256+2, 0
	DD	16839,	768,	0.0000260,	20736
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft768_1, OFFSET xfft768_2
	DPTR			OFFSET xfft768_3, OFFSET xfft768_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			96, 96
	DD			47, 23, 11, 11, 5*256+2, 0
	DD	19639,	896,	0.0000321,	24288
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft896_1, OFFSET xfft896_2
	DPTR			OFFSET xfft896_3, OFFSET xfft896_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			112, (64*2048+32)*2048+16
	DD			55, 27, 13, 13, 6*256+3, 0
	DD	22477,	1024,	0.0000349,	26016
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024_1, OFFSET xfft1024_2
	DPTR			OFFSET xfft1024_3, OFFSET xfft1024_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			128, 128
	DD			63, 31, 15, 15, 7*256+3, 0
	DD	27899,	1280,	0.0000494,	33600
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft1280_1, OFFSET xfft1280_2
	DPTR			OFFSET xfft1280_3, OFFSET xfft1280_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			160, 128*2048+32
	DD			79, 39, 9, 4*256+19, 2, 0
	DD	33289,	1536,	0.0000601,	40320
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536_1, OFFSET xfft1536_2
	DPTR			OFFSET xfft1536_3, OFFSET xfft1536_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			192, 192
	DD			95, 47, 11, 5*256+23, 2, 0
	DD	38799,	1792,	0.0000719,	47136
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft1792_1, OFFSET xfft1792_2
	DPTR			OFFSET xfft1792_3, OFFSET xfft1792_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			224, (128*2048+64)*2048+32
	DD			111, 55, 13, 6*256+27, 3, 0
	DD	44339,	2048,	0.0000773,	52128
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048_1, OFFSET xfft2048_2
	DPTR			OFFSET xfft2048_3, OFFSET xfft2048_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			256, 256
	DD			127, 63, 15, 7*256+31, 3, 0
	DD	55099,	2560,	0.000111,	68064
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft2560_1, OFFSET xfft2560_2
	DPTR			OFFSET xfft2560_3, OFFSET xfft2560_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			320, 256*2048+64
	DD			159, 79, 9*256+19, 9*256+39, 4*256+2, 0
	DD	65729,	3072,	0.000131,	81696
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072_1, OFFSET xfft3072_2
	DPTR			OFFSET xfft3072_3, OFFSET xfft3072_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			384, 384
	DD			191, 95, 11*256+23, 11*256+47, 5*256+2, 0
	DD	76559,	3584,	0.000165,	95424
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft3584_1, OFFSET xfft3584_2
	DPTR			OFFSET xfft3584_3, OFFSET xfft3584_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			448, (256*2048+128)*2048+64
	DD			223, 111, 13*256+27, 13*256+55, 6*256+3, 0
	DD	87549,	4096,	0.000172,	104352
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096_1, OFFSET xfft4096_2
	DPTR			OFFSET xfft4096_3, OFFSET xfft4096_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			512, 512
	DD			255, 127, 15*256+31, 15*256+63, 7*256+3, 0
	DD	108800,	5120,	0.000244,	135264
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft5120_1, OFFSET xfft5120_2
	DPTR			OFFSET xfft5120_3, OFFSET xfft5120_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			640, 512*2048+128
	DD			319, 159, 9*256+39, 19*256+79, 4*256+2, 0
	DD	129900,	6144,	0.000291,	162336
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft6144_1, OFFSET xfft6144_2
	DPTR			OFFSET xfft6144_3, OFFSET xfft6144_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			768, 768
	DD			383, 191, 11*256+47, 23*256+95, 5*256+2, 0
	DD	151300,	7168,	0.000381,	189504
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft7168_1, OFFSET xfft7168_2
	DPTR			OFFSET xfft7168_3, OFFSET xfft7168_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			896, (512*2048+256)*2048+128
	DD			447, 223, 13*256+55, 27*256+111, 6*256+3, 0
	DD	172700,	8192,	0.000395,	208800
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft8192_1, OFFSET xfft8192_2
	DPTR			OFFSET xfft8192_3, OFFSET xfft8192_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			1024, 1024
	DD			511, 255, 15*256+63, 31*256+127, 7*256+3, 0
	DD	214400,	10240,	0.000514,	50496
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft10K_1AMD, OFFSET xfft10K_2AMD
	DPTR			OFFSET xfft10K_3AMD, OFFSET xfft10K_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10, 1
	DD			9, 1, 8*1024+2, 1, 1, 0
	DD	214400,	10240,	0.000514,	50496
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft10K_1, OFFSET xfft10K_2
	DPTR			OFFSET xfft10K_3, OFFSET xfft10K_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10, 1
	DD			9, 1, 8*1024+2, 1, 1, 0
	DD	255300,	12288,	0.000626,	57216
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft12K_1AMD, OFFSET xfft12K_2AMD
	DPTR			OFFSET xfft12K_3AMD, OFFSET xfft12K_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12, 1
	DD			11, 1, 12, 1, 1, 0
	DD	255300,	12288,	0.000626,	57216
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft12K_1, OFFSET xfft12K_2
	DPTR			OFFSET xfft12K_3, OFFSET xfft12K_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12, 1
	DD			11, 1, 12, 1, 1, 0
	DD	297300,	14336,	0.000759,	64032
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft14K_1AMD, OFFSET xfft14K_2AMD
	DPTR			OFFSET xfft14K_3AMD, OFFSET xfft14K_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			14, 1
	DD			13, 1, (8*1024+4)*1024+2, 1, 1, 0
	DD	297300,	14336,	0.000759,	64032
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft14K_1, OFFSET xfft14K_2
	DPTR			OFFSET xfft14K_3, OFFSET xfft14K_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			14, 1
	DD			13, 1, (8*1024+4)*1024+2, 1, 1, 0
	DD	340400,	16384,	0.000857,	70464
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft16K_1AMD, OFFSET xfft16K_2AMD
	DPTR			OFFSET xfft16K_3AMD, OFFSET xfft16K_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			16, 1
	DD			15, 1, 16, 1, 1, 0
	DD	340400,	16384,	0.000857,	70464
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft16K_1, OFFSET xfft16K_2
	DPTR			OFFSET xfft16K_3, OFFSET xfft16K_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			16, 1
	DD			15, 1, 16, 1, 1, 0
	DD	423300,	20480,	0.00109,	83136
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft20K_1AMD, OFFSET xfft20K_2AMD
	DPTR			OFFSET xfft20K_3AMD, OFFSET xfft20K_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			20, 1
	DD			19, 1, 16*1024+4, 1, 1, 0
	DD	423300,	20480,	0.00109,	83136
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft20K_1, OFFSET xfft20K_2
	DPTR			OFFSET xfft20K_3, OFFSET xfft20K_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			20, 1
	DD			19, 1, 16*1024+4, 1, 1, 0
	DD	504600,	24576,	0.00135,	96384
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft24K_1AMD, OFFSET xfft24K_2AMD
	DPTR			OFFSET xfft24K_3AMD, OFFSET xfft24K_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			24, 1
	DD			23, 1, 24, 1, 1, 0
	DD	504600,	24576,	0.00135,	96384
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft24K_1, OFFSET xfft24K_2
	DPTR			OFFSET xfft24K_3, OFFSET xfft24K_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			24, 1
	DD			23, 1, 24, 1, 1, 0
	DD	587500,	28672,	0.00164,	109728
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft28K_1AMD, OFFSET xfft28K_2AMD
	DPTR			OFFSET xfft28K_3AMD, OFFSET xfft28K_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			28, 1
	DD			27, 1, (16*1024+8)*1024+4, 1, 1, 0
	DD	587500,	28672,	0.00164,	109728
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft28K_1, OFFSET xfft28K_2
	DPTR			OFFSET xfft28K_3, OFFSET xfft28K_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			28, 1
	DD			27, 1, (16*1024+8)*1024+4, 1, 1, 0
	DD	671400,	32768,	0.00179,	122688
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft32K_1AMD, OFFSET xfft32K_2AMD
	DPTR			OFFSET xfft32K_3AMD, OFFSET xfft32K_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			32, 1
	DD			31, 1, 32, 1, 1, 0
	DD	671400,	32768,	0.00179,	122688
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft32K_1, OFFSET xfft32K_2
	DPTR			OFFSET xfft32K_3, OFFSET xfft32K_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			32, 1
	DD			31, 1, 32, 1, 1, 0
	DD	835200,	40960,	0.00248,	244768
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft40K_1AMD, OFFSET xfft40K_2AMD
	DPTR			OFFSET xfft40K_3AMD, OFFSET xfft40K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			5, 1
	DD			4, 1, 4*1024+1, 1, 1, 0
	DD	835200,	40960,	0.00248,	244768
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft40K_1, OFFSET xfft40K_2
	DPTR			OFFSET xfft40K_3, OFFSET xfft40K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			5, 1
	DD			4, 1, 4*1024+1, 1, 1, 0
	DD	995500, 49152,	0.00303,	265952
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft48K_1AMD, OFFSET xfft48K_2AMD
	DPTR			OFFSET xfft48K_3AMD, OFFSET xfft48K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			6, 1
	DD			5, 1, 6, 1, 1, 0
	DD	995500, 49152,	0.00303,	265952
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft48K_1, OFFSET xfft48K_2
	DPTR			OFFSET xfft48K_3, OFFSET xfft48K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			6, 1
	DD			5, 1, 6, 1, 1, 0
	DD	1158000, 57344,	0.00363,	287232
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft56K_1AMD, OFFSET xfft56K_2AMD
	DPTR			OFFSET xfft56K_3AMD, OFFSET xfft56K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			7, 1
	DD			6, 1, (4*1024+2)*1024+1, 1, 1, 0
	DD	1158000, 57344,	0.00363,	287232
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft56K_1, OFFSET xfft56K_2
	DPTR			OFFSET xfft56K_3, OFFSET xfft56K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			7, 1
	DD			6, 1, (4*1024+2)*1024+1, 1, 1, 0
	DD	1325000, 65536,	0.00404,	308416
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft64K_1AMD, OFFSET xfft64K_2AMD
	DPTR			OFFSET xfft64K_3AMD, OFFSET xfft64K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			8, 1
	DD			7, 1, 8, 1, 1, 0
	DD	1325000, 65536,	0.00404,	308416
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft64K_1, OFFSET xfft64K_2
	DPTR			OFFSET xfft64K_3, OFFSET xfft64K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			8, 1
	DD			7, 1, 8, 1, 1, 0
	DD	1648000, 81920, 0.00533,	351552
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft80K_1AMD, OFFSET xfft80K_2AMD
	DPTR			OFFSET xfft80K_3AMD, OFFSET xfft80K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			10, 1
	DD			9, 1, 8*1024+2, 1, 1, 0
	DD	1648000, 81920, 0.00533,	351552
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft80K_1, OFFSET xfft80K_2
	DPTR			OFFSET xfft80K_3, OFFSET xfft80K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			10, 1
	DD			9, 1, 8*1024+2, 1, 1, 0
	DD	1966000, 98304, 0.00644,	394112
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft96K_1AMD, OFFSET xfft96K_2AMD
	DPTR			OFFSET xfft96K_3AMD, OFFSET xfft96K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			12, 1
	DD			11, 1, 12, 1, 1, 0
	DD	1966000, 98304, 0.00644,	394112
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft96K_1, OFFSET xfft96K_2
	DPTR			OFFSET xfft96K_3, OFFSET xfft96K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			12, 1
	DD			11, 1, 12, 1, 1, 0
	DD	2287000, 114688, 0.00771,	436768
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft112K_1AMD, OFFSET xfft112K_2AMD
	DPTR			OFFSET xfft112K_3AMD, OFFSET xfft112K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			14, 1
	DD			13, 1, (8*1024+4)*1024+2, 1, 1, 0
	DD	2287000, 114688, 0.00771,	436768
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft112K_1, OFFSET xfft112K_2
	DPTR			OFFSET xfft112K_3, OFFSET xfft112K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			14, 1
	DD			13, 1, (8*1024+4)*1024+2, 1, 1, 0
	DD	2614000, 131072, 0.00871,	479040
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft128K_1AMD, OFFSET xfft128K_2AMD
	DPTR			OFFSET xfft128K_3AMD, OFFSET xfft128K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			16, 1
	DD			15, 1, 16, 1, 1, 0
	DD	2614000, 131072, 0.00871,	479040
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft128K_1, OFFSET xfft128K_2
	DPTR			OFFSET xfft128K_3, OFFSET xfft128K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			16, 1
	DD			15, 1, 16, 1, 1, 0
	DD	3251000, 163840, 0.0111,	563392
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft160K_1AMD, OFFSET xfft160K_2AMD
	DPTR			OFFSET xfft160K_3AMD, OFFSET xfft160K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			20, 1
	DD			19, 1, 16*1024+4, 1, 1, 0
	DD	3251000, 163840, 0.0111,	563392
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft160K_1, OFFSET xfft160K_2
	DPTR			OFFSET xfft160K_3, OFFSET xfft160K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			20, 1
	DD			19, 1, 16*1024+4, 1, 1, 0
	DD	3875000, 196608, 0.0136,	648320
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft192K_1AMD, OFFSET xfft192K_2AMD
	DPTR			OFFSET xfft192K_3AMD, OFFSET xfft192K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			24, 1
	DD			23, 1, 24, 1, 1, 0
	DD	3875000, 196608, 0.0136,	648320
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft192K_1, OFFSET xfft192K_2
	DPTR			OFFSET xfft192K_3, OFFSET xfft192K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			24, 1
	DD			23, 1, 24, 1, 1, 0
	DD	4512000, 229376, 0.0160,	733344
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft224K_1AMD, OFFSET xfft224K_2AMD
	DPTR			OFFSET xfft224K_3AMD, OFFSET xfft224K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			28, 1
	DD			27, 1, (16*1024+8)*1024+4, 1, 1, 0
	DD	4512000, 229376, 0.0160,	733344
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft224K_1, OFFSET xfft224K_2
	DPTR			OFFSET xfft224K_3, OFFSET xfft224K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			28, 1
	DD			27, 1, (16*1024+8)*1024+4, 1, 1, 0
	DD	5158000, 262144, 0.0179,	817984
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft256K_1AMD, OFFSET xfft256K_2AMD
	DPTR			OFFSET xfft256K_3AMD, OFFSET xfft256K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			32, 1
	DD			31, 1, 32, 1, 1, 0
	DD	5158000, 262144, 0.0179,	817984
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft256K_1, OFFSET xfft256K_2
	DPTR			OFFSET xfft256K_3, OFFSET xfft256K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			32, 1
	DD			31, 1, 32, 1, 1, 0
	DD	6421000, 327680, 0.0239,	988896
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft320K_1AMD, OFFSET xfft320K_2AMD
	DPTR			OFFSET xfft320K_3AMD, OFFSET xfft320K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			40, 1
	DD			39, 1, 32*1024+8, 1, 1, 0
	DD	6421000, 327680, 0.0239,	988896
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft320K_1, OFFSET xfft320K_2
	DPTR			OFFSET xfft320K_3, OFFSET xfft320K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			40, 1
	DD			39, 1, 32*1024+8, 1, 1, 0
	DD	7651000, 393216, 0.0273,	1158944
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft384K_1AMD, OFFSET xfft384K_2AMD
	DPTR			OFFSET xfft384K_3AMD, OFFSET xfft384K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			48, 1
	DD			47, 1, 48, 1, 1, 0
	DD	7651000, 393216, 0.0273,	1158944
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft384K_1, OFFSET xfft384K_2
	DPTR			OFFSET xfft384K_3, OFFSET xfft384K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			48, 1
	DD			47, 1, 48, 1, 1, 0
	DD	8908000, 458752, 0.0323,	1329088
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft448K_1AMD, OFFSET xfft448K_2AMD
	DPTR			OFFSET xfft448K_3AMD, OFFSET xfft448K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			56, 1
	DD			55, 1, (32*1024+16)*1024+8, 1, 1, 0
	DD	8908000, 458752, 0.0323,	1329088
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft448K_1, OFFSET xfft448K_2
	DPTR			OFFSET xfft448K_3, OFFSET xfft448K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			56, 1
	DD			55, 1, (32*1024+16)*1024+8, 1, 1, 0
	DD	10180000, 524288, 0.0368,	1497408
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft512K_1AMD, OFFSET xfft512K_2AMD
	DPTR			OFFSET xfft512K_3AMD, OFFSET xfft512K_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			64, 1
	DD			63, 1, 64, 1, 1, 0
	DD	10180000, 524288, 0.0368,	1497408
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft512K_1, OFFSET xfft512K_2
	DPTR			OFFSET xfft512K_3, OFFSET xfft512K_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			64, 1
	DD			63, 1, 64, 1, 1, 0
	DD	12650000, 655360, 0.0464,	1838176+42240
	DD			RPFW+256*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			42240		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft640K4_1AMD, OFFSET xfft640K4_2AMD
	DPTR			OFFSET xfft640K4_3AMD, OFFSET xfft640K4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			80, 1
	DD			79, 1, 64*1024+16, 1, 1, 0
	DD	12650000, 655360, 0.0464,	1838176+42240
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			42240		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft640K4_1, OFFSET xfft640K4_2
	DPTR			OFFSET xfft640K4_3, OFFSET xfft640K4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			80, 1
	DD			79, 1, 64*1024+16, 1, 1, 0
	DD	12650000, 655360, 0.0464,	1838176+21760
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			21760		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft640K2_1AMD, OFFSET xfft640K2_2AMD
	DPTR			OFFSET xfft640K2_3AMD, OFFSET xfft640K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			80, 1
	DD			79, 1, 64*1024+16, 1, 1, 0
	DD	12650000, 655360, 0.0464,	1838176+21760
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			21760		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft640K2_1, OFFSET xfft640K2_2
	DPTR			OFFSET xfft640K2_3, OFFSET xfft640K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			80, 1
	DD			79, 1, 64*1024+16, 1, 1, 0
	DD	15070000, 786432, 0.0566,	2178080+50688
	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			50688		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft768K4_1AMD, OFFSET xfft768K4_2AMD
	DPTR			OFFSET xfft768K4_3AMD, OFFSET xfft768K4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			96, 1
	DD			95, 1, 96, 1, 1, 0
	DD	15070000, 786432, 0.0566,	2178080+50688
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			50688		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft768K4_1, OFFSET xfft768K4_2
	DPTR			OFFSET xfft768K4_3, OFFSET xfft768K4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			96, 1
	DD			95, 1, 96, 1, 1, 0
	DD	15070000, 786432, 0.0566,	2178080+26112
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			26112		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft768K2_1AMD, OFFSET xfft768K2_2AMD
	DPTR			OFFSET xfft768K2_3AMD, OFFSET xfft768K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			96, 1
	DD			95, 1, 96, 1, 1, 0
	DD	15070000, 786432, 0.0566,	2178080+26112
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			26112		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft768K2_1, OFFSET xfft768K2_2
	DPTR			OFFSET xfft768K2_3, OFFSET xfft768K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			96, 1
	DD			95, 1, 96, 1, 1, 0
	DD	17550000, 917504, 0.0693,	2518080+59136
	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			59136		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft896K4_1AMD, OFFSET xfft896K4_2AMD
	DPTR			OFFSET xfft896K4_3AMD, OFFSET xfft896K4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, (64*1024+32)*1024+16, 1, 1, 0
	DD	17550000, 917504, 0.0693,	2518080+59136
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			59136		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft896K4_1, OFFSET xfft896K4_2
	DPTR			OFFSET xfft896K4_3, OFFSET xfft896K4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, (64*1024+32)*1024+16, 1, 1, 0
	DD	17550000, 917504, 0.0693,	2518080+30464
	DD			RPFW+256*655362*2 ;; Flags, min_l2_cache, clm
	DD			30464		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft896K2_1AMD, OFFSET xfft896K2_2AMD
	DPTR			OFFSET xfft896K2_3AMD, OFFSET xfft896K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, (64*1024+32)*1024+16, 1, 1, 0
	DD	17550000, 917504, 0.0693,	2518080+30464
	DD			256*655362*2	;; Flags, min_l2_cache, clm
	DD			30464		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft896K2_1, OFFSET xfft896K2_2
	DPTR			OFFSET xfft896K2_3, OFFSET xfft896K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, (64*1024+32)*1024+16, 1, 1, 0
	DD	17550000, 917504, 0.0693,	2518080+16128
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			16128		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft896K1_1AMD, OFFSET xfft896K1_2AMD
	DPTR			OFFSET xfft896K1_3AMD, OFFSET xfft896K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, (64*1024+32)*1024+16, 1, 1, 0
	DD	17550000, 917504, 0.0693,	2518080+16128
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			16128		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft896K1_1, OFFSET xfft896K1_2
	DPTR			OFFSET xfft896K1_3, OFFSET xfft896K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			112, 1
	DD			111, 1, (64*1024+32)*1024+16, 1, 1, 0
	DD	20050000, 1048576, 0.0761,	2856256+67584
	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			67584		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024K4_1AMD, OFFSET xfft1024K4_2AMD
	DPTR			OFFSET xfft1024K4_3AMD, OFFSET xfft1024K4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 128, 1, 1, 0
	DD	20050000, 1048576, 0.0761,	2856256+67584
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			67584		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024K4_1, OFFSET xfft1024K4_2
	DPTR			OFFSET xfft1024K4_3, OFFSET xfft1024K4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 128, 1, 1, 0
	DD	20050000, 1048576, 0.0761,	2856256+34816
	DD			RPFW+256*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			34816		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024K2_1AMD, OFFSET xfft1024K2_2AMD
	DPTR			OFFSET xfft1024K2_3AMD, OFFSET xfft1024K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 128, 1, 1, 0
	DD	20050000, 1048576, 0.0761,	2856256+34816
	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
	DD			34816		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024K2_1, OFFSET xfft1024K2_2
	DPTR			OFFSET xfft1024K2_3, OFFSET xfft1024K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 128, 1, 1, 0
	DD	20050000, 1048576, 0.0761,	2856256+18432
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			18432		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024K1_1AMD, OFFSET xfft1024K1_2AMD
	DPTR			OFFSET xfft1024K1_3AMD, OFFSET xfft1024K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 128, 1, 1, 0
	DD	20050000, 1048576, 0.0761,	2856256+18432
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			18432		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024K1_1, OFFSET xfft1024K1_2
	DPTR			OFFSET xfft1024K1_3, OFFSET xfft1024K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			128, 1
	DD			127, 1, 128, 1, 1, 0
	DD	24930000, 1310720, 0.102,	3538560+84480
	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			84480		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1280K4_1AMD, OFFSET xfft1280K4_2AMD
	DPTR			OFFSET xfft1280K4_3AMD, OFFSET xfft1280K4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 128*1024+32, 1, 1, 0
	DD	24930000, 1310720, 0.102,	3538560+84480
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			84480		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1280K4_1, OFFSET xfft1280K4_2
	DPTR			OFFSET xfft1280K4_3, OFFSET xfft1280K4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 128*1024+32, 1, 1, 0
	DD	24930000, 1310720, 0.102,	3538560+43520
	DD			RPFW+256*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			43520		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1280K2_1AMD, OFFSET xfft1280K2_2AMD
	DPTR			OFFSET xfft1280K2_3AMD, OFFSET xfft1280K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 128*1024+32, 1, 1, 0
	DD	24930000, 1310720, 0.102,	3538560+43520
	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
	DD			43520		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1280K2_1, OFFSET xfft1280K2_2
	DPTR			OFFSET xfft1280K2_3, OFFSET xfft1280K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 128*1024+32, 1, 1, 0
	DD	24930000, 1310720, 0.102,	3538560+23040
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			23040		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1280K1_1AMD, OFFSET xfft1280K1_2AMD
	DPTR			OFFSET xfft1280K1_3AMD, OFFSET xfft1280K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 128*1024+32, 1, 1, 0
	DD	24930000, 1310720, 0.102,	3538560+23040
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			23040		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1280K1_1, OFFSET xfft1280K1_2
	DPTR			OFFSET xfft1280K1_3, OFFSET xfft1280K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			160, 1
	DD			159, 1, 128*1024+32, 1, 1, 0
	DD	29690000, 1572864, 0.125,	4218560+101376
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			101376		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536K4_1AMD, OFFSET xfft1536K4_2AMD
	DPTR			OFFSET xfft1536K4_3AMD, OFFSET xfft1536K4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 192, 1, 1, 0
	DD	29690000, 1572864, 0.125,	4218560+101376
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			101376		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536K4_1, OFFSET xfft1536K4_2
	DPTR			OFFSET xfft1536K4_3, OFFSET xfft1536K4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 192, 1, 1, 0
	DD	29690000, 1572864, 0.125,	4218560+52224
	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			52224		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536K2_1AMD, OFFSET xfft1536K2_2AMD
	DPTR			OFFSET xfft1536K2_3AMD, OFFSET xfft1536K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 192, 1, 1, 0
	DD	29690000, 1572864, 0.125,	4218560+52224
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			52224		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536K2_1, OFFSET xfft1536K2_2
	DPTR			OFFSET xfft1536K2_3, OFFSET xfft1536K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 192, 1, 1, 0
	DD	29690000, 1572864, 0.125,	4218560+27648
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			27648		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536K1_1AMD, OFFSET xfft1536K1_2AMD
	DPTR			OFFSET xfft1536K1_3AMD, OFFSET xfft1536K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 192, 1, 1, 0
	DD	29690000, 1572864, 0.125,	4218560+27648
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			27648		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536K1_1, OFFSET xfft1536K1_2
	DPTR			OFFSET xfft1536K1_3, OFFSET xfft1536K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			192, 1
	DD			191, 1, 192, 1, 1, 0
	DD	34560000, 1835008, 0.151,	4898656+118272
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			118272		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1792K4_1AMD, OFFSET xfft1792K4_2AMD
	DPTR			OFFSET xfft1792K4_3AMD, OFFSET xfft1792K4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, (128*1024+64)*1024+32, 1, 1, 0
	DD	34560000, 1835008, 0.151,	4898656+118272
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			118272		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1792K4_1, OFFSET xfft1792K4_2
	DPTR			OFFSET xfft1792K4_3, OFFSET xfft1792K4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, (128*1024+64)*1024+32, 1, 1, 0
	DD	34560000, 1835008, 0.151,	4898656+60928
	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			60928		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1792K2_1AMD, OFFSET xfft1792K2_2AMD
	DPTR			OFFSET xfft1792K2_3AMD, OFFSET xfft1792K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, (128*1024+64)*1024+32, 1, 1, 0
	DD	34560000, 1835008, 0.151,	4898656+60928
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			60928		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1792K2_1, OFFSET xfft1792K2_2
	DPTR			OFFSET xfft1792K2_3, OFFSET xfft1792K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, (128*1024+64)*1024+32, 1, 1, 0
	DD	34560000, 1835008, 0.151,	4898656+32256
	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			32256		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1792K1_1AMD, OFFSET xfft1792K1_2AMD
	DPTR			OFFSET xfft1792K1_3AMD, OFFSET xfft1792K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, (128*1024+64)*1024+32, 1, 1, 0
	DD	34560000, 1835008, 0.151,	4898656+32256
	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
	DD			32256		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1792K1_1, OFFSET xfft1792K1_2
	DPTR			OFFSET xfft1792K1_3, OFFSET xfft1792K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, (128*1024+64)*1024+32, 1, 1, 0
	DD	34560000, 1835008, 0.151,	4898656+32256
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			32256		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1792K0_1AMD, OFFSET xfft1792K0_2AMD
	DPTR			OFFSET xfft1792K0_3AMD, OFFSET xfft1792K0_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, (128*1024+64)*1024+32, 1, 1, 0
	DD	34560000, 1835008, 0.151,	4898656+32256
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			32256		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1792K0_1, OFFSET xfft1792K0_2
	DPTR			OFFSET xfft1792K0_3, OFFSET xfft1792K0_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			224, 1
	DD			223, 1, (128*1024+64)*1024+32, 1, 1, 0
	DD	39500000, 2097152, 0.169,	5573952+135168
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			135168		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048K4_1AMD, OFFSET xfft2048K4_2AMD
	DPTR			OFFSET xfft2048K4_3AMD, OFFSET xfft2048K4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 256, 1, 1, 0
	DD	39500000, 2097152, 0.169,	5573952+135168
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			135168		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048K4_1, OFFSET xfft2048K4_2
	DPTR			OFFSET xfft2048K4_3, OFFSET xfft2048K4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 256, 1, 1, 0
	DD	39500000, 2097152, 0.169,	5573952+69632
	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			69632		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048K2_1AMD, OFFSET xfft2048K2_2AMD
	DPTR			OFFSET xfft2048K2_3AMD, OFFSET xfft2048K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 256, 1, 1, 0
	DD	39500000, 2097152, 0.169,	5573952+69632
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			69632		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048K2_1, OFFSET xfft2048K2_2
	DPTR			OFFSET xfft2048K2_3, OFFSET xfft2048K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 256, 1, 1, 0
	DD	39500000, 2097152, 0.169,	5573952+36864
	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			36864		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048K1_1AMD, OFFSET xfft2048K1_2AMD
	DPTR			OFFSET xfft2048K1_3AMD, OFFSET xfft2048K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 256, 1, 1, 0
	DD	39500000, 2097152, 0.169,	5573952+36864
	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
	DD			36864		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048K1_1, OFFSET xfft2048K1_2
	DPTR			OFFSET xfft2048K1_3, OFFSET xfft2048K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 256, 1, 1, 0
	DD	39500000, 2097152, 0.169,	5573952+36864
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			36864		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048K0_1AMD, OFFSET xfft2048K0_2AMD
	DPTR			OFFSET xfft2048K0_3AMD, OFFSET xfft2048K0_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 256, 1, 1, 0
	DD	39500000, 2097152, 0.169,	5573952+36864
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			36864		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048K0_1, OFFSET xfft2048K0_2
	DPTR			OFFSET xfft2048K0_3, OFFSET xfft2048K0_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			256, 1
	DD			255, 1, 256, 1, 1, 0
	DD	49100000, 2621440, 0.222,	6937600+87040
	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			87040		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2560K2_1AMD, OFFSET xfft2560K2_2AMD
	DPTR			OFFSET xfft2560K2_3AMD, OFFSET xfft2560K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 256*1024+64, 1, 1, 0
	DD	49100000, 2621440, 0.222,	6937600+87040
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			87040		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2560K2_1, OFFSET xfft2560K2_2
	DPTR			OFFSET xfft2560K2_3, OFFSET xfft2560K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 256*1024+64, 1, 1, 0
	DD	49100000, 2621440, 0.222,	6937600+46080
	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			46080		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2560K1_1AMD, OFFSET xfft2560K1_2AMD
	DPTR			OFFSET xfft2560K1_3AMD, OFFSET xfft2560K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 256*1024+64, 1, 1, 0
	DD	49100000, 2621440, 0.222,	6937600+46080
	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
	DD			46080		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2560K1_1, OFFSET xfft2560K1_2
	DPTR			OFFSET xfft2560K1_3, OFFSET xfft2560K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 256*1024+64, 1, 1, 0
	DD	49100000, 2621440, 0.222,	6937600+46080
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			46080		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2560K0_1AMD, OFFSET xfft2560K0_2AMD
	DPTR			OFFSET xfft2560K0_3AMD, OFFSET xfft2560K0_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 256*1024+64, 1, 1, 0
	DD	49100000, 2621440, 0.222,	6937600+46080
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			46080		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2560K0_1, OFFSET xfft2560K0_2
	DPTR			OFFSET xfft2560K0_3, OFFSET xfft2560K0_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			320, 1
	DD			319, 1, 256*1024+64, 1, 1, 0
	DD	58520000, 3145728, 0.289,	8297408+104448
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			104448		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072K2_1AMD, OFFSET xfft3072K2_2AMD
	DPTR			OFFSET xfft3072K2_3AMD, OFFSET xfft3072K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 384, 1, 1, 0
	DD	58520000, 3145728, 0.289,	8297408+104448
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			104448		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072K2_1, OFFSET xfft3072K2_2
	DPTR			OFFSET xfft3072K2_3, OFFSET xfft3072K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 384, 1, 1, 0
	DD	58520000, 3145728, 0.289,	8297408+55296
	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			55296		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072K1_1AMD, OFFSET xfft3072K1_2AMD
	DPTR			OFFSET xfft3072K1_3AMD, OFFSET xfft3072K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 384, 1, 1, 0
	DD	58520000, 3145728, 0.289,	8297408+55296
	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
	DD			55296		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072K1_1, OFFSET xfft3072K1_2
	DPTR			OFFSET xfft3072K1_3, OFFSET xfft3072K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 384, 1, 1, 0
	DD	58520000, 3145728, 0.289,	8297408+55296
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			55296		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072K0_1AMD, OFFSET xfft3072K0_2AMD
	DPTR			OFFSET xfft3072K0_3AMD, OFFSET xfft3072K0_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 384, 1, 1, 0
	DD	58520000, 3145728, 0.289,	8297408+55296
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			55296		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072K0_1, OFFSET xfft3072K0_2
	DPTR			OFFSET xfft3072K0_3, OFFSET xfft3072K0_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			384, 1
	DD			383, 1, 384, 1, 1, 0
	DD	68130000, 3670016, 0.369,	9657312+121856
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			121856		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3584K2_1AMD, OFFSET xfft3584K2_2AMD
	DPTR			OFFSET xfft3584K2_3AMD, OFFSET xfft3584K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, (256*1024+128)*1024+64, 1, 1,0
	DD	68130000, 3670016, 0.369,	9657312+121856
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			121856		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3584K2_1, OFFSET xfft3584K2_2
	DPTR			OFFSET xfft3584K2_3, OFFSET xfft3584K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, (256*1024+128)*1024+64, 1, 1,0
	DD	68130000, 3670016, 0.369,	9657312+64512
	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			64512		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3584K1_1AMD, OFFSET xfft3584K1_2AMD
	DPTR			OFFSET xfft3584K1_3AMD, OFFSET xfft3584K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, (256*1024+128)*1024+64, 1, 1,0
	DD	68130000, 3670016, 0.369,	9657312+64512
	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
	DD			64512		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3584K1_1, OFFSET xfft3584K1_2
	DPTR			OFFSET xfft3584K1_3, OFFSET xfft3584K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, (256*1024+128)*1024+64, 1, 1,0
	DD	68130000, 3670016, 0.369,	9657312+64512
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			64512		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3584K0_1AMD, OFFSET xfft3584K0_2AMD
	DPTR			OFFSET xfft3584K0_3AMD, OFFSET xfft3584K0_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, (256*1024+128)*1024+64, 1, 1,0
	DD	68130000, 3670016, 0.369,	9657312+64512
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			64512		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3584K0_1, OFFSET xfft3584K0_2
	DPTR			OFFSET xfft3584K0_3, OFFSET xfft3584K0_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			448, 1
	DD			447, 1, (256*1024+128)*1024+64, 1, 1,0
	DD	77910000, 4194304, 0.425,	11009344+139264
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			139264		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096K2_1AMD, OFFSET xfft4096K2_2AMD
	DPTR			OFFSET xfft4096K2_3AMD, OFFSET xfft4096K2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 512, 1, 1, 0
	DD	77910000, 4194304, 0.425,	11009344+139264
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			139264		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096K2_1, OFFSET xfft4096K2_2
	DPTR			OFFSET xfft4096K2_3, OFFSET xfft4096K2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 512, 1, 1, 0
	DD	77910000, 4194304, 0.425,	11009344+73728
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			73728		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096K1_1AMD, OFFSET xfft4096K1_2AMD
	DPTR			OFFSET xfft4096K1_3AMD, OFFSET xfft4096K1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 512, 1, 1, 0
	DD	77910000, 4194304, 0.425,	11009344+73728
	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
	DD			73728		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096K1_1, OFFSET xfft4096K1_2
	DPTR			OFFSET xfft4096K1_3, OFFSET xfft4096K1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 512, 1, 1, 0
	DD	77910000, 4194304, 0.425,	11009344+73728
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			73728		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096K0_1AMD, OFFSET xfft4096K0_2AMD
	DPTR			OFFSET xfft4096K0_3AMD, OFFSET xfft4096K0_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 512, 1, 1, 0
	DD	0
	DD	77910000, 4194304, 0.425,	11009344+73728
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			73728		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096K0_1, OFFSET xfft4096K0_2
	DPTR			OFFSET xfft4096K0_3, OFFSET xfft4096K0_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			512, 1
	DD			511, 1, 512, 1, 1, 0
	DD	0
xjmptablep DD	739,	32,	0.00000111,	976
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft32p_1, OFFSET xfft32p_2
	DPTR			OFFSET xfft32p_3, OFFSET xfft32p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			4, 4
	DD			2, 1, 1, 1, 1, 0
	DD	1095,	48,	0.00000144,	1752
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft48p_1, OFFSET xfft48p_2
	DPTR			OFFSET xfft48p_3, OFFSET xfft48p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			6, 6
	DD			3, 1, 1, 1, 1, 0
	DD	1465,	64,	0.00000178,	2336
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft64p_1, OFFSET xfft64p_2
	DPTR			OFFSET xfft64p_3, OFFSET xfft64p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			8, 8
	DD			4, 1, 1, 1, 1, 0
	DD	2173,	96,	0.00000259,	3216
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft96p_1, OFFSET xfft96p_2
	DPTR			OFFSET xfft96p_3, OFFSET xfft96p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			12, 12
	DD			6, 3, 1, 1, 1, 0
	DD	2897,	128,	0.00000319,	4288
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft128p_1, OFFSET xfft128p_2
	DPTR			OFFSET xfft128p_3, OFFSET xfft128p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			16, 16
	DD			8, 4, 1, 2, 1, 0
	DD	4295,	192,	0.00000542,	6720
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft192p_1, OFFSET xfft192p_2
	DPTR			OFFSET xfft192p_3, OFFSET xfft192p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			24, 24
	DD			12, 6, 3, 3, 1, 0
	DD	5729,	256,	0.00000691,	8576
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft256p_1, OFFSET xfft256p_2
	DPTR			OFFSET xfft256p_3, OFFSET xfft256p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			32, 32
	DD			16, 8, 4, 4, 1, 0
	DD	8493,	384,	0.0000111,	13152
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft384p_1, OFFSET xfft384p_2
	DPTR			OFFSET xfft384p_3, OFFSET xfft384p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			48, 48
	DD			24, 12, 3, 6, 1, 0
	DD	11319,	512,	0.0000143,	17152
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft512p_1, OFFSET xfft512p_2
	DPTR			OFFSET xfft512p_3, OFFSET xfft512p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			64, 64
	DD			32, 16, 4, 2*256+8, 1, 0
	DD	16779,	768,	0.0000260,	26592
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft768p_1, OFFSET xfft768p_2
	DPTR			OFFSET xfft768p_3, OFFSET xfft768p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			96, 96
	DD			48, 24, 3*256+6, 3*256+12, 1, 0
	DD	22381,	1024,	0.0000349,	34304
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024p_1, OFFSET xfft1024p_2
	DPTR			OFFSET xfft1024p_3, OFFSET xfft1024p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			128, 128
	DD			64, 32, 4*256+8, 4*256+16, 1, 0
	DD	33189,	1536,	0.0000601,	52896
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536p_1, OFFSET xfft1536p_2
	DPTR			OFFSET xfft1536p_3, OFFSET xfft1536p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			192, 192
	DD			96, 48, 12, 6*256+24, 3, 0
	DD	44221,	2048,	0.0000773,	68608
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048p_1, OFFSET xfft2048p_2
	DPTR			OFFSET xfft2048p_3, OFFSET xfft2048p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			256, 256
	DD			128, 64, 16, 8*256+32, 2*256+4, 0
	DD	65519,	3072,	0.000131,	106080
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072p_1, OFFSET xfft3072p_2
	DPTR			OFFSET xfft3072p_3, OFFSET xfft3072p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			384, 384
	DD			192, 96, 6*256+24, 12*256+48, 3*256+3, 0
	DD	87271,	4096,	0.000172,	137216
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096p_1, OFFSET xfft4096p_2
	DPTR			OFFSET xfft4096p_3, OFFSET xfft4096p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			512, 512
	DD			256, 128, 8*256+32, 16*256+64, 4*256+4, 0
	DD	129600,	6144,	0.000291,	211872
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft6144p_1, OFFSET xfft6144p_2
	DPTR			OFFSET xfft6144p_3, OFFSET xfft6144p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			768, 768
	DD			384, 192, 12*256+48, 24*256+96, 6*256+3, 0
	DD	172400,	8192,	0.000395,	274432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			0		;; FFT levels done in pass2
	DPTR			OFFSET xfft8192p_1, OFFSET xfft8192p_2
	DPTR			OFFSET xfft8192p_3, OFFSET xfft8192p_4
	DPTR			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			1024, 1024
	DD			512, 256, 16*256+64, 32*256+128, 8*256+4, 0
	DD	254800,	12288,	0.000626,	51040
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft12Kp_1AMD, OFFSET xfft12Kp_2AMD
	DPTR			OFFSET xfft12Kp_3AMD, OFFSET xfft12Kp_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12, 1
	DD			12, 1, 12, 1, 1, 0
	DD	254800,	12288,	0.000626,	51040
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft12Kp_1, OFFSET xfft12Kp_2
	DPTR			OFFSET xfft12Kp_3, OFFSET xfft12Kp_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12, 1
	DD			12, 1, 12, 1, 1, 0
	DD	339100,	16384,	0.000857,	64256
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft16Kp_1AMD, OFFSET xfft16Kp_2AMD
	DPTR			OFFSET xfft16Kp_3AMD, OFFSET xfft16Kp_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			16, 1
	DD			16, 1, 16, 1, 1, 0
	DD	339100,	16384,	0.000857,	64256
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft16Kp_1, OFFSET xfft16Kp_2
	DPTR			OFFSET xfft16Kp_3, OFFSET xfft16Kp_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			16, 1
	DD			16, 1, 16, 1, 1, 0
	DD	503400,	24576,	0.00135,	91552
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft24Kp_1AMD, OFFSET xfft24Kp_2AMD
	DPTR			OFFSET xfft24Kp_3AMD, OFFSET xfft24Kp_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			24, 1
	DD			24, 1, 24, 1, 1, 0
	DD	503400,	24576,	0.00135,	91552
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft24Kp_1, OFFSET xfft24Kp_2
	DPTR			OFFSET xfft24Kp_3, OFFSET xfft24Kp_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			24, 1
	DD			24, 1, 24, 1, 1, 0
	DD	669600,	32768,	0.00179,	118272
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft32Kp_1AMD, OFFSET xfft32Kp_2AMD
	DPTR			OFFSET xfft32Kp_3AMD, OFFSET xfft32Kp_4AMD
	DPTR			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			32, 1
	DD			32, 1, 32, 1, 1, 0
	DD	669600,	32768,	0.00179,	118272
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			8		;; FFT levels done in pass2
	DPTR			OFFSET xfft32Kp_1, OFFSET xfft32Kp_2
	DPTR			OFFSET xfft32Kp_3, OFFSET xfft32Kp_4
	DPTR			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			32, 1
	DD			32, 1, 32, 1, 1, 0
	DD	992800, 49152,	0.00303,	209696
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft48Kp_1AMD, OFFSET xfft48Kp_2AMD
	DPTR			OFFSET xfft48Kp_3AMD, OFFSET xfft48Kp_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			6, 1
	DD			6, 1, 6, 1, 1, 0
	DD	992800, 49152,	0.00303,	209696
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft48Kp_1, OFFSET xfft48Kp_2
	DPTR			OFFSET xfft48Kp_3, OFFSET xfft48Kp_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			6, 1
	DD			6, 1, 6, 1, 1, 0
	DD	1320000, 65536,	0.00404,	252288
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft64Kp_1AMD, OFFSET xfft64Kp_2AMD
	DPTR			OFFSET xfft64Kp_3AMD, OFFSET xfft64Kp_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			8, 1
	DD			8, 1, 8, 1, 1, 0
	DD	1320000, 65536,	0.00404,	252288
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft64Kp_1, OFFSET xfft64Kp_2
	DPTR			OFFSET xfft64Kp_3, OFFSET xfft64Kp_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			8, 1
	DD			8, 1, 8, 1, 1, 0
	DD	1962000, 98304, 0.00644,	337760
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft96Kp_1AMD, OFFSET xfft96Kp_2AMD
	DPTR			OFFSET xfft96Kp_3AMD, OFFSET xfft96Kp_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			12, 1
	DD			12, 1, 12, 1, 1, 0
	DD	1962000, 98304, 0.00644,	337760
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft96Kp_1, OFFSET xfft96Kp_2
	DPTR			OFFSET xfft96Kp_3, OFFSET xfft96Kp_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			12, 1
	DD			12, 1, 12, 1, 1, 0
	DD	2610000, 131072, 0.00871,	422656
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft128Kp_1AMD, OFFSET xfft128Kp_2AMD
	DPTR			OFFSET xfft128Kp_3AMD, OFFSET xfft128Kp_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			16, 1
	DD			16, 1, 16, 1, 1, 0
	DD	2610000, 131072, 0.00871,	422656
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft128Kp_1, OFFSET xfft128Kp_2
	DPTR			OFFSET xfft128Kp_3, OFFSET xfft128Kp_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			16, 1
	DD			16, 1, 16, 1, 1, 0
	DD	3867000, 196608, 0.0136,	593312
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft192Kp_1AMD, OFFSET xfft192Kp_2AMD
	DPTR			OFFSET xfft192Kp_3AMD, OFFSET xfft192Kp_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			24, 1
	DD			24, 1, 24, 1, 1, 0
	DD	3867000, 196608, 0.0136,	593312
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft192Kp_1, OFFSET xfft192Kp_2
	DPTR			OFFSET xfft192Kp_3, OFFSET xfft192Kp_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			24, 1
	DD			24, 1, 24, 1, 1, 0
	DD	5146000, 262144, 0.0179,	763392
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft256Kp_1AMD, OFFSET xfft256Kp_2AMD
	DPTR			OFFSET xfft256Kp_3AMD, OFFSET xfft256Kp_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			32, 1
	DD			32, 1, 32, 1, 1, 0
	DD	5146000, 262144, 0.0179,	763392
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft256Kp_1, OFFSET xfft256Kp_2
	DPTR			OFFSET xfft256Kp_3, OFFSET xfft256Kp_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			32, 1
	DD			32, 1, 32, 1, 1, 0
	DD	7635000, 393216, 0.0273,	1104992
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft384Kp_1AMD, OFFSET xfft384Kp_2AMD
	DPTR			OFFSET xfft384Kp_3AMD, OFFSET xfft384Kp_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			48, 1
	DD			48, 1, 48, 1, 1, 0
	DD	7635000, 393216, 0.0273,	1104992
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft384Kp_1, OFFSET xfft384Kp_2
	DPTR			OFFSET xfft384Kp_3, OFFSET xfft384Kp_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			48, 1
	DD			48, 1, 48, 1, 1, 0
	DD	10150000, 524288, 0.0368,	1444864
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft512Kp_1AMD, OFFSET xfft512Kp_2AMD
	DPTR			OFFSET xfft512Kp_3AMD, OFFSET xfft512Kp_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			64, 1
	DD			64, 1, 64, 1, 1, 0
	DD	10150000, 524288, 0.0368,	1444864
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft512Kp_1, OFFSET xfft512Kp_2
	DPTR			OFFSET xfft512Kp_3, OFFSET xfft512Kp_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			64, 1
	DD			64, 1, 64, 1, 1, 0
	DD	15040000, 786432, 0.0566,	2178464
	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			50688		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft768Kp4_1AMD, OFFSET xfft768Kp4_2AMD
	DPTR			OFFSET xfft768Kp4_3AMD, OFFSET xfft768Kp4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			96, 1
	DD			96, 1, 96, 1, 1, 0
	DD	15040000, 786432, 0.0566,	2178464
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			50688		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft768Kp4_1, OFFSET xfft768Kp4_2
	DPTR			OFFSET xfft768Kp4_3, OFFSET xfft768Kp4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			96, 1
	DD			96, 1, 96, 1, 1, 0
	DD	15040000, 786432, 0.0566,	2153888
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			26112		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft768Kp2_1AMD, OFFSET xfft768Kp2_2AMD
	DPTR			OFFSET xfft768Kp2_3AMD, OFFSET xfft768Kp2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			96, 1
	DD			96, 1, 96, 1, 1, 0
	DD	15040000, 786432, 0.0566,	2153888
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			26112		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft768Kp2_1, OFFSET xfft768Kp2_2
	DPTR			OFFSET xfft768Kp2_3, OFFSET xfft768Kp2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			96, 1
	DD			96, 1, 96, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	2875392
	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			67584		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024Kp4_1AMD,OFFSET xfft1024Kp4_2AMD
	DPTR			OFFSET xfft1024Kp4_3AMD,OFFSET xfft1024Kp4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 128, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	2875392
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			67584		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024Kp4_1, OFFSET xfft1024Kp4_2
	DPTR			OFFSET xfft1024Kp4_3, OFFSET xfft1024Kp4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 128, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	2842624
	DD			RPFW+256*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			34816		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024Kp2_1AMD,OFFSET xfft1024Kp2_2AMD
	DPTR			OFFSET xfft1024Kp2_3AMD,OFFSET xfft1024Kp2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 128, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	2842624
	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
	DD			34816		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024Kp2_1, OFFSET xfft1024Kp2_2
	DPTR			OFFSET xfft1024Kp2_3, OFFSET xfft1024Kp2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 128, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	2826240
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			18432		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024Kp1_1AMD,OFFSET xfft1024Kp1_2AMD
	DPTR			OFFSET xfft1024Kp1_3AMD,OFFSET xfft1024Kp1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 128, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	2826240
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			18432		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1024Kp1_1, OFFSET xfft1024Kp1_2
	DPTR			OFFSET xfft1024Kp1_3, OFFSET xfft1024Kp1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			128, 1
	DD			128, 1, 128, 1, 1, 0
	DD	29640000, 1572864, 0.125,	4275296
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			101376		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536Kp4_1AMD,OFFSET xfft1536Kp4_2AMD
	DPTR			OFFSET xfft1536Kp4_3AMD,OFFSET xfft1536Kp4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 192, 1, 1, 0
	DD	29640000, 1572864, 0.125,	4275296
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			101376		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536Kp4_1, OFFSET xfft1536Kp4_2
	DPTR			OFFSET xfft1536Kp4_3, OFFSET xfft1536Kp4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 192, 1, 1, 0
	DD	29640000, 1572864, 0.125,	4226144
	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			52224		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536Kp2_1AMD,OFFSET xfft1536Kp2_2AMD
	DPTR			OFFSET xfft1536Kp2_3AMD,OFFSET xfft1536Kp2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 192, 1, 1, 0
	DD	29640000, 1572864, 0.125,	4226144
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			52224		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536Kp2_1, OFFSET xfft1536Kp2_2
	DPTR			OFFSET xfft1536Kp2_3, OFFSET xfft1536Kp2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 192, 1, 1, 0
	DD	29640000, 1572864, 0.125,	4201568
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			27648		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536Kp1_1AMD,OFFSET xfft1536Kp1_2AMD
	DPTR			OFFSET xfft1536Kp1_3AMD,OFFSET xfft1536Kp1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 192, 1, 1, 0
	DD	29640000, 1572864, 0.125,	4201568
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			27648		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft1536Kp1_1, OFFSET xfft1536Kp1_2
	DPTR			OFFSET xfft1536Kp1_3, OFFSET xfft1536Kp1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			192, 1
	DD			192, 1, 192, 1, 1, 0
	DD	39390000, 2097152, 0.169,	5668864
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			135168		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048Kp4_1AMD,OFFSET xfft2048Kp4_2AMD
	DPTR			OFFSET xfft2048Kp4_3AMD,OFFSET xfft2048Kp4_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 256, 1, 1, 0
	DD	39390000, 2097152, 0.169,	5668864
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			135168		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048Kp4_1, OFFSET xfft2048Kp4_2
	DPTR			OFFSET xfft2048Kp4_3, OFFSET xfft2048Kp4_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 256, 1, 1, 0
	DD	39390000, 2097152, 0.169,	5603328
	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			69632		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048Kp2_1AMD,OFFSET xfft2048Kp2_2AMD
	DPTR			OFFSET xfft2048Kp2_3AMD,OFFSET xfft2048Kp2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 256, 1, 1, 0
	DD	39390000, 2097152, 0.169,	5603328
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			69632		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048Kp2_1, OFFSET xfft2048Kp2_2
	DPTR			OFFSET xfft2048Kp2_3, OFFSET xfft2048Kp2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 256, 1, 1, 0
	DD	39390000, 2097152, 0.169,	5570560
	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			36864		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048Kp1_1AMD,OFFSET xfft2048Kp1_2AMD
	DPTR			OFFSET xfft2048Kp1_3AMD,OFFSET xfft2048Kp1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 256, 1, 1, 0
	DD	39390000, 2097152, 0.169,	5570560
	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
	DD			36864		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048Kp1_1, OFFSET xfft2048Kp1_2
	DPTR			OFFSET xfft2048Kp1_3, OFFSET xfft2048Kp1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 256, 1, 1, 0
	DD	39390000, 2097152, 0.169,	5570560
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			36864		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048Kp0_1AMD,OFFSET xfft2048Kp0_2AMD
	DPTR			OFFSET xfft2048Kp0_3AMD,OFFSET xfft2048Kp0_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 256, 1, 1, 0
	DD	39390000, 2097152, 0.169,	5570560
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			36864		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft2048Kp0_1, OFFSET xfft2048Kp0_2
	DPTR			OFFSET xfft2048Kp0_3, OFFSET xfft2048Kp0_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			256, 1
	DD			256, 1, 256, 1, 1, 0
	DD	58410000, 3145728, 0.289,	8370080
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			104448		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072Kp2_1AMD,OFFSET xfft3072Kp2_2AMD
	DPTR			OFFSET xfft3072Kp2_3AMD,OFFSET xfft3072Kp2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 384, 1, 1, 0
	DD	58410000, 3145728, 0.289,	8370080
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			104448		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072Kp2_1, OFFSET xfft3072Kp2_2
	DPTR			OFFSET xfft3072Kp2_3, OFFSET xfft3072Kp2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 384, 1, 1, 0
	DD	58410000, 3145728, 0.289,	8320928
	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			55296		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072Kp1_1AMD,OFFSET xfft3072Kp1_2AMD
	DPTR			OFFSET xfft3072Kp1_3AMD,OFFSET xfft3072Kp1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 384, 1, 1, 0
	DD	58410000, 3145728, 0.289,	8320928
	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
	DD			55296		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072Kp1_1, OFFSET xfft3072Kp1_2
	DPTR			OFFSET xfft3072Kp1_3, OFFSET xfft3072Kp1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 384, 1, 1, 0
	DD	58410000, 3145728, 0.289,	8320928
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			55296		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072Kp0_1AMD,OFFSET xfft3072Kp0_2AMD
	DPTR			OFFSET xfft3072Kp0_3AMD,OFFSET xfft3072Kp0_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 384, 1, 1, 0
	DD	58410000, 3145728, 0.289,	8320928
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			55296		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft3072Kp0_1, OFFSET xfft3072Kp0_2
	DPTR			OFFSET xfft3072Kp0_3, OFFSET xfft3072Kp0_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			384, 1
	DD			384, 1, 384, 1, 1, 0
	DD	77700000, 4194304, 0.425,	11124736
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			139264		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096Kp2_1AMD,OFFSET xfft4096Kp2_2AMD
	DPTR			OFFSET xfft4096Kp2_3AMD,OFFSET xfft4096Kp2_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			512, 1
	DD			512, 1, 512, 1, 1, 0
	DD	77700000, 4194304, 0.425,	11124736
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			139264		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096Kp2_1, OFFSET xfft4096Kp2_2
	DPTR			OFFSET xfft4096Kp2_3, OFFSET xfft4096Kp2_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			512, 1
	DD			512, 1, 512, 1, 1, 0
	DD	77700000, 4194304, 0.425,	11059200
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			73728		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096Kp1_1AMD,OFFSET xfft4096Kp1_2AMD
	DPTR			OFFSET xfft4096Kp1_3AMD,OFFSET xfft4096Kp1_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			512, 1
	DD			512, 1, 512, 1, 1, 0
	DD	77700000, 4194304, 0.425,	11059200
	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
	DD			73728		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096Kp1_1, OFFSET xfft4096Kp1_2
	DPTR			OFFSET xfft4096Kp1_3, OFFSET xfft4096Kp1_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			512, 1
	DD			512, 1, 512, 1, 1, 0
	DD	77700000, 4194304, 0.425,	11059200
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			73728		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096Kp0_1AMD,OFFSET xfft4096Kp0_2AMD
	DPTR			OFFSET xfft4096Kp0_3AMD,OFFSET xfft4096Kp0_4AMD
	DPTR			OFFSET xprctab3a ;; Table of add/sub/norm procs
	DD			512, 1
	DD			512, 1, 512, 1, 1, 0
	DD	77700000, 4194304, 0.425,	11059200
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			73728		;; scratch area size
	DD			11		;; FFT levels done in pass2
	DPTR			OFFSET xfft4096Kp0_1, OFFSET xfft4096Kp0_2
	DPTR			OFFSET xfft4096Kp0_3, OFFSET xfft4096Kp0_4
	DPTR			OFFSET xprctab3	;; Table of add/sub/norm procs
	DD			512, 1
	DD			512, 1, 512, 1, 1, 0
	DD	0

	;; Align so that other GWDATA areas are aligned on a cache line
	align 128
_GWDATA ENDS

; Macros to copy a double once, twice, or three times from rsi to rdi and rbx.
; C code has calculated the value and we must copy or copy and duplicate it in
; the x87 and XMM global variables.

copy_float_once MACRO dest
	mov	rdi, OFFSET dest
	movsd				; Copy float once
	ENDM

copy_double_once MACRO dest
	mov	rdi, OFFSET dest
	movsd				; Copy double once
	movsd
	ENDM

copy_double_twice MACRO dest
	mov	rdi, OFFSET dest
	movsd				; Copy double once
	movsd
	lea	rsi, [rsi-8]		; Copy double again
	movsd
	movsd
	ENDM

copy_double_thrice MACRO dest1, dest2
	mov	rdi, OFFSET dest1
	movsd				; Copy double once
	movsd
	lea	rsi, [rsi-8]
	mov	rdi, OFFSET dest2
	movsd				; Copy double once
	movsd
	lea	rsi, [rsi-8]		; Copy double again
	movsd
	movsd
	ENDM

copy_eight_doubles_four_times MACRO dest
	mov	rdi, OFFSET dest
	mov	rcx, 16			; Copy 8 doubles
	rep movsd
	lea	rsi, [rsi-8*8]		; Copy again
	mov	rcx, 16
	rep movsd
	lea	rsi, [rsi-8*8]		; Copy again
	mov	rcx, 16
	rep movsd
	lea	rsi, [rsi-8*8]		; Copy again
	mov	rcx, 16
	rep movsd
	ENDM

_TEXT	SEGMENT

; Return address of jmp tables for C code to examine

	PUBLIC	_gwinfo1
_gwinfo1 PROC
	mov	rax, OFFSET xjmptable	 ; P4 mersenne mod FFTs
	mov	_INFT+0*SZPTR, rax
	mov	rax, OFFSET xjmptablep	 ; P4 2^N+1 mod FFTs
	mov	_INFT+1*SZPTR, rax
	IFNDEF X86_64
	mov	_INFT+2*SZPTR, OFFSET jmptable ; x86 mersenne mod FFTs
	mov	_INFT+3*SZPTR, OFFSET jmptablep ; x86 2^N+1 mod FFTs
	ENDIF
	mov	_GWVERNUM, VERSION_NUMBER
	ret
_gwinfo1 ENDP


; Setup routine.  Copy FFT constants set by C code, return procedure pointers

	PUBLIC	_gwsetup2
_gwsetup2 PROC
	push	rbp
	push	rbx
	push	rdi
	push	rsi

; Copy pointers to 1 premultiplier, 11 sin/cos, 3 normalization tables
; 1 carries table, 1 scratch area, and 1 plus1-premultiplier table.

	mov	rsi, OFFSET _GWPROCPTRS	; C code put ptrs here
	mov	rdi, OFFSET pass2_premults ; Addr to store table pointers
	mov	ecx, 18			; Copy 18 table ptrs
cp3:	mov	rax, [rsi]		; Load routine ptr
	mov	[rdi], rax		; Store routine ptr
	lea	rsi, [rsi+SZPTR]
	lea	rdi, [rdi+SZPTR]
	dec	ecx
	jnz	short cp3

; Copy jmptable counts common to x87 and sse2 code

	mov	rsi, _INFT		; Reload table pointer
	mov	eax, [rsi+16]		; Copy clm
	and	eax, 0ffffh
	mov	cache_line_multiplier, eax
	mov	eax, [rsi+28+5*SZPTR]	; Copy add/normalize counters
	mov	addcount1, eax
	mov	eax, [rsi+28+5*SZPTR+4]
	mov	normcount1, eax
	mov	eax, [rsi+28+5*SZPTR+8]	; Copy 5 counts
	mov	count1, eax
	mov	eax, [rsi+28+5*SZPTR+12]
	mov	count2, eax
	mov	eax, [rsi+28+5*SZPTR+16]
	mov	count3, eax
	mov	eax, [rsi+28+5*SZPTR+20]
	mov	count4, eax
	mov	eax, [rsi+28+5*SZPTR+24]
	mov	count5, eax

; Copy procedure pointers for C code to use

	mov	rsi, _INFT		; Reload table pointer
	lea	rsi, [rsi+28]
	mov	rdi, OFFSET _GWPROCPTRS
	mov	ecx, 4			; Copy 4 routine ptrs
cp4:	mov	rax, [rsi]		; Load routine ptr
	mov	[rdi], rax		; Store routine ptr
	lea	rsi, [rsi+SZPTR]
	lea	rdi, [rdi+SZPTR]
	dec	ecx
	jnz	short cp4

; Next copy add/sub/copy/normalization procedure pointers for C code to use

	mov	rsi, _INFT		; Reload table pointer
	mov	rsi, [rsi+28+4*SZPTR]	; Load prctab offset
	mov	rdi, OFFSET _GWPROCPTRS+4*SZPTR
	mov	ecx, 10			; Copy 10 routine ptrs
cp5:	mov	rax, [rsi]		; Load routine ptr
	mov	[rdi], rax		; Store routine ptr
	lea	rsi, [rsi+SZPTR]
	lea	rdi, [rdi+SZPTR]
	dec	ecx
	jnz	short cp5
	cmp	_RATIONAL_FFT, 0	; Are we doing rational FFTS?
	jne	short rat		; jump if yes
	lea	rsi, [rsi+12*SZPTR]	; Use irrational set of norm routines
rat:	cmp	_ZERO_PADDED_FFT, 0	; Are we doing zero padded FFTS?
	je	short nzpad		; jump if not
	lea	rsi, [rsi+8*SZPTR]	; Use zero-pad set of norm routines
nzpad:	mov	ecx, 8
cp6:	mov	rax, [rsi]		; Load norm routine ptr
	mov	[rdi], rax		; Store norm routine ptr
	lea	rsi, [rsi+SZPTR]
	lea	rdi, [rdi+SZPTR]
	dec	ecx
	jnz	short cp6

; Copy C computed constants

	mov	rsi, _SRCARG		;; Address of C computed constants
	copy_double_thrice MINUS_CARG, XMM_MINUS_CARG
	copy_float_once BIGVAL
	copy_float_once BIGBIGVAL
	copy_double_twice XMM_BIGVAL
	copy_double_twice XMM_BIGBIGVAL
	copy_double_thrice SQRTHALF, XMM_SQRTHALF
	copy_double_once ttmp_ff_inv
	copy_double_twice XMM_NORM012_FF
	copy_double_twice XMM_K_HI
	copy_double_twice XMM_K_LO
	copy_double_twice XMM_K_HI_2
	copy_double_twice XMM_K_HI_1
	copy_eight_doubles_four_times XMM_LIMIT_INVERSE
	copy_eight_doubles_four_times XMM_LIMIT_BIGMAX
	copy_eight_doubles_four_times XMM_LIMIT_BIGMAX_NEG
	copy_double_thrice P951, XMM_P951
	copy_double_thrice P618, XMM_P618
	copy_double_thrice P309, XMM_P309
	copy_double_thrice M262, XMM_M262
	copy_double_thrice P588, XMM_P588
	copy_double_thrice M162, XMM_M162
	copy_double_thrice M809, XMM_M809
	copy_double_thrice M382, XMM_M382
	copy_double_thrice P866, XMM_P866
	copy_double_thrice P433, XMM_P433
	copy_double_thrice P577, XMM_P577
	copy_double_thrice P975, XMM_P975
	copy_double_thrice P445, XMM_P445
	copy_double_thrice P180, XMM_P180
	copy_double_thrice P623, XMM_P623
	copy_double_thrice M358, XMM_M358
	copy_double_thrice P404, XMM_P404
	copy_double_once M223
	copy_double_once M901
	copy_double_once M691

; P4 SSE2 initialization is different than non-SSE2 initialization

IFNDEF X86_64
	test	_CPU_FLAGS, 0010h	; See if we should run P4 SSE2 routines
	jnz	p4init			; Yes, go to SSE2 init code

; Calculate constants used in two pass FFTs.  Foremost is the pass 1 blkdst
; and normalize blkdst for auxillary mult routines.  The two values are the
; same except for larger FFTs which use a scratch area.

	mov	esi, _INFT		; Reload table pointer
	mov	ecx, [esi+24]		; levels done in pass2
	cmp	ecx, 0			; skip calculations for one pass FFTs
	je	not2
	mov	eax, 16			; size of a complex number
	shl	eax, cl			; data size of pass2 block
	shr	eax, 12			; divide by 4KB
	mov	normval4, eax		; Num pages, used in normalized add/sub
	imul	eax, 4096+64		; pad 64 bytes every 4KB
	add	eax, 64			; pad 64 bytes between blocks
	mov	pass1blkdst, rax
	mov	normblkdst, rax
	mov	normblkdst8, 0
	cmp	DWORD PTR [esi+20], 0	; Is a scratch area used?
	je	short notscr		; No, normblkdst is OK.
	mov	eax, 32			; Compute scratch area normblkdst
	imul	eax, cache_line_multiplier
	mov	normblkdst, rax
	cmp	DWORD PTR [esi+28+5*SZPTR], 128	; Only larger pass1's have padding
	jl	short notscr		; Jump if small pass 1
	mov	normblkdst8, 64
notscr:
	mov	eax, 4096/32		; Cache lines in a page
	sub	edx, edx		; clear edx for divide instructions
	div	cache_line_multiplier
	mov	normval1, eax		; Save count of clms in a page 
	mov	eax, cache_line_multiplier
	shl	eax, 1
	mov	edx, addcount1
	dec	edx
	imul	eax, edx
	mov	normval2, rax		; Flags ptr fudge factor in add/sub
	mov	edx, cache_line_multiplier
	shl	edx, 1
	add	eax, edx
	imul	eax, normval1
	imul	eax, normval4
	sub	edx, eax
	mov	normval3, rdx		; Flags ptr fudge factor #2 in add/sub
not2:
	jmp	pop_and_ret
ENDIF

; P4 SSE2 initialization code

p4init:

	; Calculate pass 1 blkdst and normalize blkdst for xmult3 routines.
	; Small 2-pass FFTs use the in-place v22 blkdst of (65536+4096+128).
	; Larger FFTs use the scratch area with sporadic 128 byte gaps
	; in clmblkdst.
	cmp	_FFTLEN, 524288
	jg	short lg2p		; Jump if this is a larger FFT
	mov	rax, 65536+4096+128	; V22 blkdst
	mov	pass1blkdst, rax
	mov	edx, cache_line_multiplier ; clm
	shl	rdx, 6			; clm*64
	sub	rax, rdx		; Normblkdst = pass1blkdst - clm*64
	mov	normblkdst, rax
	mov	normblkdst8, 0
	jmp	short done3
lg2p:	mov	eax, cache_line_multiplier ; Calc pass1blkdst
	shl	rax, 7
	add	rax, 65536
	mov	pass1blkdst, rax	; 65536 + clm*128
	mov	normblkdst, 0		; Pad in clmblkdst is zero
	mov	normblkdst8, 128	; Pad for clmblkdst8 is 128
done3:

; Compute more normalization counters and constants

	sub	edx, edx
	mov	eax, _FFTLEN
	cmp	eax, 8192
	jle	short onep
	div	addcount1
	shr	eax, 3			; Cache lines in a block
	div	cache_line_multiplier
	mov	normval1, eax		; Used in normalized add/sub
	mov	eax, cache_line_multiplier
	shl	eax, 2
	mov	edx, addcount1
	dec	edx
	imul	eax, edx
	mov	normval2, rax		; Flags ptr fudge factor in add/sub
	mov	edx, cache_line_multiplier
	shl	edx, 2
	add	eax, edx
	imul	eax, normval1
	sub	rdx, rax
	mov	normval3, rdx		; Flags ptr fudge factor #2 in add/sub
onep:

; Return

pop_and_ret:
	pop	rsi
	pop	rdi
	pop	rbx
	pop	rbp
	ret
_gwsetup2 ENDP

;
; Set various globals when the multiplication constant changes
;

	PUBLIC	_eset_mul_const
_eset_mul_const PROC
	push	rsi
	push	rdi
	mov	rsi, _SRCARG		; Adderss of doubles to copy
	copy_double_twice XMM_MULCONST
	copy_double_twice XMM_MINUS_CARG_TIMES_MULCONST
	copy_double_twice XMM_K_TIMES_MULCONST_HI
	copy_double_twice XMM_K_TIMES_MULCONST_LO
	copy_double_twice XMM_K_TIMES_MULCONST_HI_2
	copy_double_twice XMM_K_TIMES_MULCONST_HI_1
	pop	rdi
	pop	rsi
	ret
_eset_mul_const ENDP

_TEXT ENDS
END
