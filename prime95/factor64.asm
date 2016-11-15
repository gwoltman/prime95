; Copyright 1995-2016 Mersenne Research, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements fast trial factoring of Mersenne numbers
;
; This only runs on 64-bit CPUs.  For 32-bit CPUs see factor32.asm
;

TITLE   factor64

INCLUDE	yarch.mac
INCLUDE	ybasics.mac
INCLUDE	unravel.mac
INCLUDE	factor64.mac

; In 32-bit mode we are so starved for registers that we are forced to
; use the stack pointer to access the asm_data.  In 64-bit mode we
; use one of the extra 8 registers.  We can't use the rsp trick because
; that violates the Window's exception handling stack unwind mechanism.

AD_BASE		EQU	<r11>

;
; Define the offsets into the C / assembly communication structure
;

p			EQU	QWORD PTR [AD_BASE+0*8]		; Exponent of Mersenne number being trial factored
facdists		EQU	QWORD PTR [AD_BASE+1*8]		; 64 distances between sieve factors
facdist64		EQU	QWORD PTR [AD_BASE+65*8]	
facdist12K		EQU	QWORD PTR [AD_BASE+66*8]	; Distance to first factor in next sieve
savefac1		EQU	QWORD PTR [AD_BASE+67*8]	; LSW of the first factor in the sieve
savefac0		EQU	QWORD PTR [AD_BASE+68*8]	; MSW of the first factor in the sieve
sieve			EQU	[AD_BASE+69*8]			; Area to sieve, or an already sieved area to TF
primearray		EQU	[AD_BASE+70*8]			; Array of primes
offsetarray		EQU	[AD_BASE+71*8]			; Array of bit-to-clear offsets
initsieve		EQU	[AD_BASE+72*8]			; Array used to initialize sieve
TWO_TO_FACSIZE_PLUS_62	EQU	QWORD PTR [AD_BASE+73*8]

FACHSW			EQU	DWORD PTR [AD_BASE+74*8+0*4]
FACMSW			EQU	DWORD PTR [AD_BASE+74*8+1*4]
FACLSW			EQU	DWORD PTR [AD_BASE+74*8+2*4]
CPU_FLAGS		EQU	DWORD PTR [AD_BASE+74*8+3*4]
SSE2_LOOP_COUNTER	EQU	DWORD PTR [AD_BASE+74*8+4*4]
FMA_INITVAL_TYPE	EQU	DWORD PTR [AD_BASE+74*8+5*4]
FMA_64_MINUS_LOW_BITS	EQU	DWORD PTR [AD_BASE+74*8+6*4]
initstart		EQU	DWORD PTR [AD_BASE+74*8+7*4]	; First byte in initsieve to copy

;; Common constants set in C code for SSE2, AVX2 when using
;; factoring algorithm that uses 30-bit integer arithmetic.

;;pad to next cache line to match facasm_data structure defined in commonb.c
YMM_INITVAL		EQU	DWORD PTR [AD_BASE+80*8]
XMM_INITVAL		EQU	DWORD PTR [AD_BASE+80*8]
XMM_INIT120BS		EQU	DWORD PTR [AD_BASE+84*8]
XMM_INITBS		EQU	DWORD PTR [AD_BASE+85*8]
XMM_BS			EQU	DWORD PTR [AD_BASE+86*8]
XMM_SHIFTER		EQU	DWORD PTR [AD_BASE+87*8]	; 32 quad word shifter values

;; SSE2 supports 128-bit operations on integers (four factors simultaneously)
;; Use the XMM registers and the memory area below
;; AVX2 supports 256-bit operations on integers (eight trial factors simultaneouly)
;; Use the YMM registers and the memory area below

XMM_INVFAC		EQU	DWORD PTR [AD_BASE+128*8]	;; Start at next cache line
XMM_INVFACa		EQU	DWORD PTR [AD_BASE+130*8]
XMM_I1			EQU	DWORD PTR [AD_BASE+132*8]
XMM_I1a			EQU	DWORD PTR [AD_BASE+134*8]
XMM_I2			EQU	DWORD PTR [AD_BASE+136*8]
XMM_I2a			EQU	DWORD PTR [AD_BASE+138*8]
XMM_F1			EQU	DWORD PTR [AD_BASE+140*8]
XMM_F1a			EQU	DWORD PTR [AD_BASE+142*8]
XMM_F2			EQU	DWORD PTR [AD_BASE+144*8]
XMM_F2a			EQU	DWORD PTR [AD_BASE+146*8]
XMM_F3			EQU	DWORD PTR [AD_BASE+148*8]
XMM_F3a			EQU	DWORD PTR [AD_BASE+150*8]
XMM_TWO_120_MODF1	EQU	DWORD PTR [AD_BASE+152*8]
XMM_TWO_120_MODF1a	EQU	DWORD PTR [AD_BASE+154*8]
XMM_TWO_120_MODF2	EQU	DWORD PTR [AD_BASE+156*8]
XMM_TWO_120_MODF2a	EQU	DWORD PTR [AD_BASE+158*8]
XMM_TWO_120_MODF3	EQU	DWORD PTR [AD_BASE+160*8]
XMM_TWO_120_MODF3a	EQU	DWORD PTR [AD_BASE+162*8]

YMM_INVFAC		EQU	DWORD PTR [AD_BASE+128*8]
YMM_INVFACa		EQU	DWORD PTR [AD_BASE+132*8]
YMM_I1			EQU	DWORD PTR [AD_BASE+136*8]
YMM_I1a			EQU	DWORD PTR [AD_BASE+140*8]
YMM_I2			EQU	DWORD PTR [AD_BASE+144*8]
YMM_I2a			EQU	DWORD PTR [AD_BASE+148*8]
YMM_F1			EQU	DWORD PTR [AD_BASE+152*8]
YMM_F1a			EQU	DWORD PTR [AD_BASE+156*8]
YMM_F2			EQU	DWORD PTR [AD_BASE+160*8]
YMM_F2a			EQU	DWORD PTR [AD_BASE+164*8]
YMM_F3			EQU	DWORD PTR [AD_BASE+168*8]
YMM_F3a			EQU	DWORD PTR [AD_BASE+172*8]
YMM_TWO_120_MODF1	EQU	DWORD PTR [AD_BASE+176*8]
YMM_TWO_120_MODF1a	EQU	DWORD PTR [AD_BASE+180*8]
YMM_TWO_120_MODF2	EQU	DWORD PTR [AD_BASE+184*8]
YMM_TWO_120_MODF2a	EQU	DWORD PTR [AD_BASE+188*8]
YMM_TWO_120_MODF3	EQU	DWORD PTR [AD_BASE+192*8]
YMM_TWO_120_MODF3a	EQU	DWORD PTR [AD_BASE+196*8]

;; Now the variables needed when using factoring algorithm that uses two doubles.
;;
;; AVX-FMA supports 256-bit operations on doubles (12 trial factors simultaneouly)
;; Use the YMM registers and the memory area below
;; AVX512-FMA supports 512-bit operations on doubles (32 trial factors simultaneouly)
;; Use the YMM registers and the memory area below

;;pad to next cache line to match facasm_data structure defined in commonb.c
FMA_SHIFTER			EQU	DWORD PTR [AD_BASE+80*8]	; Up to 32 dword shifter values
YMM_FMA_LOW_MASK		EQU	QWORD PTR [AD_BASE+96*8]
YMM_FMA_RND_CONST0		EQU	QWORD PTR [AD_BASE+100*8]	; 3 * 2^(51 + SPLIT_POINT)
YMM_FMA_RND_CONST1		EQU	QWORD PTR [AD_BASE+104*8]	; 3 * 2^(51 + BITS_IN_LO_WORD*2)
YMM_FMA_RND_CONST2		EQU	QWORD PTR [AD_BASE+108*8]	; 3 * 2^(51 + BITS_IN_LO_WORD)
YMM_FMA_RND_CONST4		EQU	QWORD PTR [AD_BASE+112*8]	; 3 * 2^(51 - SPLIT_POINT)
YMM_FMA_RND_CONST5		EQU	QWORD PTR [AD_BASE+116*8]	; 3 * 2^(51 + BITS_IN_LO_WORD - SPLIT_POINT)
YMM_FMA_INITVAL			EQU	QWORD PTR [AD_BASE+120*8]	; Initial squared value (a power of two)
YMM_FMA_TWOPOW_SPLIT		EQU	QWORD PTR [AD_BASE+124*8]	; 2^SPLIT_POINT
YMM_FMA_TWO_TO_LO		EQU	QWORD PTR [AD_BASE+128*8]	; 2^LOW_BITS
YMM_FMA_INVFAC			EQU	QWORD PTR [AD_BASE+132*8]
YMM_FMA_INVFACa			EQU	QWORD PTR [AD_BASE+136*8]
YMM_FMA_INVFACb			EQU	QWORD PTR [AD_BASE+140*8]
YMM_FMA_FAC_HI			EQU	QWORD PTR [AD_BASE+144*8]
YMM_FMA_FAC_HIa			EQU	QWORD PTR [AD_BASE+148*8]
YMM_FMA_FAC_HIb			EQU	QWORD PTR [AD_BASE+152*8]
YMM_FMA_FAC_LO			EQU	QWORD PTR [AD_BASE+156*8]
YMM_FMA_FAC_LOa			EQU	QWORD PTR [AD_BASE+160*8]
YMM_FMA_FAC_LOb			EQU	QWORD PTR [AD_BASE+164*8]
YMM_FMA_TWOPOW1_MODF_HI		EQU	QWORD PTR [AD_BASE+168*8]
YMM_FMA_TWOPOW2_MODF_HI		EQU	QWORD PTR [AD_BASE+172*8]
YMM_FMA_TWOPOW1_MODF_HIa	EQU	QWORD PTR [AD_BASE+176*8]
YMM_FMA_TWOPOW2_MODF_HIa	EQU	QWORD PTR [AD_BASE+180*8]
YMM_FMA_TWOPOW1_MODF_HIb	EQU	QWORD PTR [AD_BASE+184*8]
YMM_FMA_TWOPOW2_MODF_HIb	EQU	QWORD PTR [AD_BASE+188*8]
YMM_FMA_TWOPOW1_MODF_LO		EQU	QWORD PTR [AD_BASE+192*8]
YMM_FMA_TWOPOW2_MODF_LO		EQU	QWORD PTR [AD_BASE+196*8]
YMM_FMA_TWOPOW1_MODF_LOa	EQU	QWORD PTR [AD_BASE+200*8]
YMM_FMA_TWOPOW2_MODF_LOa	EQU	QWORD PTR [AD_BASE+204*8]
YMM_FMA_TWOPOW1_MODF_LOb	EQU	QWORD PTR [AD_BASE+208*8]
YMM_FMA_TWOPOW2_MODF_LOb	EQU	QWORD PTR [AD_BASE+212*8]
YMM_FMA_FAC_HI_INT		EQU	QWORD PTR [AD_BASE+216*8]
YMM_FMA_FAC_HI_INTa		EQU	QWORD PTR [AD_BASE+220*8]
YMM_FMA_FAC_HI_INTb		EQU	QWORD PTR [AD_BASE+224*8]
YMM_FMA_FAC_LO_INT		EQU	QWORD PTR [AD_BASE+228*8]
YMM_FMA_FAC_LO_INTa		EQU	QWORD PTR [AD_BASE+232*8]
YMM_FMA_FAC_LO_INTb		EQU	QWORD PTR [AD_BASE+236*8]
YMM_FMA_SAVE			EQU	QWORD PTR [AD_BASE+240*8]
YMM_FMA_SAVEa			EQU	QWORD PTR [AD_BASE+244*8]
YMM_FMA_SAVEb			EQU	QWORD PTR [AD_BASE+248*8]
YMM_FMA_SAVE2			EQU	QWORD PTR [AD_BASE+252*8]
YMM_FMA_SAVE2a			EQU	QWORD PTR [AD_BASE+256*8]
YMM_FMA_SAVE2b			EQU	QWORD PTR [AD_BASE+260*8]

;;pad to next cache line to match facasm_data structure defined in commonb.c
FMA_SHIFTER			EQU	DWORD PTR [AD_BASE+80*8]	; Up to 32 dword shifter values
ZMM_FMA_LOW_MASK		EQU	QWORD PTR [AD_BASE+96*8]
ZMM_FMA_RND_CONST0		EQU	QWORD PTR [AD_BASE+104*8]	; 3 * 2^(51 + SPLIT_POINT)
ZMM_FMA_RND_CONST1		EQU	QWORD PTR [AD_BASE+112*8]	; 3 * 2^(51 + BITS_IN_LO_WORD*2)
ZMM_FMA_RND_CONST2		EQU	QWORD PTR [AD_BASE+120*8]	; 3 * 2^(51 + BITS_IN_LO_WORD)
ZMM_FMA_RND_CONST4		EQU	QWORD PTR [AD_BASE+128*8]	; 3 * 2^(51 - SPLIT_POINT)
ZMM_FMA_RND_CONST5		EQU	QWORD PTR [AD_BASE+136*8]	; 3 * 2^(51 + BITS_IN_LO_WORD - SPLIT_POINT)
ZMM_FMA_INITVAL			EQU	QWORD PTR [AD_BASE+144*8]	; Initial squared value (a power of two)
ZMM_FMA_TWOPOW_SPLIT		EQU	QWORD PTR [AD_BASE+152*8]	; 2^SPLIT_POINT
ZMM_FMA_TWO_TO_LO		EQU	QWORD PTR [AD_BASE+160*8]	; 2^LOW_BITS
ZMM_FMA_INVFAC			EQU	QWORD PTR [AD_BASE+168*8]
ZMM_FMA_INVFACa			EQU	QWORD PTR [AD_BASE+176*8]
ZMM_FMA_INVFACb			EQU	QWORD PTR [AD_BASE+184*8]
ZMM_FMA_INVFACc			EQU	QWORD PTR [AD_BASE+192*8]
ZMM_FMA_FAC_HI			EQU	QWORD PTR [AD_BASE+200*8]
ZMM_FMA_FAC_HIa			EQU	QWORD PTR [AD_BASE+208*8]
ZMM_FMA_FAC_HIb			EQU	QWORD PTR [AD_BASE+216*8]
ZMM_FMA_FAC_HIc			EQU	QWORD PTR [AD_BASE+224*8]
ZMM_FMA_FAC_LO			EQU	QWORD PTR [AD_BASE+232*8]
ZMM_FMA_FAC_LOa			EQU	QWORD PTR [AD_BASE+240*8]
ZMM_FMA_FAC_LOb			EQU	QWORD PTR [AD_BASE+248*8]
ZMM_FMA_FAC_LOc			EQU	QWORD PTR [AD_BASE+256*8]
ZMM_FMA_TWOPOW1_MODF_HI		EQU	QWORD PTR [AD_BASE+264*8]
ZMM_FMA_TWOPOW2_MODF_HI		EQU	QWORD PTR [AD_BASE+272*8]
ZMM_FMA_TWOPOW1_MODF_HIa	EQU	QWORD PTR [AD_BASE+280*8]
ZMM_FMA_TWOPOW2_MODF_HIa	EQU	QWORD PTR [AD_BASE+288*8]
ZMM_FMA_TWOPOW1_MODF_HIb	EQU	QWORD PTR [AD_BASE+296*8]
ZMM_FMA_TWOPOW2_MODF_HIb	EQU	QWORD PTR [AD_BASE+304*8]
ZMM_FMA_TWOPOW1_MODF_HIc	EQU	QWORD PTR [AD_BASE+312*8]
ZMM_FMA_TWOPOW2_MODF_HIc	EQU	QWORD PTR [AD_BASE+320*8]
ZMM_FMA_TWOPOW1_MODF_LO		EQU	QWORD PTR [AD_BASE+328*8]
ZMM_FMA_TWOPOW2_MODF_LO		EQU	QWORD PTR [AD_BASE+336*8]
ZMM_FMA_TWOPOW1_MODF_LOa	EQU	QWORD PTR [AD_BASE+344*8]
ZMM_FMA_TWOPOW2_MODF_LOa	EQU	QWORD PTR [AD_BASE+352*8]
ZMM_FMA_TWOPOW1_MODF_LOb	EQU	QWORD PTR [AD_BASE+360*8]
ZMM_FMA_TWOPOW2_MODF_LOb	EQU	QWORD PTR [AD_BASE+368*8]
ZMM_FMA_TWOPOW1_MODF_LOc	EQU	QWORD PTR [AD_BASE+376*8]
ZMM_FMA_TWOPOW2_MODF_LOc	EQU	QWORD PTR [AD_BASE+384*8]
ZMM_FMA_FAC_HI_INT		EQU	QWORD PTR [AD_BASE+392*8]
ZMM_FMA_FAC_HI_INTa		EQU	QWORD PTR [AD_BASE+400*8]
ZMM_FMA_FAC_HI_INTb		EQU	QWORD PTR [AD_BASE+408*8]
ZMM_FMA_FAC_HI_INTc		EQU	QWORD PTR [AD_BASE+416*8]
ZMM_FMA_FAC_LO_INT		EQU	QWORD PTR [AD_BASE+424*8]
ZMM_FMA_FAC_LO_INTa		EQU	QWORD PTR [AD_BASE+432*8]
ZMM_FMA_FAC_LO_INTb		EQU	QWORD PTR [AD_BASE+440*8]
ZMM_FMA_FAC_LO_INTc		EQU	QWORD PTR [AD_BASE+448*8]

;; Non-SSE2 and non-AVX2 and non-AVX512 queues can share the SSE2/AVX2/AVX512/FMA memory

fac1			EQU	QWORD PTR [AD_BASE+128*8]	; Queued factors to test
fac2			EQU	QWORD PTR [AD_BASE+129*8]
fac3			EQU	QWORD PTR [AD_BASE+130*8]
fac4			EQU	QWORD PTR [AD_BASE+131*8]
fac1hi			EQU	QWORD PTR [AD_BASE+132*8]	; High word of queued factors to test
fac2hi			EQU	QWORD PTR [AD_BASE+133*8]
fac3hi			EQU	QWORD PTR [AD_BASE+134*8]
fac4hi			EQU	QWORD PTR [AD_BASE+135*8]

;; Other variables (start after the AVX512 variables)

shifter			EQU	QWORD PTR [AD_BASE+472*8]
shift66			EQU	QWORD PTR [AD_BASE+473*8]	; Two shift counts used in 66-bit factoring
initval1		EQU	QWORD PTR [AD_BASE+474*8]	; Lower 64 bits of initial value for squarer
initval0		EQU	QWORD PTR [AD_BASE+475*8]	; Upper 64 bits of initial value for squarer
initshift		EQU	QWORD PTR [AD_BASE+476*8]	; Initial shift count to make first quotient
initshift2		EQU	QWORD PTR [AD_BASE+477*8]	; Initial shift count to make first quotient
initdiv0		EQU	QWORD PTR [AD_BASE+478*8]	; Value to compute 1 / factor
shift_count		EQU	QWORD PTR [AD_BASE+479*8]	; Shift count used in squaring loop
temp			EQU	QWORD PTR [AD_BASE+480*8]
memshifter		EQU	QWORD PTR [AD_BASE+481*8]	; tlp66 case doesn't have a free register to store shifter
SAVED_REG1		EQU	[AD_BASE+482*8]
SAVED_REG2		EQU	[AD_BASE+483*8]
SAVED_REG3		EQU	[AD_BASE+484*8]
SAVED_REG4		EQU	[AD_BASE+485*8]
SAVED_REG5		EQU	[AD_BASE+486*8]
SAVED_REG6		EQU	[AD_BASE+487*8]
facinv2_in_memory	EQU	QWORD PTR [AD_BASE+488*8]

;
; Global variables
;

_GWDATA SEGMENT PAGE
	align 8
smfacendpt		DQ	100000000000h					; Limit for brute force factoring code (2^44)

masks19			DQ	(1 SHL 0) OR (1 SHL 19) OR (1 SHL 38) OR (1 SHL 57)
			DQ	(1 SHL 1) OR (1 SHL 20) OR (1 SHL 39) OR (1 SHL 58)
			DQ	(1 SHL 2) OR (1 SHL 21) OR (1 SHL 40) OR (1 SHL 59)
			DQ	(1 SHL 3) OR (1 SHL 22) OR (1 SHL 41) OR (1 SHL 60)
			DQ	(1 SHL 4) OR (1 SHL 23) OR (1 SHL 42) OR (1 SHL 61)
			DQ	(1 SHL 5) OR (1 SHL 24) OR (1 SHL 43) OR (1 SHL 62)
			DQ	(1 SHL 6) OR (1 SHL 25) OR (1 SHL 44) OR (1 SHL 63)
			DQ	(1 SHL 7) OR (1 SHL 26) OR (1 SHL 45)
			DQ	(1 SHL 8) OR (1 SHL 27) OR (1 SHL 46)
			DQ	(1 SHL 9) OR (1 SHL 28) OR (1 SHL 47)
			DQ	(1 SHL 10) OR (1 SHL 29) OR (1 SHL 48)
			DQ	(1 SHL 11) OR (1 SHL 30) OR (1 SHL 49)
			DQ	(1 SHL 12) OR (1 SHL 31) OR (1 SHL 50)
			DQ	(1 SHL 13) OR (1 SHL 32) OR (1 SHL 51)
			DQ	(1 SHL 14) OR (1 SHL 33) OR (1 SHL 52)
			DQ	(1 SHL 15) OR (1 SHL 34) OR (1 SHL 53)
			DQ	(1 SHL 16) OR (1 SHL 35) OR (1 SHL 54)
			DQ	(1 SHL 17) OR (1 SHL 36) OR (1 SHL 55)
			DQ	(1 SHL 18) OR (1 SHL 37) OR (1 SHL 56)

masks23			DQ	(1 SHL 0) OR (1 SHL 23) OR (1 SHL 46)
			DQ	(1 SHL 1) OR (1 SHL 24) OR (1 SHL 47)
			DQ	(1 SHL 2) OR (1 SHL 25) OR (1 SHL 48)
			DQ	(1 SHL 3) OR (1 SHL 26) OR (1 SHL 49)
			DQ	(1 SHL 4) OR (1 SHL 27) OR (1 SHL 50)
			DQ	(1 SHL 5) OR (1 SHL 28) OR (1 SHL 51)
			DQ	(1 SHL 6) OR (1 SHL 29) OR (1 SHL 52)
			DQ	(1 SHL 7) OR (1 SHL 30) OR (1 SHL 53)
			DQ	(1 SHL 8) OR (1 SHL 31) OR (1 SHL 54)
			DQ	(1 SHL 9) OR (1 SHL 32) OR (1 SHL 55)
			DQ	(1 SHL 10) OR (1 SHL 33) OR (1 SHL 56)
			DQ	(1 SHL 11) OR (1 SHL 34) OR (1 SHL 57)
			DQ	(1 SHL 12) OR (1 SHL 35) OR (1 SHL 58)
			DQ	(1 SHL 13) OR (1 SHL 36) OR (1 SHL 59)
			DQ	(1 SHL 14) OR (1 SHL 37) OR (1 SHL 60)
			DQ	(1 SHL 15) OR (1 SHL 38) OR (1 SHL 61)
			DQ	(1 SHL 16) OR (1 SHL 39) OR (1 SHL 62)
			DQ	(1 SHL 17) OR (1 SHL 40) OR (1 SHL 63)
			DQ	(1 SHL 18) OR (1 SHL 41)
			DQ	(1 SHL 19) OR (1 SHL 42)
			DQ	(1 SHL 20) OR (1 SHL 43)
			DQ	(1 SHL 21) OR (1 SHL 44)
			DQ	(1 SHL 22) OR (1 SHL 45)

masks29			DQ	(1 SHL 0) OR (1 SHL 29) OR (1 SHL 58)
			DQ	(1 SHL 1) OR (1 SHL 30) OR (1 SHL 59)
			DQ	(1 SHL 2) OR (1 SHL 31) OR (1 SHL 60)
			DQ	(1 SHL 3) OR (1 SHL 32) OR (1 SHL 61)
			DQ	(1 SHL 4) OR (1 SHL 33) OR (1 SHL 62)
			DQ	(1 SHL 5) OR (1 SHL 34) OR (1 SHL 63)
			DQ	(1 SHL 6) OR (1 SHL 35)
			DQ	(1 SHL 7) OR (1 SHL 36)
			DQ	(1 SHL 8) OR (1 SHL 37)
			DQ	(1 SHL 9) OR (1 SHL 38)
			DQ	(1 SHL 10) OR (1 SHL 39)
			DQ	(1 SHL 11) OR (1 SHL 40)
			DQ	(1 SHL 12) OR (1 SHL 41)
			DQ	(1 SHL 13) OR (1 SHL 42)
			DQ	(1 SHL 14) OR (1 SHL 43)
			DQ	(1 SHL 15) OR (1 SHL 44)
			DQ	(1 SHL 16) OR (1 SHL 45)
			DQ	(1 SHL 17) OR (1 SHL 46)
			DQ	(1 SHL 18) OR (1 SHL 47)
			DQ	(1 SHL 19) OR (1 SHL 48)
			DQ	(1 SHL 20) OR (1 SHL 49)
			DQ	(1 SHL 21) OR (1 SHL 50)
			DQ	(1 SHL 22) OR (1 SHL 51)
			DQ	(1 SHL 23) OR (1 SHL 52)
			DQ	(1 SHL 24) OR (1 SHL 53)
			DQ	(1 SHL 25) OR (1 SHL 54)
			DQ	(1 SHL 26) OR (1 SHL 55)
			DQ	(1 SHL 27) OR (1 SHL 56)
			DQ	(1 SHL 28) OR (1 SHL 57)

masks31			DQ	(1 SHL 0) OR (1 SHL 31) OR (1 SHL 62)
			DQ	(1 SHL 1) OR (1 SHL 32) OR (1 SHL 63)
			DQ	(1 SHL 2) OR (1 SHL 33)
			DQ	(1 SHL 3) OR (1 SHL 34)
			DQ	(1 SHL 4) OR (1 SHL 35)
			DQ	(1 SHL 5) OR (1 SHL 36)
			DQ	(1 SHL 6) OR (1 SHL 37)
			DQ	(1 SHL 7) OR (1 SHL 38)
			DQ	(1 SHL 8) OR (1 SHL 39)
			DQ	(1 SHL 9) OR (1 SHL 40)
			DQ	(1 SHL 10) OR (1 SHL 41)
			DQ	(1 SHL 11) OR (1 SHL 42)
			DQ	(1 SHL 12) OR (1 SHL 43)
			DQ	(1 SHL 13) OR (1 SHL 44)
			DQ	(1 SHL 14) OR (1 SHL 45)
			DQ	(1 SHL 15) OR (1 SHL 46)
			DQ	(1 SHL 16) OR (1 SHL 47)
			DQ	(1 SHL 17) OR (1 SHL 48)
			DQ	(1 SHL 18) OR (1 SHL 49)
			DQ	(1 SHL 19) OR (1 SHL 50)
			DQ	(1 SHL 20) OR (1 SHL 51)
			DQ	(1 SHL 21) OR (1 SHL 52)
			DQ	(1 SHL 22) OR (1 SHL 53)
			DQ	(1 SHL 23) OR (1 SHL 54)
			DQ	(1 SHL 24) OR (1 SHL 55)
			DQ	(1 SHL 25) OR (1 SHL 56)
			DQ	(1 SHL 26) OR (1 SHL 57)
			DQ	(1 SHL 27) OR (1 SHL 58)
			DQ	(1 SHL 28) OR (1 SHL 59)
			DQ	(1 SHL 29) OR (1 SHL 60)
			DQ	(1 SHL 30) OR (1 SHL 61)

masks37			DQ	(1 SHL 0) OR (1 SHL 37)
masks41			DQ	(1 SHL 0) OR (1 SHL 41)
masks43			DQ	(1 SHL 0) OR (1 SHL 43)
masks47			DQ	(1 SHL 0) OR (1 SHL 47)
masks53			DQ	(1 SHL 0) OR (1 SHL 53)
masks59			DQ	(1 SHL 0) OR (1 SHL 59)
masks61			DQ	(1 SHL 0) OR (1 SHL 61)

	align 16
XMM_LOWONE		DD	1,0,0,0
XMM_HIGHONE		DD	0,0,1,0
XMM_BITS28		DD	0FFFFFFFh,0,0FFFFFFFh,0
XMM_BITS30		DD	3FFFFFFFh,0,3FFFFFFFh,0
TWO_TO_30		DQ	1073741824.0				;; 2^30
	align 32
YMM_BITS28		DD	0FFFFFFFh,0,0FFFFFFFh,0,0FFFFFFFh,0,0FFFFFFFh,0
YMM_BITS30		DD	3FFFFFFFh,0,3FFFFFFFh,0,3FFFFFFFh,0,3FFFFFFFh,0
YMM_28TH_BIT		DD	08000000h,0,08000000h,0,08000000h,0,08000000h,0
YMM_ONE			DD	1,0,1,0,1,0,1,0
YMM_FMA_RND_CONST3	DQ 	6755399441055744.0,6755399441055744.0,6755399441055744.0,6755399441055744.0	; 3*2^51
YMM_FMA_ONE_OR_TWO TEXTEQU <YMM_FMA_ONE>
YMM_FMA_ONE		DQ	1.0,1.0,1.0,1.0
YMM_FMA_TWO		DQ	2.0,2.0,2.0,2.0
YMM_FMA_ZERO		DQ	0.0,0.0,0.0,0.0
	align 64
ZMM_FMA_RND_CONST3	DQ 	6755399441055744.0,6755399441055744.0,6755399441055744.0,6755399441055744.0	; 3*2^51
			DQ 	6755399441055744.0,6755399441055744.0,6755399441055744.0,6755399441055744.0
ZMM_FMA_ONE_OR_TWO TEXTEQU <ZMM_FMA_ONE>
ZMM_FMA_ONE		DQ	1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
ZMM_FMA_TWO		DQ	2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0
ZMM_FMA_ZERO		DQ	0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

;FOO1			DQ	15				;; The easy way to debug - enter the known factor here
;FOO2			DQ	1577855855906095289

	;; Align so that other GWDATA areas are aligned on a cache line
	align 128
_GWDATA ENDS


sievesize	EQU	00003000h	; 12KB sieve

;; Debugging macros

debug	MACRO ops:vararg
	IFDEF GDEBUG
	&ops
	ENDIF
	ENDM

assert	MACRO jcond
	LOCAL	ok
	IFDEF GDEBUG
	jcond	short ok
	sub	rsp, rsp
	pop	rax
ok:
	ENDIF
	ENDM


_TEXT	SEGMENT

; factor64_pass_setup (struct facasm_data *)
;	Do setup required when starting a new mod 120 pass
; Windows 64-bit (factor64_pass_setup)
;	Parameter ptr = rcx
; Linux 64-bit (factor64_pass_setup)
;	Parameter ptr = rdi

PROCF	factor64_pass_setup
	ad_prolog 0,0,rbx,rbp,rsi,rdi

; Pre-square the number as much as possible.  We have separate code
; paths for single-word and double-word trial factors.

	mov	rbx, savefac0		; Load factor
	mov	rcx, savefac1
	and	rbx, rbx		; Test MSW
	jz	short psq64		; If zero, presquaring 64-bits or less

; This is the 65-bits and above presquare code.
; Calculate number of bits in the trial factor

	bsr	rax, rbx		; Look for highest set bit
	mov	rdx, 63			; Create shift count
	sub	rdx, rax
	mov	shift_count, rdx	; Save shift count for squaring loop

; Compute the two shift counts used in 66-bit factoring.  The counts for
; 66-bits is (33-1)*256+1, for 67-bits is (33-2)*256+2, etc.

	mov	rbx, 33
	sub	rbx, rax
	shl	rbx, 8
	add	rbx, rax
	mov	shift66, rbx

; Pre-square the number as much as possible.  We could push this a little
; further because the tlp65 and tlp66 code doesn't really care that the
; first remainder is larger than the trial factor.

	add	rax, 147		; Compute maximum size of initval
	mov	rcx, p			; Load p
	sub	rbx, rbx		; This will be the shifter
setlp:	shrd	rbx, rcx, 1		; Move bits from p to shifter
	shr	rcx, 1			; Remove bit from p
	cmp	rcx, rax		; Is p > max initval
	jg	short setlp		; Yes, keep shifting
	mov	shifter, rbx		; Save shifter
	sub	rcx, 64			; Compute initval = 2^(rcx-64)
	jmp	short psqdn		; Join common presquare code

; Same code as above, but for 64-bit factors
; Calculate number of bits in the trial factor

psq64:	bsr	rax, rcx		; Look for highest set bit
	mov	rdx, 63			; Create shift count
	sub	rdx, rax
	mov	shift_count, rdx	; Save shift count for squaring loop

; Pre-square the number as much as possible.  That is, make sure initval
; divided by trail factor will not result in a 65-bit quotient.  The 64-bit
; factoring code can handle one more bit than 63-bit and less code.

	add	rax, 63			; Compute maximum size of initval
	sub	rbx, rbx		; Set rbx to 1 if 64-bit factors
	shld	rbx, rcx, 1
	add	rax, rbx		; Inc max initval size if 64-bit factor
	mov	rcx, p			; Load p
	sub	rbx, rbx		; This will be the shifter
setlp64:shrd	rbx, rcx, 1		; Move bits from p to shifter
	shr	rcx, 1			; Remove bit from p
	cmp	rcx, rax		; Is p > max initval
	jg	short setlp64		; Yes, keep shifting
	mov	shifter, rbx		; Save shifter

; Common pre-square cleanup code

psqdn:	mov	rax, 127		; Compute initial shift count
	sub	rax, rcx
	sub	rax, shift_count
	mov	initshift, rax		; Save initial shift count
	mov	rdx, shift66		; Compute alternate initial shift count
	shr	rdx, 8
	add	rax, rdx
	mov	initshift2, rax		; Save alternate initial shift count
	mov	rax, 1			; Compute initval = 2^rcx
	xor	rdx, rdx
	and	rcx, rcx
	jz	short shfdn
shflp:	shld	rdx, rax, 1
	shl	rax, 1
	dec	rcx
	jnz	short shflp
shfdn:	mov	initval0, rdx		; Save initval
	mov	initval1, rax
	mov	rcx, shift_count	; Compute initdiv = 2^(63-shift_count)
	mov	rbx, 1
	shl	rbx, 63
	shr	rbx, cl
	mov	initdiv0, rbx		; Save initdiv
	inc	rcx			; Increase shift_count because only
	mov	shift_count, rcx	; 2^63 / factor fits in 64 bits.  This
					; extra shift will give us the right
					; number of bits in a quotient
	ad_epilog 0,0,rbx,rbp,rsi,rdi
factor64_pass_setup ENDP


; factor64_small (struct facasm_data *)
;	Trial divide for small factors without using a sieve.  Recommended for TF below 2^44.
; Windows 64-bit (factor64_small)
;	Parameter ptr = rcx
; Linux 64-bit (factor64_small)
;	Parameter ptr = rdi

PROCF	factor64_small
	ad_prolog 0,0,rbx,rbp,rsi,rdi

;
; Brute force code to find small factors of 2**p - 1.  It works
; for trial factors of 63-bits or less.  However, the sieving algorithm
; is better if you are planning on testing a lot of trial factors.  So we
; only use this for 44-bit factors.
;

	mov	rcx, p
smslp:	add	rcx, rcx		; Shift until top bit on
	jns	short smslp
	mov	r8, rcx			; Compute shifter in r8
	shl	r8, 6
	shr	rcx, 58			; Compute initial value(?) in r9
	mov	r9, 1
	shl	r9, cl

; First trial factor is 2p + 1

	mov	rcx, p
	add	rcx, rcx
	inc	rcx

; If factor = 3 or 5 mod 8, then it can't be a factor of 2**p - 1

	mov	rdi, smfacendpt
smtest:	mov	rax, rcx
	and	al, 6
	jz	short smdoit
	cmp	al, 6
	jnz	short smnext

; Square the number until we computed 2**p MOD factor

smdoit:	mov	rbx, r8
	mov	rax, r9
	sub	rdx, rdx
	div	rcx
smloop:	mov	rax, rdx		; Square remainder
	mul	rax
	div	rcx			; Divide squared rem by trial factor
	add	rbx, rbx
	jnc	short smloop
	jz	short smexit
	add	rdx, rdx		; Double remainder
	jmp	short smloop

; Multiply remainder by two one last time (for the last carry out of shifter)
; If result = 1 mod factor, then we found a divisor of 2**p - 1

smexit:	add	rdx, rdx
	dec	rdx
	cmp	rdx, rcx
	jz	short smwin

; Try next possible factor

smnext:	add	rcx, p
	add	rcx, p
	cmp	rcx, rdi
	jl	short smtest

; No small factor found - return for ESC check

	mov	rax, 2			; Return for ESC check
	shr	rdi, 32
	mov	FACMSW, edi		; Restart after endpt
	jmp	short smdone

; Divisor found, return TRUE

smwin:	mov	rax, 1
	mov	FACLSW, ecx
	shr	rcx, 32
	mov	FACMSW, ecx
	jmp	short smdone

; Pop registers and return

smdone:	ad_epilog 0,0,rbx,rbp,rsi,rdi

factor64_small ENDP


; factor64_sieve (struct facasm_data *)
;	Sieve one 12KB chunk to eliminate potential factors that are not prime
; Windows 64-bit (factor64_sieve)
;	Parameter ptr = rcx
; Linux 64-bit (factor64_sieve)
;	Parameter ptr = rdi

PROCF	factor64_sieve
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r12,r13,r14,r15

;
; Init the sieve
;

	mov	rdi, sieve		; Load sieve address
	mov	rsi, initsieve		; Load sieve initialization bits addr
	mov	edx, initstart		; Load qword offset into initsieve
	add	rsi, rdx		; Compute copy address
	mov	rcx, sievesize / 8	; Count of qwords to copy
	rep	movsq			; Copy the qwords

;
; Loop through the small prime array, clearing sieve bits for very, very small primes
;

	mov	rdi, sieve		; Load address of sieve bits
	add	rdi, sievesize		; Compute end-of-sieve address
	mov	r9, offsetarray		; Ptr to first offset in array
	lea	rsi, masks19		; Load address to small prime masks

	mov	rbp, primearray		; Ptr to first prime in array
	cmp	DWORD PTR [rbp], 19	; Does the small prime array start at 19 or 23?
	jne	sv23andup		; Jump if sieving from 23 onwards

; Sieve primes 19 to 61

	mov	edx, [r9]		; Load bit offset to clear for prime 19
	mov	r12d, [r9+4]		; Load bit offset to clear for prime 23
	mov	r13d, [r9+8]		; Load bit offset to clear for prime 29
	mov	r14d, [r9+12]		; Load bit offset to clear for prime 31
	mov	r15d, [r9+16]		; Load bit offset to clear for prime 37
	mov	r8d, [r9+20]		; Load bit offset to clear for prime 41
	mov	r10d, [r9+24]		; Load bit offset to clear for prime 43
	mov	ebp, [r9+28]		; Load bit offset to clear for prime 47
	mov	bx, [r9+36]		; Load bit offset to clear for prime 59
	shl	rbx, 16
	mov	bx, [r9+32]		; Load bit offset to clear for prime 53
	mov	r9d, [r9+40]		; Load bit offset to clear for prime 61

	push	r11			; Save asm_data ptr
	mov	r11, -sievesize		; Load loop counter and index into sieve
smsvlp19:mov	rax, [rsi][rdx*8]	; Load mask for prime 19
	lea	rcx, [rdx+12]		; Calc next bit to clear by either adding (19 - 64 % 19)
	add	rdx, 12-19		; or adding (19 - 64 % 19) - 19
	cmovs	rdx, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	or	rax, [rsi+19*8][r12*8]	; Apply mask for prime 23
	lea	rcx, [r12+5]		; Calc next bit to clear by either adding (23 - 64 % 23)
	add	r12, 5-23		; or adding (23 - 64 % 23) - 23
	cmovs	r12, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	or	rax, [rsi+(19+23)*8][r13*8] ; Apply mask for prime 29
	lea	rcx, [r13+23]		; Calc next bit to clear by either adding (29 - 64 % 29)
	add	r13, 23-29		; or adding (29 - 64 % 29) - 29
	cmovs	r13, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	or	rax, [rsi+(19+23+29)*8][r14*8] ; Apply mask for prime 31
	lea	rcx, [r14+29]		; Calc next bit to clear by either adding (31 - 64 % 31)
	add	r14, 29-31		; or adding (31 - 64 % 31) - 31
	cmovs	r14, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, r15		; Apply mask for prime 37
	mov	r15, [rsi+(19+23+29+31)*8] ; Load mask if bit 0 is cleared
	shl	r15, cl			; Adjust the mask
	or	rax, r15		; Apply the mask
	lea	r15, [rcx+10]		; Calc next bit to clear by either adding (37 - 64 % 37)
	add	rcx, 10-37		; or adding (37 - 64 % 37) - 37
	cmovns	r15, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, r8			; Apply mask for prime 41
	mov	r8, [rsi+(19+23+29+31+1)*8] ; Load mask if bit 0 is cleared
	shl	r8, cl			; Adjust the mask
	or	rax, r8			; Apply the mask
	lea	r8, [rcx+18]		; Calc next bit to clear by either adding (41 - 64 % 41)
	add	rcx, 18-41		; or adding (41 - 64 % 41) - 41
	cmovns	r8, rcx			; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, r10		; Apply mask for prime 43
	mov	r10, [rsi+(19+23+29+31+2)*8] ; Load mask if bit 0 is cleared
	shl	r10, cl			; Adjust the mask
	or	rax, r10		; Apply the mask
	lea	r10, [rcx+22]		; Calc next bit to clear by either adding (43 - 64 % 43)
	add	rcx, 22-43		; or adding (43 - 64 % 43) - 43
	cmovns	r10, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, rbp		; Apply mask for prime 47
	mov	rbp, [rsi+(19+23+29+31+3)*8] ; Load mask if bit 0 is cleared
	shl	rbp, cl			; Adjust the mask
	or	rax, rbp		; Apply the mask
	lea	rbp, [rcx+30]		; Calc next bit to clear by either adding (47 - 64 % 47)
	add	rcx, 30-47		; or adding (47 - 64 % 47) - 47
	cmovns	rbp, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, r9			; Apply mask for prime 61
	mov	r9, [rsi+(19+23+29+31+6)*8] ; Load mask if bit 0 is cleared
	shl	r9, cl			; Adjust the mask
	or	rax, r9			; Apply the mask
	lea	r9, [rcx+58]		; Calc next bit to clear by either adding (61 - 64 % 61)
	add	rcx, 58-61		; or adding (61 - 64 % 61) - 61
	cmovns	r9, rcx			; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, rbx		; Apply mask for prime 53
	mov	rbx, [rsi+(19+23+29+31+4)*8] ; Load mask if bit 0 is cleared
	shl	rbx, cl			; Adjust the mask
	or	rax, rbx		; Apply the mask

	ror	rcx, 16			; Apply mask for prime 59
	mov	rbx, [rsi+(19+23+29+31+5)*8] ; Load mask if bit 0 is cleared
	shl	rbx, cl			; Adjust the mask
	or	rax, rbx		; Apply the mask

	not	rax			; Invert the mask (which is a list of bits to clear)
	and	QWORD PTR [rdi][r11], rax ; Apply the mask

	mov	rax, rcx		; Isolate bit to clear for prime 53
	shr	rax, 48
	lea	rbx, [rax+42]		; Calc next bit to clear by either adding (53 - 64 % 53)
	add	rax, 42-53		; or adding (53 - 64 % 53) - 53
	cmovs	rax, rbx		; Pick whichever is the smaller positive number as the bit to clear for next word

	and	rcx, 0FFh		; Isolate bit to clear for prime 59
	lea	rbx, [rcx+54]		; Calc next bit to clear by either adding (59 - 64 % 59)
	add	rcx, 54-59		; or adding (59 - 64 % 59) - 59
	cmovns	rbx, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word
	shl	rbx, 16			; Merge bits to clear for prime 53 and 59
	add	rbx, rax

	add	r11, 8
	js	smsvlp19

	pop	r11			; Restore asm_data ptr
	mov	rsi, offsetarray	; Ptr to first offset in array
	mov	[rsi], edx		; Save bit clear offset for prime 19 next time
	mov	[rsi+4], r12d		; Save bit clear offset for prime 23 next time
	mov	[rsi+8], r13d		; Save bit clear offset for prime 29 next time
	mov	[rsi+12], r14d		; Save bit clear offset for prime 31 next time
	mov	[rsi+16], r15d		; Save bit clear offset for prime 37 next time
	mov	[rsi+20], r8d		; Save bit clear offset for prime 41 next time
	mov	[rsi+24], r10d		; Save bit clear offset for prime 43 next time
	mov	[rsi+28], ebp		; Save bit clear offset for prime 47 next time
	mov	[rsi+32], bx		; Save bit clear offset for prime 53 next time
	shr	rbx, 16
	mov	[rsi+36], bx		; Save bit clear offset for prime 59 next time
	mov	[rsi+40], r9d		; Save bit clear offset for prime 61 next time

	mov	rbp, primearray		; Ptr to first prime in array
	add	rbp, 44			; Next prime in primearray
	add	rsi, 44			; Next offset in offsetarray
	jmp	sv67andup		; Go do the one small prime at a time sieving

; Sieve for small primes 23 to 61

sv23andup:
	mov	r12d, [r9]		; Load bit offset to clear for prime 23
	mov	r13d, [r9+4]		; Load bit offset to clear for prime 29
	mov	r14d, [r9+8]		; Load bit offset to clear for prime 31
	mov	r15d, [r9+12]		; Load bit offset to clear for prime 37
	mov	r8d, [r9+16]		; Load bit offset to clear for prime 41
	mov	r10d, [r9+20]		; Load bit offset to clear for prime 43
	mov	ebp, [r9+24]		; Load bit offset to clear for prime 47
	mov	ebx, [r9+28]		; Load bit offset to clear for prime 53
	mov	edx, [r9+32]		; Load bit offset to clear for prime 59
	mov	r9d, [r9+36]		; Load bit offset to clear for prime 61

	push	r11			; Save asm_data ptr
	mov	r11, -sievesize		; Load loop counter and index into sieve
smsvlp23:mov	rax, [rsi+19*8][r12*8]	; Apply mask for prime 23
	lea	rcx, [r12+5]		; Calc next bit to clear by either adding (23 - 64 % 23)
	add	r12, 5-23		; or adding (23 - 64 % 23) - 23
	cmovs	r12, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	or	rax, [rsi+(19+23)*8][r13*8] ; Apply mask for prime 29
	lea	rcx, [r13+23]		; Calc next bit to clear by either adding (29 - 64 % 29)
	add	r13, 23-29		; or adding (29 - 64 % 29) - 29
	cmovs	r13, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	or	rax, [rsi+(19+23+29)*8][r14*8] ; Apply mask for prime 31
	lea	rcx, [r14+29]		; Calc next bit to clear by either adding (31 - 64 % 31)
	add	r14, 29-31		; or adding (31 - 64 % 31) - 31
	cmovs	r14, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, r15		; Apply mask for prime 37
	mov	r15, [rsi+(19+23+29+31)*8] ; Load mask if bit 0 is cleared
	shl	r15, cl			; Adjust the mask
	or	rax, r15		; Apply the mask
	lea	r15, [rcx+10]		; Calc next bit to clear by either adding (37 - 64 % 37)
	add	rcx, 10-37		; or adding (37 - 64 % 37) - 37
	cmovns	r15, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, r8			; Apply mask for prime 41
	mov	r8, [rsi+(19+23+29+31+1)*8] ; Load mask if bit 0 is cleared
	shl	r8, cl			; Adjust the mask
	or	rax, r8			; Apply the mask
	lea	r8, [rcx+18]		; Calc next bit to clear by either adding (41 - 64 % 41)
	add	rcx, 18-41		; or adding (41 - 64 % 41) - 41
	cmovns	r8, rcx			; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, r10		; Apply mask for prime 43
	mov	r10, [rsi+(19+23+29+31+2)*8] ; Load mask if bit 0 is cleared
	shl	r10, cl			; Adjust the mask
	or	rax, r10		; Apply the mask
	lea	r10, [rcx+22]		; Calc next bit to clear by either adding (43 - 64 % 43)
	add	rcx, 22-43		; or adding (43 - 64 % 43) - 43
	cmovns	r10, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, rbp		; Apply mask for prime 47
	mov	rbp, [rsi+(19+23+29+31+3)*8] ; Load mask if bit 0 is cleared
	shl	rbp, cl			; Adjust the mask
	or	rax, rbp		; Apply the mask
	lea	rbp, [rcx+30]		; Calc next bit to clear by either adding (47 - 64 % 47)
	add	rcx, 30-47		; or adding (47 - 64 % 47) - 47
	cmovns	rbp, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, rbx		; Apply mask for prime 53
	mov	rbx, [rsi+(19+23+29+31+4)*8] ; Load mask if bit 0 is cleared
	shl	rbx, cl			; Adjust the mask
	or	rax, rbx		; Apply the mask
	lea	rbx, [rcx+42]		; Calc next bit to clear by either adding (53 - 64 % 53)
	add	rcx, 42-53		; or adding (53 - 64 % 53) - 53
	cmovns	rbx, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, rdx		; Apply mask for prime 59
	mov	rdx, [rsi+(19+23+29+31+5)*8] ; Load mask if bit 0 is cleared
	shl	rdx, cl			; Adjust the mask
	or	rax, rdx		; Apply the mask
	lea	rdx, [rcx+54]		; Calc next bit to clear by either adding (59 - 64 % 59)
	add	rcx, 54-59		; or adding (59 - 64 % 59) - 59
	cmovns	rdx, rcx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rcx, r9			; Apply mask for prime 61
	mov	r9, [rsi+(19+23+29+31+6)*8] ; Load mask if bit 0 is cleared
	shl	r9, cl			; Adjust the mask
	or	rax, r9			; Apply the mask
	lea	r9, [rcx+58]		; Calc next bit to clear by either adding (61 - 64 % 61)
	add	rcx, 58-61		; or adding (61 - 64 % 61) - 61
	cmovns	r9, rcx			; Pick whichever is the smaller positive number as the bit to clear for next word

	not	rax			; Invert the mask (which is a list of bits to clear)
	and	QWORD PTR [rdi][r11], rax ; Apply the mask

	add	r11, 8
	js	smsvlp23

	pop	r11			; Restore asm_data ptr
	mov	rsi, offsetarray	; Ptr to first offset in array
	mov	[rsi], r12d		; Save bit clear offset for prime 23 next time
	mov	[rsi+4], r13d		; Save bit clear offset for prime 29 next time
	mov	[rsi+8], r14d		; Save bit clear offset for prime 31 next time
	mov	[rsi+12], r15d		; Save bit clear offset for prime 37 next time
	mov	[rsi+16], r8d		; Save bit clear offset for prime 41 next time
	mov	[rsi+20], r10d		; Save bit clear offset for prime 43 next time
	mov	[rsi+24], ebp		; Save bit clear offset for prime 47 next time
	mov	[rsi+28], ebx		; Save bit clear offset for prime 53 next time
	mov	[rsi+32], edx		; Save bit clear offset for prime 59 next time
	mov	[rsi+36], r9d		; Save bit clear offset for prime 61 next time

	mov	rbp, primearray		; Ptr to first prime in array
	add	rbp, 40			; Next prime in primearray
	add	rsi, 40			; Next offset in offsetarray

;
; Sieve small primes 67 to 101
; On Skylake there was no advantage to going above primes 101 (we can push registers
; rbp, rsi, r11 to free up three registers to do primes 103, 107, and 109)
; However, accurate Skylake timings are hard to come by, so it may well pay
; to do more small primes.
;

sv67andup:
	mov	r12d, [rsi]		; Load bit offset to clear for prime 67
	mov	r13d, [rsi+4]		; Load bit offset to clear for prime 71
	mov	r14d, [rsi+8]		; Load bit offset to clear for prime 73
	mov	r15d, [rsi+12]		; Load bit offset to clear for prime 79
	mov	r8d, [rsi+16]		; Load bit offset to clear for prime 83
	mov	r10d, [rsi+20]		; Load bit offset to clear for prime 89
	mov	r9d, [rsi+24]		; Load bit offset to clear for prime 97
	mov	ebx, [rsi+28]		; Load bit offset to clear for prime 101
;	push	r11			; Save asm_data ptr
;	mov	r11d, [rsi+32]		; Load bit offset to clear for prime 103
;	push	rbp			; Save primearray ptr
;	mov	ebp, [rsi+36]		; Load bit offset to clear for prime 107
;	push	rsi			; Save offsetarray ptr
;	mov	esi, [rsi+40]		; Load bit offset to clear for prime 109

	mov	rcx, -sievesize		; Load loop counter and index into sieve
smsvlp67:
	mov	rax, [rdi][rcx]		; Get the sieve word

	mov	rdx, rax		; Copy sieve word
	btr	rax, r12		; Tentatively clear bit for prime 67
	sub	r12, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [r12+67]		; Alternative bit to clear adds the small prime 67
	cmovs	r12, rdx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rdx, rax		; Copy sieve word
	btr	rax, r13		; Tentatively clear bit for prime 71
	sub	r13, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [r13+71]		; Alternative bit to clear adds the small prime 71
	cmovs	r13, rdx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rdx, rax		; Copy sieve word
	btr	rax, r14		; Tentatively clear bit for prime 73
	sub	r14, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [r14+73]		; Alternative bit to clear adds the small prime 73
	cmovs	r14, rdx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rdx, rax		; Copy sieve word
	btr	rax, r15		; Tentatively clear bit for prime 79
	sub	r15, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [r15+79]		; Alternative bit to clear adds the small prime 79
	cmovs	r15, rdx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rdx, rax		; Copy sieve word
	btr	rax, r8			; Tentatively clear bit for prime 83
	sub	r8, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [r8+83]		; Alternative bit to clear adds the small prime 83
	cmovs	r8, rdx			; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rdx, rax		; Copy sieve word
	btr	rax, r10		; Tentatively clear bit for prime 89
	sub	r10, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [r10+89]		; Alternative bit to clear adds the small prime 89
	cmovs	r10, rdx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rdx, rax		; Copy sieve word
	btr	rax, r9			; Tentatively clear bit for prime 97
	sub	r9, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [r9+97]		; Alternative bit to clear adds the small prime 97
	cmovs	r9, rdx			; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rdx, rax		; Copy sieve word
	btr	rax, rbx		; Tentatively clear bit for prime 101
	sub	rbx, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [rbx+101]		; Alternative bit to clear adds the small prime 101
	cmovs	rbx, rdx		; Pick whichever is the smaller positive number as the bit to clear for next word

IFDEF DO_THREE_MORE_SMALL_PRIMES
	mov	rdx, rax		; Copy sieve word
	btr	rax, r11		; Tentatively clear bit for prime 103
	sub	r11, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [r11+103]		; Alternative bit to clear adds the small prime 103
	cmovs	r11, rdx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rdx, rax		; Copy sieve word
	btr	rax, rbp		; Tentatively clear bit for prime 107
	sub	rbp, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [rbp+107]		; Alternative bit to clear adds the small prime 107
	cmovs	rbp, rdx		; Pick whichever is the smaller positive number as the bit to clear for next word

	mov	rdx, rax		; Copy sieve word
	btr	rax, rsi		; Tentatively clear bit for prime 109
	sub	rsi, 64			; Subtract 64 from bit-to-clear offset
	cmovns	rax, rdx		; If negative then we want to use the sieve word with the tenatively cleared bit
	lea	rdx, [rsi+109]		; Alternative bit to clear adds the small prime 109
	cmovs	rsi, rdx		; Pick whichever is the smaller positive number as the bit to clear for next word
ENDIF

	mov	[rdi][rcx], rax		; Save the new sieve word
	add	rcx, 8
	js	smsvlp67

;	mov	rcx, r11		; Save bit-to-clr for prime 103
;	mov	rdx, rbp		; Save bit-to-clr for prime 107
;	mov	rax, rsi		; Save bit-to-clr for prime 109
;	pop	rsi			; Restore offset array ptr
;	pop	rbp			; Restore prime array ptr
;	pop	r11			; Restore asm_data ptr
	mov	[rsi], r12d		; Save bit clear offset for prime 67 next time
	mov	[rsi+4], r13d		; Save bit clear offset for prime 71 next time
	mov	[rsi+8], r14d		; Save bit clear offset for prime 73 next time
	mov	[rsi+12], r15d		; Save bit clear offset for prime 79 next time
	mov	[rsi+16], r8d		; Save bit clear offset for prime 83 next time
	mov	[rsi+20], r10d		; Save bit clear offset for prime 89 next time
	mov	[rsi+24], r9d		; Save bit clear offset for prime 97 next time
	mov	[rsi+28], ebx		; Save bit clear offset for prime 101 next time
;	mov	[rsi+32], ecx		; Save bit clear offset for prime 103 next time
;	mov	[rsi+36], edx		; Save bit clear offset for prime 107 next time
;	mov	[rsi+40], eax		; Save bit clear offset for prime 109 next time

	add	rbp, 32;+12		; Next prime in primearray
	add	rsi, 32;+12		; Next offset in offsetarray

;
; Loop through the rest of small prime array, clearing sieve bits.
;

sievelp:mov	ebx, [rbp]		; Load small prime
	and	rbx, rbx		; See if primearray fully processed
	jz	short sievedn		; Yes, break out of loop
	mov	eax, [rsi]		; Load bit offset to clear
	sub	rax, sievesize*8	; Make offset relative to end of the sieve
	jns	short clrdn		; If offset is positive there are no bits to clear
clrlp:	btr	[rdi], rax		; Clear the sieve bit
	add	rax, rbx		; Next bit to clear
	js	short clrlp		; Past the end of the sieve when rax goes positive
clrdn:	mov	[rsi], eax		; Save bit clear offset for next time
	lea	rbp, [rbp+4]		; Next primearray address
	lea	rsi, [rsi+4]		; Next offsetarray address
	jmp	short sievelp		; Work on next small prime

; Bump the first factor in sieve for the caller (caller looks at the value to decide if more sieving is needed)

sievedn:
	mov	rax, facdist12K
	add	savefac1, rax
	adc	savefac0, 0

; Pop registers and return

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r12,r13,r14,r15

factor64_sieve ENDP


; factor64_tf (struct facasm_data *)
;	Take a pre-sieved array and test potential factors
; Windows 64-bit (factor64_tf)
;	Parameter ptr = rcx
; Linux 64-bit (factor64_tf)
;	Parameter ptr = rdi

PROCF	factor64_tf
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r12,r13,r14,r15,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

; Branch to proper code to do the TF

	mov	rbx, savefac0		; Load trial factor corresponding
	mov	rcx, savefac1		; to first sieve bit
	test	CPU_FLAGS, 8000h	; Is this an FMA-capable machine?
	jz	short notfma		; No, use slower methods

	test	CPU_FLAGS, 100000h	; Is this an AVX512 machine?
	jnz	zftlp96			; Yes, use 45-92 bit (or more) AVX512-FMA code
	jmp	ftlp96			; No, use 45-92 bit (or more) AVX-FMA code

notfma:
	and	rbx, rbx		; Are we testing 65+ bit factors?
	jz	short oneword		; No, go do 64 bits or less
	test	CPU_FLAGS, 20000h	; Is this an AVX2 machine?
	jz	short notavx2		; No, use non-AVX2 code
	cmp	rbx, 400h		; Are we testing 75+ bit numbers?
	jge	atlp86			; Yes, use 75-86 bit AVX2 code (obsoleted by FMA code)
	jmp	atlp74			; No, use 65-74 bit AVX2 code (obsoleted by FMA code)

notavx2:
	test	CPU_FLAGS, 200h		; Is this an SSE2 machine?
	jz	short notsse2		; No, use non-SSE2 code
	cmp	rbx, 400h		; Are we testing 75+ bit numbers?
	jge	tlp86			; Yes, use 75-86 bit SSE2 code
	jmp	tlp74			; No, use 65-74 bit SSE2 code

notsse2:
	cmp	rbx, 1			; Are we testing 66+ bit factors?
	jg	tlp66			; Yes, go do it
	je	tlp65			; Jump if testing 65 bit factors?

oneword:
	mov	rax, rcx
	shr	rax, 58			; Testing 58 bit or less factors?
	jz	tlp58			; Yes, go test 58-bit factors
	shr	rax, 2			; Testing 60 bit or less factors?
	jz	tlp60			; Yes, go test 60-bit factors
	shr	rax, 1			; Testing 61 bit or less factors?
	jz	tlp61			; Yes, go test 61-bit factors
	shr	rax, 1			; Testing 62 bit or less factors?
	jz	tlp62			; Yes, go test 62-bit factors
	shr	rax, 1			; Testing 63 bit or less factors?
	jz	tlp63			; Yes, go test 63-bit factors
	jmp	tlp64			; No, go test 64-bit factors

;
; Make register assignments for one word remainders
;

facinv1	EQU	r13
facinv2	EQU	facinv2_in_memory	;; Formerly r11
facinv3	EQU	r15
temp3	EQU	rdi
rem3	EQU	r8
fac3reg	EQU	r9
temp2	EQU	r10
rem2	EQU	r14
fac2reg	EQU	r12
temp1	EQU	rbp
rem1	EQU	rbx
fac1reg	EQU	rsi

;
; Make different register assignments for two word remainders
;

rem3hi	EQU	rem3
rem3lo	EQU	fac3reg
rem2hi	EQU	rem2
rem2lo	EQU	fac2reg
rem1hi	EQU	rem1
rem1lo	EQU	fac1reg

;***********************************************************************
; Test up to 58 bit factors (actually handles most 59 bit factors too)
;***********************************************************************

tlp58:	mov	rsi, sieve		; Sieve address
	lea	rbp, [rsi+sievesize]	; Sieve end address
	xor	rdi, rdi		; Count of queued factors to be tested
svlp58:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf58:	bsf	rdx, rax		; Look for a set bit
	jnz	short test58		; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	cmp	rsi, rbp		; End of sieve?
	jl	short svlp58		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem58		; Yes, test them

; Save state

	mov	savefac1, rcx		; Save for the restart or more sieving

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem58:	mov	r8, fac1		; Copy first trial factor in queue
	mov	fac1[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	short rem58		; No, go copy another
	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
					; it will be re-added when we branch back to bsf58
	jmp	short do58		; Go test the potential factors

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test58:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	add	r8, facdists[rdx*8]	; Determine factor to test
	mov	fac1[rdi*8], r8		; Save the factor to test
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	bsf58			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

do58:	mov	SAVED_REG1, rsi		; Save sieve testing registers
	mov	SAVED_REG2, rax
	mov	SAVED_REG3, rcx
	mov	SAVED_REG4, rbp

; Precompute shifted 1 / factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac1reg, fac1
	div	fac1reg			; Compute 1 / trial_factor
	mov	facinv1, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac2reg, fac2
	div	fac2reg			; Compute 1 / trial_factor
	mov	facinv2, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac3reg, fac3
	div	fac3reg			; Compute 1 / trial_factor
	mov	facinv3, rax		; Save 1 / trial_factor

;
; Perform a division on the initial value to get started.
;

	mov	rcx, initshift		; Load initial shift count

	mov	rem1, initval1		; Load initial value
	mov	rem2, rem1		; Load initial value
	mov	rem3, rem1		; Load initial value

	mov	rax, facinv1		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac1reg			; Compute quotient * factor
	sub	rem1, rax		; Subtract from initval
debug	mov	temp1, initval0
debug	sbb	temp1, rdx
assert	jz

	mov	rax, facinv2		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac2reg			; Compute quotient * factor
	sub	rem2, rax		; Subtract from initval
debug	mov	temp2, initval0
debug	sbb	temp2, rdx
assert	jz

	mov	rax, facinv3		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac3reg			; Compute quotient * factor
	sub	rem3, rax		; Subtract from initval
debug	mov	temp3, initval0
debug	sbb	temp3, rdx
assert	jz

;
; Square remainder and get new remainder.  Input remainder is between
; zero and 6 * factor.  The squared remainder is up to 36 * factor^2.  The 
; quotient is up to 36 * factor.  For the shrd value to fit in 64 bits, the
; maximum factor is 2^63 / 36 = 2^57.  However, this estimate is too
; pessimistic.  It can be shown that the quotient is off by 2 only when facinv
; is large.  When facinv is large, the factor is small, and there is one fewer
; bit in the squared remainder, allowing 58 bits to work safely. 
;

	mov	temp1, shifter		; Load shifter
;	mov	rcx, shift_count	; Load squared value shift count
;	inc	rcx			; Compute an extra quotient bit
					; for better branch prediction
mov rcx,64
sub rcx,shift_count
sqlp58:	mov	rax, rem1		; Load remainder
	imul	rax			; Square remainder
	mov	rem1, rax		; Save squared remainder
	shrd	rax, rdx, cl		; Shift squared_upper
debug	shr	rdx, cl
assert	jz
	mul	facinv1
	imul	rdx, fac1reg		; Compute quotient * factor
	sub	rem1, rdx		; Subtract from squared remainder

	mov	rax, rem2		; Load remainder
	imul	rax			; Square remainder
	mov	rem2, rax		; Save squared remainder
	shrd	rax, rdx, cl		; Shift squared_upper
debug	shr	rdx, cl
assert	jz
	mul	facinv2
	imul	rdx, fac2reg		; Compute quotient * factor
	sub	rem2, rdx		; Subtract from squared remainder

	mov	rax, rem3		; Load remainder
	imul	rax			; Square remainder
	mov	rem3, rax		; Save squared remainder
	shrd	rax, rdx, cl		; Shift squared_upper
debug	shr	rdx, cl
assert	jz
	mul	facinv3
	imul	rdx, fac3reg		; Compute quotient * factor
	sub	rem3, rdx		; Subtract from squared remainder

;
; At this point, since quotient can be off by 2, the remainder is between
; zero and 3 * factor.   Multiply by two if necessary, making the remainder
; between zero and 6 * factor.  Test for end of squaring loop
;

	add	temp1, temp1		; One squaring completed, shift
	jnc	sqlp58			; If mul by 2 not needed, jump to loop
	jz	short exit58		; Are we done squaring?
	add	rem1, rem1		; Multiply remainder by 2
	add	rem2, rem2		; Multiply remainder by 2
	add	rem3, rem3		; Multiply remainder by 2
	jmp	sqlp58			; Do next iteration

; Make remainder between 0 and factor.  Multiply remainder by two one last
; time (for the last carry out of shifter).  If result = 1 mod factor, then
; we found a divisor of 2^p - 1

exit58:	mov	rax, rem1		; Make rem1 between 0 and 2*fac
	sub	rax, fac1reg
	cmovns	rem1, rax
	mov	rcx, rem2		; Make rem2 between 0 and 2*fac
	sub	rcx, fac2reg
	cmovns	rem2, rcx
	mov	rdx, rem3		; Make rem3 between 0 and 2*fac
	sub	rdx, fac3reg
	cmovns	rem3, rdx
	sub	rax, fac1reg		; Make rem1 between 0 and fac
	cmovns	rem1, rax
	sub	rcx, fac2reg		; Make rem2 between 0 and fac
	cmovns	rem2, rcx
	sub	rdx, fac3reg		; Make rem3 between 0 and fac
	cmovns	rem3, rdx

debug	cmp	rem1, fac1reg		; Remainder should be < factor
assert	jb
debug	cmp	rem2, fac2reg		; Remainder should be < factor
assert	jb
debug	cmp	rem3, fac3reg		; Remainder should be < factor
assert	jb

	add	rem1, rem1		; Double the remainder
	add	rem2, rem2		; Double the remainder
	add	rem3, rem3		; Double the remainder

	sub	rem1, fac1reg		; If rem == fac+1, we found a factor!
	dec	rem1
	jz	short win58
	sub	rem2, fac2reg		; If rem == fac+1, we found a factor!
	dec	rem2
	jz	short win58a
	sub	rem3, fac3reg		; If rem == fac+1, we found a factor!
	dec	rem3
	jz	short win58b

	sub	rdi, rdi		; Clear queued factors count

	mov	rbp, SAVED_REG4		; Restore sieve testing register
	mov	rcx, SAVED_REG3
	mov	rax, SAVED_REG2
	mov	rsi, SAVED_REG1
	jmp	bsf58			; Test next factor from sieve

win58:	mov	rax, fac1reg
	jmp	short win58c
win58a:	mov	rax, fac2reg
	jmp	short win58c
win58b:	mov	rax, fac3reg
win58c:	mov	FACLSW, eax		; Factor found!!! Return TRUE
	shr	rax, 32
	mov	FACMSW, eax
	mov	rax, 1
	jmp	done

;***********************************************************************
; Test up to 60 bit factors
;***********************************************************************

tlp60:	mov	rsi, sieve		; Sieve address
	lea	rbp, [rsi+sievesize]	; Sieve end address
	mov	rbx, 1			; Compute maximum trial factor
	shl	rbx, 60
	xor	rdi, rdi		; Count of queued factors to be tested
svlp60:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf60:	bsf	rdx, rax		; Look for a set bit
	jnz	short test60		; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	cmp	rsi, rbp		; End of sieve?
	jl	short svlp60		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem60a		; Yes, test them

; Save state

	mov	savefac1, rcx		; Save for the restart or more sieving

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle 60-bit overflow

oflow60:and	rdi, rdi		; Are there untested factors?
	jnz	short rem60		; Yes, go do remaining trial factors
	mov	FACMSW, 10000000h	; Return end point
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem60a:	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
					; it will be re-added when we branch back to bsf60
rem60:	mov	r8, fac1		; Copy first trial factor in queue
	mov	fac1[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	short rem60		; No, go copy another
	jmp	short do60		; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test60:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	add	r8, facdists[rdx*8]	; Determine factor to test
	cmp	r8, rbx			; Test for 60-bit overflow
	jae	short oflow60		; Jump if 60-bit overflow
	mov	fac1[rdi*8], r8		; Save the factor to test
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	bsf60			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

do60:	mov	SAVED_REG1, rsi		; Save sieve testing registers
	mov	SAVED_REG2, rax
	mov	SAVED_REG3, rbx
	mov	SAVED_REG4, rcx
	mov	SAVED_REG5, rbp

; Precompute shifted 1 / factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac1reg, fac1
	div	fac1reg			; Compute 1 / trial_factor
	mov	facinv1, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac2reg, fac2
	div	fac2reg			; Compute 1 / trial_factor
	mov	facinv2, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac3reg, fac3
	div	fac3reg			; Compute 1 / trial_factor
	mov	facinv3, rax		; Save 1 / trial_factor

;
; Perform a division on the initial value to get started.
;

	mov	rcx, initshift		; Load initial shift count

	mov	rem1, initval1		; Load initial value
	mov	rem2, rem1		; Load initial value
	mov	rem3, rem1		; Load initial value

	mov	rax, facinv1		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac1reg			; Compute quotient * factor
	sub	rem1, rax		; Subtract from initval
debug	mov	temp1, initval0
debug	sbb	temp1, rdx
assert	jz

	mov	rax, facinv2		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac2reg			; Compute quotient * factor
	sub	rem2, rax		; Subtract from initval
debug	mov	temp2, initval0
debug	sbb	temp2, rdx
assert	jz

	mov	rax, facinv3		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac3reg			; Compute quotient * factor
	sub	rem3, rax		; Subtract from initval
debug	mov	temp3, initval0
debug	sbb	temp3, rdx
assert	jz

;
; Square remainder and get new remainder.  Input remainder is between
; -2 * factor and 4 * factor.  The squared remainder is up to 16 * factor^2.
; quotient is up to 16 * factor.  For the shrd value to fit in 64 bits, the
; maximum factor is 2^63 / 16 = 2^59.  However, this estimate is too
; pessimistic.  It can be shown that the quotient is off by 2 only when facinv
; is large.  When facinv is large, the factor is small, there is one fewer
; bit in the squared remainder, allowing 60 bits to work safely. 
;

	mov	temp1, shifter		; Load shifter
;	mov	rcx, shift_count	; Load squared value shift count
;	inc	rcx			; Compute an extra quotient bit
					; for better branch prediction
mov rcx,64
sub rcx,shift_count
sqlp60:	mov	rax, rem1		; Load remainder
	imul	rax			; Square remainder
	mov	rem1, rax		; Save squared remainder
	shrd	rax, rdx, cl		; Shift squared_upper
debug	shr	rdx, cl
assert	jz
	mul	facinv1
	imul	rdx, fac1reg		; Compute quotient * factor
	sub	rem1, rdx		; Subtract from squared remainder

	mov	rax, rem2		; Load remainder
	imul	rax			; Square remainder
	mov	rem2, rax		; Save squared remainder
	shrd	rax, rdx, cl		; Shift squared_upper
debug	shr	rdx, cl
assert	jz
	mul	facinv2
	imul	rdx, fac2reg		; Compute quotient * factor
	sub	rem2, rdx		; Subtract from squared remainder

	mov	rax, rem3		; Load remainder
	imul	rax			; Square remainder
	mov	rem3, rax		; Save squared remainder
	shrd	rax, rdx, cl		; Shift squared_upper
debug	shr	rdx, cl
assert	jz
	mul	facinv3
	imul	rdx, fac3reg		; Compute quotient * factor
	sub	rem3, rdx		; Subtract from squared remainder

; At this point, since quotient can be off by 2, the remainder is between
; zero and 3 * factor.  Test for end of squaring loop.  Loop if doubling is
; not required.

	add	temp1, temp1		; One squaring completed, shift
	jnc	sqlp60			; Loop if mul by 2 not needed
	jz	short exit60		; Jump if we are done squaring

; Subtract factor so that remainder is between -1 * factor and 2 * factor.
; Multiply by two, making the remainder between -2 * factor and 4 * factor.

	sub	rem1, fac1reg		; Subtract factor
	sub	rem2, fac2reg
	sub	rem3, fac3reg
	add	rem1, rem1		; Multiply remainder by 2
	add	rem2, rem2
	add	rem3, rem3
	jmp	sqlp60			; Do next iteration

; Make remainder between 0 and factor.  Multiply remainder by two one last
; time (for the last carry out of shifter).  If result = 1 mod factor, then
; we found a divisor of 2^p - 1

exit60:	mov	rax, rem1		; Make rem1 between 0 and 2*fac
	sub	rax, fac1reg
	cmovns	rem1, rax
	mov	rcx, rem2		; Make rem2 between 0 and 2*fac
	sub	rcx, fac2reg
	cmovns	rem2, rcx
	mov	rdx, rem3		; Make rem3 between 0 and 2*fac
	sub	rdx, fac3reg
	cmovns	rem3, rdx
	sub	rax, fac1reg		; Make rem1 between 0 and fac
	cmovns	rem1, rax
	sub	rcx, fac2reg		; Make rem2 between 0 and fac
	cmovns	rem2, rcx
	sub	rdx, fac3reg		; Make rem3 between 0 and fac
	cmovns	rem3, rdx

debug	cmp	rem1, fac1reg		; Remainder should be < factor
assert	jb
debug	cmp	rem2, fac2reg		; Remainder should be < factor
assert	jb
debug	cmp	rem3, fac3reg		; Remainder should be < factor
assert	jb

	add	rem1, rem1		; Double the remainder
	add	rem2, rem2		; Double the remainder
	add	rem3, rem3		; Double the remainder

	sub	rem1, fac1reg		; If rem == fac+1, we found a factor!
	dec	rem1
	jz	short win60
	sub	rem2, fac2reg		; If rem == fac+1, we found a factor!
	dec	rem2
	jz	short win60a
	sub	rem3, fac3reg		; If rem == fac+1, we found a factor!
	dec	rem3
	jz	short win60b

	sub	rdi, rdi		; Clear queued factors count
	mov	rbp, SAVED_REG5		; Restore sieve testing register
	mov	rcx, SAVED_REG4
	mov	rbx, SAVED_REG3
	mov	rax, SAVED_REG2
	mov	rsi, SAVED_REG1
	jmp	bsf60			; Test next factor from sieve

win60:	mov	rax, fac1reg
	jmp	short win60c
win60a:	mov	rax, fac2reg
	jmp	short win60c
win60b:	mov	rax, fac3reg
win60c:	mov	FACLSW, eax		; Factor found!!! Return TRUE
	shr	rax, 32
	mov	FACMSW, eax
	mov	rax, 1
	jmp	done

;***********************************************************************
; Test 61 bit factors
;***********************************************************************

tlp61:	mov	rsi, sieve		; Sieve address
	lea	rbp, [rsi+sievesize]	; Sieve end address
	mov	rbx, 1			; Compute maximum trial factor
	shl	rbx, 61
	xor	rdi, rdi		; Count of queued factors to be tested
svlp61:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf61:	bsf	rdx, rax		; Look for a set bit
	jnz	short test61		; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	cmp	rsi, rbp		; End of sieve?
	jl	short svlp61		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem61a		; Yes, test them

; Save state

	mov	savefac1, rcx		; Save for the restart or more sieving

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle 61-bit overflow

oflow61:and	rdi, rdi		; Are there untested factors?
	jnz	short rem61		; Yes, go do remaining trial factors
	mov	FACMSW, 20000000h	; Return end point
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem61a:	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
					; it will be re-added when we branch back to bsf61
rem61:	mov	r8, fac1		; Copy first trial factor
	mov	fac1[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	short rem61		; No, go copy another
	jmp	short do61		; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test61:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	add	r8, facdists[rdx*8]	; Determine factor to test
	cmp	r8, rbx			; Test for 61-bit overflow
	jae	short oflow61		; Jump if 61-bit overflow
	mov	fac1[rdi*8], r8		; Save the factor to test
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	bsf61			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

do61:	mov	SAVED_REG1, rsi		; Save sieve testing registers
	mov	SAVED_REG2, rax
	mov	SAVED_REG3, rbx
	mov	SAVED_REG4, rcx
	mov	SAVED_REG5, rbp

; Precompute shifted 1 / factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac1reg, fac1
	div	fac1reg			; Compute 1 / trial_factor
	mov	facinv1, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac2reg, fac2
	div	fac2reg			; Compute 1 / trial_factor
	mov	facinv2, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac3reg, fac3
	div	fac3reg			; Compute 1 / trial_factor
	mov	facinv3, rax		; Save 1 / trial_factor

;
; Perform a division on the initial value to get started.
;

	mov	rcx, initshift		; Load initial shift count

	mov	rem1, initval1		; Load initial value
	mov	rem2, rem1		; Load initial value
	mov	rem3, rem1		; Load initial value

	mov	rax, facinv1		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac1reg			; Compute quotient * factor
	sub	rem1, rax		; Subtract from initval
debug	mov	temp1, initval0
debug	sbb	temp1, rdx
assert	jz

	mov	rax, facinv2		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac2reg			; Compute quotient * factor
	sub	rem2, rax		; Subtract from initval
debug	mov	temp2, initval0
debug	sbb	temp2, rdx
assert	jz

	mov	rax, facinv3		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac3reg			; Compute quotient * factor
	sub	rem3, rax		; Subtract from initval
debug	mov	temp3, initval0
debug	sbb	temp3, rdx
assert	jz

;
; Square remainder and get new remainder.  Input remainder is between
; -2 * factor and 2 * factor.  The squared remainder is up to 4 * factor^2.
; The quotient is up to 4 * factor.  For the quotient to fit in 63 bits (the
; maximum possible with our multiply by reciprocal approach), the
; maximum factor is 2^63 / 4, or 2^61.
;

	mov	temp1, shifter		; Load shifter
;	mov	rcx, shift_count	; Load squared value shift count
;	inc	rcx			; Compute an extra quotient bit
					; for better branch prediction
mov rcx,64
sub rcx,shift_count
sqlp61:	mov	rax, rem1		; Load remainder
	imul	rax			; Square remainder
	mov	rem1, rax		; Save squared remainder
	shrd	rax, rdx, cl		; Shift squared_upper
debug	shr	rdx, cl
assert	jz
	mul	facinv1
	imul	rdx, fac1reg		; Compute quotient * factor
	sub	rem1, rdx		; Subtract from squared remainder

	mov	rax, rem2		; Load remainder
	imul	rax			; Square remainder
	mov	rem2, rax		; Save squared remainder
	shrd	rax, rdx, cl		; Shift squared_upper
debug	shr	rdx, cl
assert	jz
	mul	facinv2
	imul	rdx, fac2reg		; Compute quotient * factor
	sub	rem2, rdx		; Subtract from squared remainder

	mov	rax, rem3		; Load remainder
	imul	rax			; Square remainder
	mov	rem3, rax		; Save squared remainder
	shrd	rax, rdx, cl		; Shift squared_upper
debug	shr	rdx, cl
assert	jz
	mul	facinv3
	imul	rdx, fac3reg		; Compute quotient * factor
	sub	rem3, rdx		; Subtract from squared remainder

; At this point, since quotient can be off by 2, the remainder is between
; zero and 3 * factor.  Make remainder between 0 and 2 * factor.

	mov	rax, rem1		; Make rem1 between 0 and 2*fac
	sub	rax, fac1reg
	cmovns	rem1, rax
	mov	rdx, rem2		; Make rem2 between 0 and 2*fac
	sub	rdx, fac2reg
	cmovns	rem2, rdx
	mov	rax, rem3		; Make rem3 between 0 and 2*fac
	sub	rax, fac3reg
	cmovns	rem3, rax

; Test for end of squaring loop.  Loop if doubling is not required.

	add	temp1, temp1		; One squaring completed, shift
	jnc	sqlp61			; Loop if mul by 2 not needed
	jz	short exit61		; Jump if we are done squaring

; Subtract factor so that remainder is between -1 * factor and 1 * factor.
; Multiply by two, making the remainder between -2 * factor and 2 * factor.

	sub	rem1, fac1reg		; Subtract factor
	sub	rem2, fac2reg
	sub	rem3, fac3reg
	add	rem1, rem1		; Multiply remainder by 2
	add	rem2, rem2
	add	rem3, rem3
	jmp	sqlp61			; Do next iteration

; Make remainder between 0 and factor.  Multiply remainder by two one last
; time (for the last carry out of shifter).  If result = 1 mod factor, then
; we found a divisor of 2^p - 1

exit61:	mov	rax, rem1		; Make rem1 between 0 and fac
	sub	rax, fac1reg
	cmovns	rem1, rax
	mov	rcx, rem2		; Make rem2 between 0 and fac
	sub	rcx, fac2reg
	cmovns	rem2, rcx
	mov	rdx, rem3		; Make rem3 between 0 and fac
	sub	rdx, fac3reg
	cmovns	rem3, rdx

debug	cmp	rem1, fac1reg		; Remainder should be < factor
assert	jb
debug	cmp	rem2, fac2reg		; Remainder should be < factor
assert	jb
debug	cmp	rem3, fac3reg		; Remainder should be < factor
assert	jb

	add	rem1, rem1		; Double the remainder
	add	rem2, rem2		; Double the remainder
	add	rem3, rem3		; Double the remainder

	sub	rem1, fac1reg		; If rem == fac+1, we found a factor!
	dec	rem1
	jz	short win61
	sub	rem2, fac2reg		; If rem == fac+1, we found a factor!
	dec	rem2
	jz	short win61a
	sub	rem3, fac3reg		; If rem == fac+1, we found a factor!
	dec	rem3
	jz	short win61b

	sub	rdi, rdi		; Clear queued factors count
	mov	rbp, SAVED_REG5		; Restore sieve testing register
	mov	rcx, SAVED_REG4
	mov	rbx, SAVED_REG3
	mov	rax, SAVED_REG2
	mov	rsi, SAVED_REG1
	jmp	bsf61			; Test next factor from sieve

win61:	mov	rax, fac1reg
	jmp	short win61c
win61a:	mov	rax, fac2reg
	jmp	short win61c
win61b:	mov	rax, fac3reg
win61c:	mov	FACLSW, eax		; Factor found!!! Return TRUE
	shr	rax, 32
	mov	FACMSW, eax
	mov	rax, 1
	jmp	done

;***********************************************************************
; Test 62 bit factors
;***********************************************************************

tlp62:	mov	rsi, sieve		; Sieve address
	lea	rbp, [rsi+sievesize]	; Sieve end address
	mov	rbx, 1			; Compute maximum trial factor
	shl	rbx, 62
	xor	rdi, rdi		; Count of queued factors to be tested
svlp62:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf62:	bsf	rdx, rax		; Look for a set bit
	jnz	short test62		; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	cmp	rsi, rbp		; End of sieve?
	jl	short svlp62		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem62a		; Yes, test them

; Save state

	mov	savefac1, rcx		; Save for the restart or more sieving

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle 62-bit overflow

oflow62:and	rdi, rdi		; Are there untested factors?
	jnz	short rem62		; Yes, go do remaining trial factors
	mov	FACMSW, 40000000h	; Return end point
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem62a:	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
					; it will be re-added when we branch back to bsf62
rem62:	mov	r8, fac1		; Copy first trial factor
	mov	fac1[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	short rem62		; No, go copy another
	jmp	short do62		; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test62:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	add	r8, facdists[rdx*8]	; Determine factor to test
	cmp	r8, rbx			; Test for 62-bit overflow
	jae	short oflow62		; Jump if 62-bit overflow
	mov	fac1[rdi*8], r8		; Save the factor to test
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	bsf62			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

do62:	mov	SAVED_REG1, rsi		; Save sieve testing registers
	mov	SAVED_REG2, rax
	mov	SAVED_REG3, rbx
	mov	SAVED_REG4, rcx
	mov	SAVED_REG5, rbp

; Precompute shifted 1 / factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac1reg, fac1
	div	fac1reg			; Compute 1 / trial_factor
	mov	facinv1, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac2reg, fac2
	div	fac2reg			; Compute 1 / trial_factor
	mov	facinv2, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac3reg, fac3
	div	fac3reg			; Compute 1 / trial_factor
	mov	facinv3, rax		; Save 1 / trial_factor

;
; Perform a division on the initial value to get started.
;

	mov	rcx, initshift		; Load initial shift count

	mov	rem1, initval1		; Load initial value
	mov	rem2, rem1		; Load initial value
	mov	rem3, rem1		; Load initial value

	mov	rax, facinv1		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac1reg			; Compute quotient * factor
	sub	rem1, rax		; Subtract from initval
debug	mov	temp1, initval0
debug	sbb	temp1, rdx
assert	jz

	mov	rax, facinv2		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac2reg			; Compute quotient * factor
	sub	rem2, rax		; Subtract from initval
debug	mov	temp2, initval0
debug	sbb	temp2, rdx
assert	jz

	mov	rax, facinv3		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac3reg			; Compute quotient * factor
	sub	rem3, rax		; Subtract from initval
debug	mov	temp3, initval0
debug	sbb	temp3, rdx
assert	jz

;
; Square remainder and get new remainder.  Input remainder is between
; -1 * factor and 1 * factor.  The squared remainder is up to 1 * factor^2.
; The quotient is up to 1 * factor.  For the quotient to fit in 63 bits (the
; maximum possible with our multiply by reciprocal approach), the
; maximum factor is 2^63 / 1, or 2^63.  So that is not a limitation.
; Instead we are limited by the fact that quotient can be off by 2, giving
; us problems fitting the remainder in 64 bits.  The maximum factor is
; 2^64 / 3 = 2^62.
;

	mov	temp1, shifter		; Load shifter
sqlp62:	mov	rax, rem1		; Load remainder
	imul	rax			; Square remainder
	mov	rem1, rax		; Save squared remainder
	shrd	rax, rdx, 61		; Shift squared_upper
debug	shr	rdx, 61
assert	jz
	mul	facinv1			; Compute quotient
	imul	rdx, fac1reg 		; Compute quotient * factor
	sub	rem1, rdx		; Subtract from squared remainder

	mov	rax, rem2		; Load remainder
	imul	rax			; Square remainder
	mov	rem2, rax		; Save squared remainder
	shrd	rax, rdx, 61		; Shift squared_upper
debug	shr	rdx, 61
assert	jz
	mul	facinv2			; Compute quotient
	imul	rdx, fac2reg 		; Compute quotient * factor
	sub	rem2, rdx		; Subtract from squared remainder

	mov	rax, rem3		; Load remainder
	imul	rax			; Square remainder
	mov	rem3, rax		; Save squared remainder
	shrd	rax, rdx, 61		; Shift squared_upper
debug	shr	rdx, 61
assert	jz
	mul	facinv3			; Compute quotient
	imul	rdx, fac3reg 		; Compute quotient * factor
	sub	rem3, rdx		; Subtract from squared remainder

; At this point, since quotient can be off by 2, the remainder is between
; zero and 3 * factor.  Make remainder between 0 and factor.

	mov	rax, rem1		; Make rem1 between 0 and 2*fac
	sub	rax, fac1reg
	cmovns	rem1, rax
	sub	rax, fac1reg		; Make rem1 between 0 and fac
	cmovns	rem1, rax
	mov	rdx, rem2		; Make rem2 between 0 and 2*fac
	sub	rdx, fac2reg
	cmovns	rem2, rdx
	sub	rdx, fac2reg		; Make rem2 between 0 and fac
	cmovns	rem2, rdx
	mov	rax, rem3		; Make rem3 between 0 and 2*fac
	sub	rax, fac3reg
	cmovns	rem3, rax
	sub	rax, fac3reg		; Make rem3 between 0 and fac
	cmovns	rem3, rax

; Test for end of squaring loop.  Loop if doubling is not required.

	add	temp1, temp1		; One squaring completed, shift
	jnc	sqlp62			; Loop if mul by 2 not needed
	jz	short exit62		; Jump if we are done squaring

; Multiply by two, making the remainder between 0 and 2 * factor.
; Subtract factor so that remainder is between -1 * factor and 1 * factor.

	add	rem1, rem1		; Multiply remainder by 2
	add	rem2, rem2
	add	rem3, rem3
	sub	rem1, fac1reg		; Subtract factor
	sub	rem2, fac2reg
	sub	rem3, fac3reg
	jmp	sqlp62			; Do next iteration

; Multiply remainder by two one last time (for the last carry out of shifter).
; If result = 1 mod factor, then we found a divisor of 2^p - 1

exit62:

debug	cmp	rem1, fac1reg		; Remainder should be < factor
assert	jb
debug	cmp	rem2, fac2reg		; Remainder should be < factor
assert	jb
debug	cmp	rem3, fac3reg		; Remainder should be < factor
assert	jb

	add	rem1, rem1		; Double the remainder
	add	rem2, rem2		; Double the remainder
	add	rem3, rem3		; Double the remainder

	sub	rem1, fac1reg		; If rem == fac+1, we found a factor!
	dec	rem1
	jz	short win62
	sub	rem2, fac2reg		; If rem == fac+1, we found a factor!
	dec	rem2
	jz	short win62a
	sub	rem3, fac3reg		; If rem == fac+1, we found a factor!
	dec	rem3
	jz	short win62b

	sub	rdi, rdi		; Clear queued factors count
	mov	rbp, SAVED_REG5		; Restore sieve testing register
	mov	rcx, SAVED_REG4
	mov	rbx, SAVED_REG3
	mov	rax, SAVED_REG2
	mov	rsi, SAVED_REG1
	jmp	bsf62			; Test next factor from sieve

win62:	mov	rax, fac1reg
	jmp	short win62c
win62a:	mov	rax, fac2reg
	jmp	short win62c
win62b:	mov	rax, fac3reg
win62c:	mov	FACLSW, eax		; Factor found!!! Return TRUE
	shr	rax, 32
	mov	FACMSW, eax
	mov	rax, 1
	jmp	done

;***********************************************************************
; Test 63 bit factors
;***********************************************************************

tlp63:	mov	rsi, sieve		; Sieve address
	lea	rbp, [rsi+sievesize]	; Sieve end address
	xor	rdi, rdi		; Count of queued factors to be tested
svlp63:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf63:	bsf	rdx, rax		; Look for a set bit
	jnz	short test63		; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	cmp	rsi, rbp		; End of sieve?
	jl	short svlp63		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem63a		; Yes, test them

; Save state

	mov	savefac1, rcx		; Save for the restart or more sieving

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle 63-bit overflow

oflow63:and	rdi, rdi		; Are there untested factors?
	jnz	short rem63		; Yes, go do remaining trial factors
	mov	FACMSW, 80000000h	; Return end point
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem63a:	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
					; it will be re-added when we branch back to bsf63
rem63:	mov	r8, fac1		; Copy first trial factor
	mov	fac1[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	short rem63		; No, go copy another
	jmp	short do63		; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test63:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	add	r8, facdists[rdx*8]	; Determine factor to test
	js	short oflow63		; Jump if 63-bit overflow
	mov	fac1[rdi*8], r8		; Save the factor to test
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	bsf63			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

do63:	mov	SAVED_REG1, rsi		; Save sieve testing registers
	mov	SAVED_REG2, rax
	mov	SAVED_REG3, rcx
	mov	SAVED_REG4, rbp

; Precompute shifted 1 / factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac1reg, fac1
	div	fac1reg			; Compute 1 / trial_factor
	mov	facinv1, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac2reg, fac2
	div	fac2reg			; Compute 1 / trial_factor
	mov	facinv2, rax		; Save 1 / trial_factor

	mov	rdx, initdiv0
	sub	rax, rax
	mov	fac3reg, fac3
	div	fac3reg			; Compute 1 / trial_factor
	mov	facinv3, rax		; Save 1 / trial_factor

;
; Perform a division on the initial value to get started.
;

	mov	rcx, initshift		; Load initial shift count

	mov	rem1, initval1		; Load initial value
	mov	rem2, rem1		; Load initial value
	mov	rem3, rem1		; Load initial value

	mov	rax, facinv1		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac1reg			; Compute quotient * factor
	sub	rem1, rax		; Subtract from initval
debug	mov	temp1, initval0
debug	sbb	temp1, rdx
assert	jz

	mov	rax, facinv2		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac2reg			; Compute quotient * factor
	sub	rem2, rax		; Subtract from initval
debug	mov	temp2, initval0
debug	sbb	temp2, rdx
assert	jz

	mov	rax, facinv3		; Quotient = initmul0 * 1 / factor
	shr	rax, cl
	mul	fac3reg			; Compute quotient * factor
	sub	rem3, rax		; Subtract from initval
debug	mov	temp3, initval0
debug	sbb	temp3, rdx
assert	jz

;
; Square remainder and get new remainder.  Input remainder is between
; -1 * factor and 1 * factor.  The squared remainder is up to 1 * factor^2.
; The quotient is up to 1 * factor.  For the quotient to fit in 63 bits (the
; maximum possible with our multiply by reciprocal approach), the
; maximum factor is 2^63 / 1, or 2^63.  So that is not a limitation.
; Instead we are limited by the fact that quotient can be off by 2, giving
; us problems fitting the remainder in 64 bits.  The maximum factor is
; 2^64 / 3.  So unfortunately we must calculate full double word remainders.
;

	mov	rcx, shifter		; Load shifter
sqlp63:	mov	rax, rem1		; Load remainder
	imul	rax			; Square remainder
	mov	rem1, rax		; Save squared remainder
	mov	temp1, rdx
	shrd	rax, rdx, 62		; Shift squared_upper
debug	shr	rdx, 62
assert	jz
	mul	facinv1			; Compute quotient
	mov	rax, fac1reg
	mul	rdx	 		; Compute quotient * factor
	sub	rem1, rax		; Subtract from squared remainder
	sbb	temp1, rdx

	mov	rax, rem2		; Load remainder
	imul	rax			; Square remainder
	mov	rem2, rax		; Save squared remainder
	mov	temp2, rdx
	shrd	rax, rdx, 62		; Shift squared_upper
debug	shr	rdx, 62
assert	jz
	mul	facinv2			; Compute quotient
	mov	rax, fac2reg
	mul	rdx	 		; Compute quotient * factor
	sub	rem2, rax		; Subtract from squared remainder
	sbb	temp2, rdx

	mov	rax, rem3		; Load remainder
	imul	rax			; Square remainder
	mov	rem3, rax		; Save squared remainder
	mov	temp3, rdx
	shrd	rax, rdx, 62		; Shift squared_upper
debug	shr	rdx, 62
assert	jz
	mul	facinv3			; Compute quotient
	mov	rax, fac3reg
	mul	rdx	 		; Compute quotient * factor
	sub	rem3, rax		; Subtract from squared remainder
	sbb	temp3, rdx

; At this point, since quotient can be off by 2, the remainder is between
; zero and 3 * factor.  Make remainder between 0 and factor.

	mov	rax, rem1		; Make rem1 between 0 and 2*fac
	sub	rax, fac1reg
	sbb	temp1, 0
	cmovns	rem1, rax
	sub	rax, fac1reg		; Make rem1 between 0 and fac
	sbb	temp1, 0
	cmovns	rem1, rax
	mov	rdx, rem2		; Make rem2 between 0 and 2*fac
	sub	rdx, fac2reg
	sbb	temp2, 0
	cmovns	rem2, rdx
	sub	rdx, fac2reg		; Make rem2 between 0 and fac
	sbb	temp2, 0
	cmovns	rem2, rdx
	mov	rax, rem3		; Make rem3 between 0 and 2*fac
	sub	rax, fac3reg
	sbb	temp3, 0
	cmovns	rem3, rax
	sub	rax, fac3reg		; Make rem3 between 0 and fac
	sbb	temp3, 0
	cmovns	rem3, rax

; Test for end of squaring loop.  Loop if doubling is not required.

	add	rcx, rcx		; One squaring completed, shift
	jnc	sqlp63			; Loop if mul by 2 not needed
	jz	short exit63		; Jump if we are done squaring

; Multiply by two, making the remainder between 0 and 2 * factor.
; Subtract factor so that remainder is between -1 * factor and 1 * factor.

	add	rem1, rem1		; Multiply remainder by 2
	add	rem2, rem2
	add	rem3, rem3
	sub	rem1, fac1reg		; Subtract factor
	sub	rem2, fac2reg
	sub	rem3, fac3reg
	jmp	sqlp63			; Do next iteration

; Multiply remainder by two one last time (for the last carry out of shifter).
; If result = 1 mod factor, then we found a divisor of 2^p - 1

exit63:

debug	cmp	rem1, fac1reg		; Remainder should be < factor
assert	jb
debug	cmp	rem2, fac2reg		; Remainder should be < factor
assert	jb
debug	cmp	rem3, fac3reg		; Remainder should be < factor
assert	jb

	add	rem1, rem1		; Double the remainder
	add	rem2, rem2		; Double the remainder
	add	rem3, rem3		; Double the remainder

	sub	rem1, fac1reg		; If rem == fac+1, we found a factor!
	dec	rem1
	jz	short win63
	sub	rem2, fac2reg		; If rem == fac+1, we found a factor!
	dec	rem2
	jz	short win63a
	sub	rem3, fac3reg		; If rem == fac+1, we found a factor!
	dec	rem3
	jz	short win63b

	sub	rdi, rdi		; Clear queued factors count
	mov	rbp, SAVED_REG4		; Restore sieve testing register
	mov	rcx, SAVED_REG3
	mov	rax, SAVED_REG2
	mov	rsi, SAVED_REG1
	jmp	bsf63			; Test next factor from sieve

win63:	mov	rax, fac1reg
	jmp	short win63c
win63a:	mov	rax, fac2reg
	jmp	short win63c
win63b:	mov	rax, fac3reg
win63c:	mov	FACLSW, eax		; Factor found!!! Return TRUE
	shr	rax, 32
	mov	FACMSW, eax
	mov	rax, 1
	jmp	done

;***********************************************************************
; For 64-bit factors
;***********************************************************************

;
; Check all the bits in the sieve looking for a factor to test
;

tlp64:	mov	rsi, sieve		; Sieve address
	lea	rbp, [rsi+sievesize]	; Sieve end address
	xor	rdi, rdi		; Count of queued factors to be tested
	mov	rdx, 1			; Compute 1 / factor for initial
	shl	rdx, 63			; estimate in Newton's method
	sub	rax, rax
	div	rcx
	mov	facinv3, rax
svlp64:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf64:	bsf	rdx, rax		; Look for a set bit
	jnz	test64			; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	cmp	rsi, rbp		; End of sieve?
	jl	short svlp64		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem64a		; Yes, test them

; Save state

	mov	savefac1, rcx		; Save for the restart or more sieving

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle 64-bit overflow

oflow64:and	rdi, rdi		; Are there untested factors?
	jnz	short rem64		; Yes, go do remaining trial factors
	mov	FACMSW, 0		; Return end point
	mov	FACHSW, 1
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem64a:	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
					; it will be re-added when we branch back to bsf64
rem64:	mov	r8, fac1		; Copy first trial factor
	mov	fac1[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	short rem64		; No, go copy another
	jmp	short do64		; Go test queued factors

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test64:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	add	r8, facdists[rdx*8]	; Determine factor to test
	jns	short oflow64		; Jump if 64-bit overflow
	mov	fac1[rdi*8], r8		; Save the factor to test
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	bsf64			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

do64:	mov	SAVED_REG1, rsi		; Save sieve testing registers
	mov	SAVED_REG2, rax
	mov	SAVED_REG3, rcx
	mov	SAVED_REG4, rbp

; Precompute 1 / factor using Newton's method and the last 1 / factor
; as our initial estimate.  I think this gives us about 45-ish digits
; of precision which is not enough for our needs.  For example, we miss
; the factor 18586270485220505521 of M530539861 when multithreading.
; Note we round up several of the calculations so that the computed
; 1 / factor is less than the actual 1 / factor.

	mov	rax, facinv3		; Load last 1 / factor
	mul	rax
	lea	temp1, [rdx+1]		; Save ceiling (est^2)

	mov	facinv1, facinv3	; Copy last 1 / factor
	mov	facinv2, facinv3

	mov	rax, fac1
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv1, rdx		; est - f * est^2

	mov	rax, fac2
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv2, rdx		; est - f * est^2

	mov	rax, fac3
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv3, rdx		; est - f * est^2

	add	facinv1, facinv1	; 1/f = 2 (est - f * est^2)
	shl	facinv2, 1		; 1/f = 2 (est - f * est^2)
	add	facinv3, facinv3	; 1/f = 2 (est - f * est^2)

; Refine our 1 / factor estimates using another iteration of Newton's method.
; Note we round up several of the calculations so that the computed 1 / factor
; is less than the actual 1 / factor.

	mov	rax, facinv1		; Load last 1 / factor
	mul	rax
	lea	temp1, [rdx+1]		; Save ceiling (est^2)
	mov	rax, facinv2		; Load last 1 / factor
	mul	rax
	lea	temp2, [rdx+1]		; Save ceiling (est^2)
	mov	rax, facinv3		; Load last 1 / factor
	mul	rax
	lea	temp3, [rdx+1]		; Save ceiling (est^2)

	mov	rax, fac1
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv1, rdx		; est - f * est^2

	mov	rax, fac2
	mul	temp2
	inc	rdx			; ceiling (f * est^2)
	sub	facinv2, rdx		; est - f * est^2

	mov	rax, fac3
	mul	temp3
	inc	rdx			; ceiling (f * est^2)
	sub	facinv3, rdx		; est - f * est^2

	add	facinv1, facinv1	; 1/f = 2 (est - f * est^2)
	shl	facinv2, 1		; 1/f = 2 (est - f * est^2)
	add	facinv3, facinv3	; 1/f = 2 (est - f * est^2)

;
; Perform a division on the initial value to get started.
;

	mov	rcx, initshift		; Load initial shift count

debug	mov	rem1lo, initval1
debug	and	rem1lo, rem1lo
assert	jz

	sub	rem1lo, rem1lo		; Load initial value
	sub	rem2lo, rem2lo
	sub	rem3lo, rem3lo

; Assuming the inverse factor is accurate to 40 bits, we must jump to
; the code that computes the quotient in two parts if the shift count
; is less than 24.

	cmp	cl, 24			; Will shift count yield < 40-bit quot?
	jge	short fast64		; Yes, do the fast first iteration

; This is the slow code to compute the first quotient and remainder.

	mov	temp1, facinv1		; Quotient = val * 1 / factor
	shr	temp1, cl
	mov	temp2, facinv2		; Quotient = val * 1 / factor
	shr	temp2, cl
	mov	temp3, facinv3		; Quotient = val * 1 / factor
	shr	temp3, cl

	mov	rem1hi, initval0	; Load high part of initval
	mov	rem2hi, rem1hi
	mov	rem3hi, rem1hi

	mov	rcx, shifter		; Load shifter
	shr	rcx, 1			; First iteration does not double
	jmp	slow64			; Jump into middle of squaring loop

; This is the fast code to compute the quotients and remainders

fast64:	mov	rax, facinv1		; Quotient = val * 1 / factor
	shr	rax, cl
	mul	fac1			; Compute quotient * factor
	sub	rem1lo, rax		; Subtract from initval
debug	mov	rem1hi, initval0
debug	sbb	rem1hi, rdx
assert	jz

	mov	rax, facinv2		; Quotient = val * 1 / factor
	shr	rax, cl
	mul	fac2			; Compute quotient * factor
	sub	rem2lo, rax		; Subtract from initval
debug	mov	rem2hi, initval0
debug	sbb	rem2hi, rdx
assert	jz

	mov	rax, facinv3		; Quotient = val * 1 / factor
	shr	rax, cl
	mul	fac3			; Compute quotient * factor
	sub	rem3lo, rax		; Subtract from initval
debug	mov	rem3hi, initval0
debug	sbb	rem3hi, rdx
assert	jz

;
; Square remainders and get new remainders.  Input remainder must be less than
; factor.
;

	mov	rcx, shifter		; Load shifter
sqlp64:	mov	rax, rem1lo		; Load remainder
	mul	rem1lo			; Square remainder
	mov	rem1lo, rax		; Save squared remainder
	mov	rem1hi, rdx

	mov	rax, rem2lo		; Load remainder
	mul	rem2lo			; Square remainder
	mov	rem2lo, rax		; Save squared remainder
	mov	rem2hi, rdx

	mov	rax, rem3lo		; Load remainder
	mul	rem3lo			; Square remainder
	mov	rem3lo, rax		; Save squared remainder
	mov	rem3hi, rdx

	mov	rax, rem1hi		; Quotient = val * 1 / factor
	mul	facinv1
	lea	temp1, [rdx+rdx]	; Save quotient

	mov	rax, rem2hi		; Quotient = val * 1 / factor
	mul	facinv2
	lea	temp2, [rdx+rdx]	; Save quotient

	mov	rax, rem3hi		; Quotient = val * 1 / factor
	mul	facinv3
	lea	temp3, [rdx+rdx]	; Save quotient

slow64:	mov	rax, fac1
	mul	temp1			; Compute quotient * factor
	sub	rem1lo, rax		; Subtract from squared remainder
	sbb	rem1hi, rdx
assert	jns

	mov	rax, fac2
	mul	temp2			; Compute quotient * factor
	sub	rem2lo, rax		; Subtract from squared remainder
	sbb	rem2hi, rdx
assert	jns

	mov	rax, fac3
	mul	temp3			; Compute quotient * factor
	sub	rem3lo, rax		; Subtract from squared remainder
	sbb	rem3hi, rdx
assert	jns

	mov	temp1, rem1hi		; Shift for computing accurate quotient
	shld	temp1, rem1lo, 31
	mov	temp2, rem2hi
	shld	temp2, rem2lo, 31
	mov	temp3, rem3hi
	shld	temp3, rem3lo, 31

	mov	rax, temp1		; Quotient = val * 1 / factor
	mul	facinv1
	shr	rdx, 30			; Adjust quotient
	mov	rax, fac1
	mul	rdx			; Compute quotient * factor
	sub	rem1lo, rax		; Subtract from squared remainder
	sbb	rem1hi, rdx
	jz	short qok64		; Most of the time branch will be taken
assert	jns
	sub	rem1lo, fac1		; Handle rare quotient too small case
debug	sbb	rem1hi, 0
assert	jz
qok64:

	mov	rax, temp2		; Quotient = val * 1 / factor
	mul	facinv2
	shr	rdx, 30			; Adjust quotient
	mov	rax, fac2
	mul	rdx			; Compute quotient * factor
	sub	rem2lo, rax		; Subtract from squared remainder
	sbb	rem2hi, rdx
	jz	short qok64a		; Most of the time branch will be taken
assert	jns
	sub	rem2lo, fac2		; Handle rare quotient too small case
debug	sbb	rem2hi, 0
assert	jz
qok64a:

	mov	rax, temp3		; Quotient = val * 1 / factor
	mul	facinv3
	shr	rdx, 30			; Adjust quotient
	mov	rax, fac3
	mul	rdx			; Compute quotient * factor
	sub	rem3lo, rax		; Subtract from squared remainder
	sbb	rem3hi, rdx
	jz	short qok64b		; Most of the time branch will be taken
assert	jns
	sub	rem3lo, fac3		; Handle rare quotient too small case
debug	sbb	rem3hi, 0
assert	jz
qok64b:

; Remainder might still be too large due to the rare quotient too small case.
; If so fix the remainder.

	cmp	rem1lo, fac1		; Fix rem1 if necessary
	jb	ok64
	sub	rem1lo, fac1
ok64:	cmp	rem2lo, fac2		; Fix rem2 if necessary
	jb	ok64a
	sub	rem2lo, fac2
ok64a:	cmp	rem3lo, fac3		; Fix rem3 if necessary
	jb	ok64b
	sub	rem3lo, fac3
ok64b:

;
; Multiply by two if necessary, test for end of squaring loop
;

	add	rcx, rcx		; One squaring completed, shift
	jnc	sqlp64			; Loop if a mul by 2 is not needed
	jz	short exit64		; Are we done squaring?

	lea	temp1, [rem1lo+rem1lo]	; Multiply remainder by 2
	sub	temp1, fac1
	add	rem1lo, rem1lo
	cmovc	rem1lo, temp1		; If mul by 2 overflowed, use rem*2-fac
	cmp	rem1lo, fac1
	cmovae	rem1lo, temp1		; or if rem > fac, use rem*2-fac

	lea	temp2, [rem2lo+rem2lo]	; Multiply remainder by 2
	sub	temp2, fac2
	add	rem2lo, rem2lo	
	cmovc	rem2lo, temp2		; If mul by 2 overflowed, use rem*2-fac
	cmp	rem2lo, fac2
	cmovae	rem2lo, temp2		; or if rem > fac, use rem*2-fac

	lea	temp3, [rem3lo+rem3lo]	; Multiply remainder by 2
	sub	temp3, fac3
	add	rem3lo, rem3lo
	cmovc	rem3lo, temp3		; If mul by 2 overflowed, use rem*2-fac
	cmp	rem3lo, fac3
	cmovae	rem3lo, temp3		; or if rem > fac, use rem*2-fac

	jmp	sqlp64			; Do next squaring

; Multiply remainder by two one last time (for the last carry out of shifter)
; If result = 1 mod factor, then we found a divisor of 2**p - 1

exit64:

debug	cmp	rem1, fac1		; Remainder should be < factor
assert	jb
debug	cmp	rem2, fac2		; Remainder should be < factor
assert	jb
debug	cmp	rem3, fac3		; Remainder should be < factor
assert	jb

	add	rem1lo, rem1lo		; Double the remainder
	add	rem2lo, rem2lo		; Double the remainder
	add	rem3lo, rem3lo		; Double the remainder

	mov	rax, fac1		; Load factor #1
	sub	rem1lo, rax		; Subtract factor
	dec	rem1lo
	jz	short win64
	mov	rax, fac2		; Load factor #2
	sub	rem2lo, rax		; Subtract factor
	dec	rem2lo
	jz	short win64
	mov	rax, fac3		; Load factor #3
	sub	rem3lo, rax		; Subtract factor
	dec	rem3lo
	jz	short win64

	sub	rdi, rdi		; Clear queued factors count
	mov	rbp, SAVED_REG4		; Restore sieve testing register
	mov	rcx, SAVED_REG3
	mov	rax, SAVED_REG2
	mov	rsi, SAVED_REG1
	jmp	bsf64			; Test next factor from sieve

win64:	mov	FACLSW, eax		; Factor found!
	shr	rax, 32
	mov	FACMSW, eax
	mov	rax, 1			; Return TRUE
	jmp	done

;***********************************************************************
; For 65-bit factors
;***********************************************************************

;
; Check all the bits in the sieve looking for a factor to test
;

tlp65:	mov	rsi, sieve		; Sieve address
	lea	rbp, [rsi+sievesize]	; Sieve end address
	xor	rdi, rdi		; Count of queued factors to be tested
	mov	rdx, 1			; Compute 1 / factor for initial
	shl	rdx, 63			; estimate in Newton's method
	sub	rax, rax
	ror	rcx, 1			; Form high 64-bits of factor
	inc	rcx			; ceiling (high 64-bits of factor)
	div	rcx
	dec	rcx
	rol	rcx, 1			; Restore lower 64-bits of factor
	mov	facinv3, rax
svlp65:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf65:	bsf	rdx, rax		; Look for a set bit
	jnz	test65			; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	jc	short oflow65a		; Jump if overflow of 65 bits
	cmp	rsi, rbp		; End of sieve?
	jl	short svlp65		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem65a		; Yes, test them

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle 65-bit overflow

oflow65a:mov	rcx, -1			; Reset rcx so that it will generate a carry if we test
					; factors already in the queue and branch to bsf65
oflow65:and	rdi, rdi		; Are there untested factors?
	jnz	short rem65		; Yes, go do remaining trial factors
	mov	FACMSW, 0		; Return end point
	mov	FACHSW, 2
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem65a:	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
					; it will be re-added when we branch back to bsf65
rem65:	mov	r8, fac1		; Copy first trial factor
	mov	fac1[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	short rem65		; No, go copy another
	jmp	short do65		; Go test queued factors

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test65:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	add	r8, facdists[rdx*8]	; Determine factor to test
	jc	short oflow65		; Jump if 65-bit overflow
	mov	fac1[rdi*8], r8		; Save the factor to test
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	bsf65			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

do65:	mov	SAVED_REG1, rsi		; Save sieve testing registers
	mov	SAVED_REG2, rax
	mov	SAVED_REG3, rcx
	mov	SAVED_REG4, rbp

; Precompute 1 / factor using Newton's method and the last 1 / factor
; as our initial estimate.  I think this gives us about 45-ish digits
; of precision which is not enough for our needs.  For example, we miss
; the factor 18586270485220505521 of M530539861 when multithreading.
; Note we round up several of the calculations so that the computed
; 1 / factor is less than the actual 1 / factor.

	mov	rax, facinv3		; Load last 1 / factor
	mul	rax
	lea	temp1, [rdx+1]		; Save ceiling (est^2)

	mov	facinv1, facinv3	; Copy last 1 / factor
	mov	facinv2, facinv3

	mov	rax, fac1
	ror	rax, 1			; Form upper 64 bits of factor
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv1, rdx		; est - f * est^2

	mov	rax, fac2
	ror	rax, 1			; Form upper 64 bits of factor
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv2, rdx		; est - f * est^2

	mov	rax, fac3
	ror	rax, 1			; Form upper 64 bits of factor
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv3, rdx		; est - f * est^2

	add	facinv1, facinv1	; 1/f = 2 (est - f * est^2)
	shl	facinv2, 1		; 1/f = 2 (est - f * est^2)
	add	facinv3, facinv3	; 1/f = 2 (est - f * est^2)

; Refine our 1 / factor estimates using another iteration of Newton's method.
; Note we round up several of the calculations so that the computed 1 / factor
; is less than the actual 1 / factor.

	mov	rax, facinv1		; Load last 1 / factor
	mul	rax
	lea	temp1, [rdx+1]		; Save ceiling (est^2)
	mov	rax, facinv2		; Load last 1 / factor
	mul	rax
	lea	temp2, [rdx+1]		; Save ceiling (est^2)
	mov	rax, facinv3		; Load last 1 / factor
	mul	rax
	lea	temp3, [rdx+1]		; Save ceiling (est^2)

	mov	rax, fac1
	ror	rax, 1			; Form upper 64 bits of factor
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv1, rdx		; est - f * est^2

	mov	rax, fac2
	ror	rax, 1			; Form upper 64 bits of factor
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp2
	inc	rdx			; ceiling (f * est^2)
	sub	facinv2, rdx		; est - f * est^2

	mov	rax, fac3
	ror	rax, 1			; Form upper 64 bits of factor
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp3
	inc	rdx			; ceiling (f * est^2)
	sub	facinv3, rdx		; est - f * est^2

	add	facinv1, facinv1	; 1/f = 2 (est - f * est^2)
	shl	facinv2, 1		; 1/f = 2 (est - f * est^2)
	add	facinv3, facinv3	; 1/f = 2 (est - f * est^2)

;
; Perform a division on the initial value to get started.
;

	mov	rcx, initshift		; Load initial shifted count

	and	rcx, rcx		; If quotient will be small,
	jns	short fast65		; jump to go compute it.

	sub	rem1hi, rem1hi		; Load initial value
	sub	rem1lo, rem1lo
	sub	rem2hi, rem2hi
	sub	rem2lo, rem2lo
	sub	rem3hi, rem3hi
	sub	rem3lo, rem3lo

	mov	rcx, initshift2		; Compute first quotients
	mov	temp1, facinv1
	shr	temp1, cl
	mov	temp2, facinv2
	shr	temp2, cl
	mov	temp3, facinv3
	shr	temp3, cl

	mov	rcx, shifter		; Load shifter
	shr	rcx, 1
	jmp	slow65			; Jump into middle of squaring loop

fast65:	mov	rem1hi, initval1	; Load initial value
	sub	rem1lo, rem1lo
	mov	rem2hi, rem1hi
	sub	rem2lo, rem2lo
	mov	rem3hi, rem1hi
	sub	rem3lo, rem3lo

	mov	rax, facinv1		; Quotient = val * 1 / factor
	shr	rax, cl
	sub	rem1hi, rax		; Subtract quotient * factor_hi
	mul	fac1			; Compute quotient * factor_lo
	sub	rem1lo, rax		; Subtract from initval
	sbb	rem1hi, rdx
assert	jns

	mov	rax, facinv2		; Quotient = val * 1 / factor
	shr	rax, cl
	sub	rem2hi, rax		; Subtract quotient * factor_hi
	mul	fac2			; Compute quotient * factor_lo
	sub	rem2lo, rax		; Subtract from initval
	sbb	rem2hi, rdx
assert	jns

	mov	rax, facinv3		; Quotient = val * 1 / factor
	shr	rax, cl
	sub	rem3hi, rax		; Subtract quotient * factor_hi
	mul	fac3			; Compute quotient * factor_lo
	sub	rem3lo, rax		; Subtract from initval
	sbb	rem3hi, rdx
assert	jns

;
; Square remainders and get new remainders
;

	mov	rcx, shifter		; Load shifter
sqlp65:	mov	temp1, rem1hi		; Load remainder_hi
	imul	temp1, temp1		; Square remainder_hi

	mov	temp2, rem2hi		; Load remainder_hi
	imul	temp2, temp2		; Square remainder_hi

	mov	temp3, rem3hi		; Load remainder_hi
	imul	temp3, temp3		; Square remainder_hi

	add	rem1hi, rem1hi		; Double remainder_hi
	add	rem2hi, rem2hi		; Double remainder_hi
	add	rem3hi, rem3hi		; Double remainder_hi

	mov	rax, rem1lo		; Load remainder_lo
	mul	rem1lo			; Square remainder_lo
	xchg	rem1lo, rax		; Save squared remainder_lo
	xchg	rem1hi, rdx
	mul	rdx			; Calc 2 * remainder_lo * remainder_hi
	add	rem1hi, rax		; Add in dbl'ed mul'ed remainder
	adc	temp1, rdx

	mov	rax, rem2lo		; Load remainder_lo
	mul	rem2lo			; Square remainder_lo
	xchg	rem2lo, rax		; Save squared remainder_lo
	xchg	rem2hi, rdx
	mul	rdx			; Calc 2 * remainder_lo * remainder_hi
	add	rem2hi, rax		; Add in dbl'ed mul'ed remainder
	adc	temp2, rdx

	mov	rax, rem3lo		; Load remainder_lo
	mul	rem3lo			; Square remainder_lo
	xchg	rem3lo, rax		; Save squared remainder_lo
	xchg	rem3hi, rdx
	mul	rdx			; Calc 2 * remainder_lo * remainder_hi
	add	rem3hi, rax		; Add in dbl'ed mul'ed remainder
	adc	temp3, rdx

	shld	temp1, rem1hi, 31	; Shift for computing quotient
	shld	temp2, rem2hi, 31
	shld	temp3, rem3hi, 31

	mov	rax, facinv1		; Quotient = val * 1 / factor
	mul	temp1
	mov	temp1, rdx		; Save quotient

	mov	rax, facinv2		; Quotient = val * 1 / factor
	mul	temp2
	mov	temp2, rdx		; Save quotient

	mov	rax, facinv3		; Quotient = val * 1 / factor
	mul	temp3
	mov	temp3, rdx		; Save quotient

slow65:	mov	rax, fac1
	mul	temp1			; Compute quotient * factor_lo
	add	rdx, temp1		; Add in quotient * factor_hi
	shld	rdx, rax, 33
	shl	rax, 33
	sub	rem1lo, rax		; Subtract from squared remainder
	sbb	rem1hi, rdx
assert	jns

	mov	rax, fac2
	mul	temp2			; Compute quotient * factor_lo
	add	rdx, temp2		; Add in quotient * factor_hi
	shld	rdx, rax, 33
	shl	rax, 33
	sub	rem2lo, rax		; Subtract from squared remainder
	sbb	rem2hi, rdx
assert	jns

	mov	rax, fac3
	mul	temp3			; Compute quotient * factor_lo
	add	rdx, temp3		; Add in quotient * factor_hi
	shld	rdx, rax, 33
	shl	rax, 33
	sub	rem3lo, rax		; Subtract from squared remainder
	sbb	rem3hi, rdx
assert	jns

	mov	rax, rem1hi		; Quotient = val * 1 / factor
	mul	facinv1
	sub	rem1hi, rdx		; Subtract quotient * factor_hi
	mov	rax, fac1
	mul	rdx			; Compute quotient * factor_lo
	sub	rem1lo, rax		; Subtract from squared remainder
	sbb	rem1hi, rdx
assert	jns

	mov	rax, rem2hi		; Quotient = val * 1 / factor
	mul	facinv2
	sub	rem2hi, rdx		; Subtract quotient * factor_hi
	mov	rax, fac2
	mul	rdx			; Compute quotient * factor_lo
	sub	rem2lo, rax		; Subtract from squared remainder
	sbb	rem2hi, rdx
assert	jns

	mov	rax, rem3hi		; Quotient = val * 1 / factor
	mul	facinv3
	sub	rem3hi, rdx		; Subtract quotient * factor_hi
	mov	rax, fac3
	mul	rdx			; Compute quotient * factor_lo
	sub	rem3lo, rax		; Subtract from squared remainder
	sbb	rem3hi, rdx
assert	jns

;
; Multiply by two if necessary, test for end of squaring loop
;

	add	rcx, rcx		; One squaring completed, shift
	jnc	sqlp65			; Loop if a mul by 2 is not needed
	jz	short exit65		; Are we done squaring?
	add	rem1lo, rem1lo		; Multiply remainder by 2
	adc	rem1hi, rem1hi
	add	rem2lo, rem2lo		; Multiply remainder by 2
	adc	rem2hi, rem2hi
	add	rem3lo, rem3lo		; Multiply remainder by 2
	adc	rem3hi, rem3hi
	jmp	sqlp65			; Do next squaring

;
; Multiply remainder by two one last time (for the last carry out of shifter)
; If result = 1 mod factor, then we found a divisor of 2**p - 1
;

exit65:	add	rem1lo, rem1lo		; Double the remainder
	adc	rem1hi, rem1hi
	add	rem2lo, rem2lo		; Double the remainder
	adc	rem2hi, rem2hi
	add	rem3lo, rem3lo		; Double the remainder
	adc	rem3hi, rem3hi
	mov	rax, fac1		; Load factor #1
sub65:	sub	rem1lo, rax		; Subtract factor
	sbb	rem1hi, 1
	jns	short sub65		; Subtract until result is negative
	add	rem1lo, rax		; Make it positive
	adc	rem1hi, 1
	dec	rem1lo			; If one, its a factor
	jnz	short fail65
	and	rem1hi, rem1hi
	jz	short win65
fail65:	mov	rax, fac2		; Load factor #2
sub65a:	sub	rem2lo, rax		; Subtract factor
	sbb	rem2hi, 1
	jns	short sub65a		; Subtract until result is negative
	add	rem2lo, rax		; Make it positive
	adc	rem2hi, 1
	dec	rem2lo			; If one, its a factor
	jnz	short fail65a
	and	rem2hi, rem2hi
	jz	short win65
fail65a:mov	rax, fac3		; Load factor #3
sub65b:	sub	rem3lo, rax		; Subtract factor
	sbb	rem3hi, 1
	jns	short sub65b		; Subtract until result is negative
	add	rem3lo, rax		; Make it positive
	adc	rem3hi, 1
	dec	rem3lo			; If one, its a factor
	jnz	short next65
	and	rem3hi, rem3hi
	jz	short win65

next65:	sub	rdi, rdi		; Clear queued factors count
	mov	rbp, SAVED_REG4		; Restore sieve testing register
	mov	rcx, SAVED_REG3
	mov	rax, SAVED_REG2
	mov	rsi, SAVED_REG1
	jmp	bsf65			; Test next factor from sieve

win65:	mov	FACLSW, eax		; Factor found!
	shr	rax, 32
	mov	FACMSW, eax
	mov	FACHSW, 1
	mov	rax, 1			; Return TRUE
	jmp	done


;********************************************
; For 66-bit factors and above (up to 2^79.5)
;********************************************

;
; Check all the bits in the sieve looking for a factor to test
;

tlp66:	mov	rsi, sieve		; Sieve address
	lea	rbp, [rsi+sievesize]	; Sieve end address
	xor	rdi, rdi		; Count of queued factors to be tested
	mov	rdx, 1			; Compute 1 / factor for initial
	shl	rdx, 63			; estimate in Newton's method
	sub	rax, rax
	mov	r8, rcx			; Save factor
	mov	rcx, shift_count
	dec	rcx
	mov	r14, 1			; Form trial factor limit
	ror	r14, cl
	shld	rbx, r8, cl		; Form high 64-bits of factor
	inc	rbx			; ceiling (high 64-bits of factor)
	div	rbx
	dec	rbx
	shr	rbx, cl			; Restore factor
	mov	rcx, r8
	mov	facinv3, rax
svlp66:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf66:	bsf	rdx, rax		; Look for a set bit
	jnz	test66			; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	adc	rbx, 0
	cmp	rsi, rbp		; End of sieve?
	jl	short svlp66		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem66a		; Yes, test them

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	FACHSW, ebx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle overflow

oflow66:and	rdi, rdi		; Are there untested factors?
	jnz	short rem66		; Yes, go do remaining trial factors
	mov	FACHSW, r14d		; Return end point
	mov	FACMSW, 0
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem66a:	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
	sbb	rbx, 0			; it will be re-added when we branch back to bsf66
rem66:	mov	r8, fac1		; Copy first trial factor
	mov	r9, fac1hi
	mov	fac1[rdi*8], r8
	mov	fac1hi[rdi*8], r9
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	short rem66		; No, go copy another
	jmp	short do66		; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test66:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	mov	r9, rbx
	add	r8, facdists[rdx*8]	; Determine factor to test
	adc	r9, 0
	cmp	r9, r14			; Test for overflow
	jae	short oflow66		; Jump if overflow
	mov	fac1[rdi*8], r8		; Save the factor to test
	mov	fac1hi[rdi*8], r9
	inc	rdi			; One more factor queued up
	cmp	rdi, 3			; Have enough been queued up?
	jne	bsf66			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

do66:	mov	SAVED_REG1, rsi		; Save sieve testing registers
	mov	SAVED_REG2, rax
	mov	SAVED_REG3, rbx
	mov	SAVED_REG4, rcx
	mov	SAVED_REG5, rbp
	mov	SAVED_REG6, r14

; Precompute 1 / factor using Newton's method and the last 1 / factor
; as our initial estimate.  I think this gives us about 45-ish digits
; of precision which is not enough for our needs.  For example, we run
; into problems with "Factor=990000793,65,66" and a large number of threads.
; Note we round up several of the calculations so that the computed
; 1 / factor is less than the actual 1 / factor.

	mov	rax, facinv3		; Load last 1 / factor
	mul	rax
	lea	temp1, [rdx+1]		; Save ceiling (est^2)

	mov	facinv1, facinv3	; Copy last 1 / factor
	mov	facinv2, facinv3

	mov	rcx, shift_count
	dec	rcx

	mov	rax, fac1hi		; Form upper 64 bits of factor
	mov	rdx, fac1
	shld	rax, rdx, cl
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv1, rdx		; est - f * est^2

	mov	rax, fac2hi		; Form upper 64 bits of factor
	mov	rdx, fac2
	shld	rax, rdx, cl
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv2, rdx		; est - f * est^2

	mov	rax, fac3hi		; Form upper 64 bits of factor
	mov	rdx, fac3
	shld	rax, rdx, cl
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv3, rdx		; est - f * est^2

	add	facinv1, facinv1	; 1/f = 2 (est - f * est^2)
	shl	facinv2, 1		; 1/f = 2 (est - f * est^2)
	add	facinv3, facinv3	; 1/f = 2 (est - f * est^2)

; Refine our 1 / factor estimates using another iteration of Newton's method.
; Note we round up several of the calculations so that the computed 1 / factor
; is less than the actual 1 / factor.

	mov	rax, facinv1		; Load last 1 / factor
	mul	rax
	lea	temp1, [rdx+1]		; Save ceiling (est^2)
	mov	rax, facinv2		; Load last 1 / factor
	mul	rax
	lea	temp2, [rdx+1]		; Save ceiling (est^2)
	mov	rax, facinv3		; Load last 1 / factor
	mul	rax
	lea	temp3, [rdx+1]		; Save ceiling (est^2)

	mov	rax, fac1hi		; Form upper 64 bits of factor
	mov	rdx, fac1
	shld	rax, rdx, cl
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp1
	inc	rdx			; ceiling (f * est^2)
	sub	facinv1, rdx		; est - f * est^2

	mov	rax, fac2hi		; Form upper 64 bits of factor
	mov	rdx, fac2
	shld	rax, rdx, cl
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp2
	inc	rdx			; ceiling (f * est^2)
	sub	facinv2, rdx		; est - f * est^2

	mov	rax, fac3hi		; Form upper 64 bits of factor
	mov	rdx, fac3
	shld	rax, rdx, cl
	inc	rax			; ceiling (upper 64 bits of factor)
	mul	temp3
	inc	rdx			; ceiling (f * est^2)
	sub	facinv3, rdx		; est - f * est^2

	add	facinv1, facinv1	; 1/f = 2 (est - f * est^2)
	shl	facinv2, 1		; 1/f = 2 (est - f * est^2)
	add	facinv3, facinv3	; 1/f = 2 (est - f * est^2)

;
; Perform a division on the initial value to get started.
;

	mov	rcx, initshift		; Load initial shifted count

	and	rcx, rcx		; If quotient will be small,
	jns	short fast66		; jump to go compute it.

	sub	rem1hi, rem1hi		; Load initial value
	sub	rem1lo, rem1lo
	sub	rem2hi, rem2hi
	sub	rem2lo, rem2lo
	sub	rem3hi, rem3hi
	sub	rem3lo, rem3lo

	mov	rcx, initshift2		; Compute first quotients
	mov	temp1, facinv1
	shr	temp1, cl
	mov	temp2, facinv2
	shr	temp2, cl
	mov	temp3, facinv3
	shr	temp3, cl

	mov	rcx, shifter		; Load shifter
	shr	rcx, 1
	mov	memshifter, rcx		; Save shifter
	mov	rcx, shift66
	jmp	slow66			; Jump into middle of squaring loop

fast66:	mov	rem1hi, initval1	; Load initial value
	sub	rem1lo, rem1lo
	mov	rem2hi, rem1hi
	sub	rem2lo, rem2lo
	mov	rem3hi, rem1hi
	sub	rem3lo, rem3lo

	mov	rax, facinv1		; Quotient = val * 1 / factor
	shr	rax, cl
	mov	temp1, fac1hi		; Subtract quotient * factor_hi
	imul	temp1, rax
	sub	rem1hi, temp1
	mul	fac1			; Compute quotient * factor_lo
	sub	rem1lo, rax		; Subtract from initval
	sbb	rem1hi, rdx
assert	jns

	mov	rax, facinv2		; Quotient = val * 1 / factor
	shr	rax, cl
	mov	temp2, fac2hi		; Subtract quotient * factor_hi
	imul	temp2, rax
	sub	rem2hi, temp2
	mul	fac2			; Compute quotient * factor_lo
	sub	rem2lo, rax		; Subtract from initval
	sbb	rem2hi, rdx
assert	jns

	mov	rax, facinv3		; Quotient = val * 1 / factor
	shr	rax, cl
	mov	temp3, fac3hi		; Subtract quotient * factor_hi
	imul	temp3, rax
	sub	rem3hi, temp3
	mul	fac3			; Compute quotient * factor_lo
	sub	rem3lo, rax		; Subtract from initval
	sbb	rem3hi, rdx
assert	jns

;
; Square remainders and get new remainders
;

	mov	rcx, shifter		; Load shifter
	mov	memshifter, rcx		; Save shifter
	mov	rcx, shift66
sqlp66:	mov	temp1, rem1hi		; Load remainder_hi
	imul	temp1, temp1		; Square remainder_hi

	mov	temp2, rem2hi		; Load remainder_hi
	imul	temp2, temp2		; Square remainder_hi

	mov	temp3, rem3hi		; Load remainder_hi
	imul	temp3, temp3		; Square remainder_hi

	add	rem1hi, rem1hi		; Double remainder_hi
	add	rem2hi, rem2hi		; Double remainder_hi
	add	rem3hi, rem3hi		; Double remainder_hi

	mov	rax, rem1lo		; Load remainder_lo
	mul	rem1lo			; Square remainder_lo
	xchg	rem1lo, rax		; Save squared remainder_lo
	xchg	rem1hi, rdx
	mul	rdx			; Calc 2 * remainder_lo * remainder_hi
	add	rem1hi, rax		; Add in dbl'ed mul'ed remainder
	adc	temp1, rdx

	mov	rax, rem2lo		; Load remainder_lo
	mul	rem2lo			; Square remainder_lo
	xchg	rem2lo, rax		; Save squared remainder_lo
	xchg	rem2hi, rdx
	mul	rdx			; Calc 2 * remainder_lo * remainder_hi
	add	rem2hi, rax		; Add in dbl'ed mul'ed remainder
	adc	temp2, rdx

	mov	rax, rem3lo		; Load remainder_lo
	mul	rem3lo			; Square remainder_lo
	xchg	rem3lo, rax		; Save squared remainder_lo
	xchg	rem3hi, rdx
	mul	rdx			; Calc 2 * remainder_lo * remainder_hi
	add	rem3hi, rax		; Add in dbl'ed mul'ed remainder
	adc	temp3, rdx

	shld	temp1, rem1hi, 31	; Shift for computing quotient
	shld	temp2, rem2hi, 31
	shld	temp3, rem3hi, 31

	mov	rax, facinv1		; Quotient = val * 1 / factor
	mul	temp1
	mov	temp1, rdx		; Save quotient

	mov	rax, facinv2		; Quotient = val * 1 / factor
	mul	temp2
	mov	temp2, rdx		; Save quotient

	mov	rax, facinv3		; Quotient = val * 1 / factor
	mul	temp3
	mov	temp3, rdx		; Save quotient

slow66:	ror	rcx, 8			; Use shift count in ch

	mov	rax, fac1
	mul	temp1			; Compute quotient * factor_lo
	imul	temp1, fac1hi		; Compute quotient * factor_hi
	add	rdx, temp1		; Add in quotient * factor_hi
	shld	rdx, rax, cl
	shl	rax, cl
	sub	rem1lo, rax		; Subtract from squared remainder
	sbb	rem1hi, rdx
assert	jns

	mov	rax, fac2
	mul	temp2			; Compute quotient * factor_lo
	imul	temp2, fac2hi		; Compute quotient * factor_hi
	add	rdx, temp2		; Add in quotient * factor_hi
	shld	rdx, rax, cl
	shl	rax, cl
	sub	rem2lo, rax		; Subtract from squared remainder
	sbb	rem2hi, rdx
assert	jns

	mov	rax, fac3
	mul	temp3			; Compute quotient * factor_lo
	imul	temp3, fac3hi		; Compute quotient * factor_hi
	add	rdx, temp3		; Add in quotient * factor_hi
	shld	rdx, rax, cl
	shl	rax, cl
	sub	rem3lo, rax		; Subtract from squared remainder
	sbb	rem3hi, rdx
assert	jns

	rol	rcx, 8			; Restore shift count

	mov	rax, rem1hi		; Quotient = val * 1 / factor
	mul	facinv1
	shr	rdx, cl
	mov	temp1, fac1hi		; Compute quotient * factor_hi
	imul	temp1, rdx
	sub	rem1hi, temp1		; Subtract quotient * factor_hi
	mov	rax, fac1
	mul	rdx			; Compute quotient * factor_lo
	sub	rem1lo, rax		; Subtract from squared remainder
	sbb	rem1hi, rdx
assert	jns

	mov	rax, rem2hi		; Quotient = val * 1 / factor
	mul	facinv2
	shr	rdx, cl
	mov	temp2, fac2hi		; Compute quotient * factor_hi
	imul	temp2, rdx
	sub	rem2hi, temp2		; Subtract quotient * factor_hi
	mov	rax, fac2
	mul	rdx			; Compute quotient * factor_lo
	sub	rem2lo, rax		; Subtract from squared remainder
	sbb	rem2hi, rdx
assert	jns

	mov	rax, rem3hi		; Quotient = val * 1 / factor
	mul	facinv3
	shr	rdx, cl
	mov	temp3, fac3hi		; Compute quotient * factor_hi
	imul	temp3, rdx
	sub	rem3hi, temp3		; Subtract quotient * factor_hi
	mov	rax, fac3
	mul	rdx			; Compute quotient * factor_lo
	sub	rem3lo, rax		; Subtract from squared remainder
	sbb	rem3hi, rdx
assert	jns

;
; Multiply by two if necessary, test for end of squaring loop
;

	mov	rax, memshifter		; One squaring completed, shift
	add	rax, rax
	mov	memshifter, rax
	jnc	sqlp66			; Loop if a mul by 2 is not needed
	jz	short exit66		; Are we done squaring?
	add	rem1lo, rem1lo		; Multiply remainder by 2
	adc	rem1hi, rem1hi
	add	rem2lo, rem2lo		; Multiply remainder by 2
	adc	rem2hi, rem2hi
	add	rem3lo, rem3lo		; Multiply remainder by 2
	adc	rem3hi, rem3hi
	jmp	sqlp66			; Do next squaring

;
; Multiply remainder by two one last time (for the last carry out of shifter)
; If result = 1 mod factor, then we found a divisor of 2**p - 1
;

exit66:	add	rem1lo, rem1lo		; Double the remainder
	adc	rem1hi, rem1hi
	add	rem2lo, rem2lo		; Double the remainder
	adc	rem2hi, rem2hi
	add	rem3lo, rem3lo		; Double the remainder
	adc	rem3hi, rem3hi
	mov	rax, fac1		; Load factor #1
	mov	rdx, fac1hi
sub66:	sub	rem1lo, rax		; Subtract factor
	sbb	rem1hi, rdx
	jns	short sub66		; Subtract until result is negative
	add	rem1lo, rax		; Make it positive
	adc	rem1hi, rdx
	jnz	short fail66
	dec	rem1lo			; If one, its a factor
	jz	short win66
fail66:	mov	rax, fac2		; Load factor #2
	mov	rdx, fac2hi
sub66a:	sub	rem2lo, rax		; Subtract factor
	sbb	rem2hi, rdx
	jns	short sub66a		; Subtract until result is negative
	add	rem2lo, rax		; Make it positive
	adc	rem2hi, rdx
	jnz	short fail66a
	dec	rem2lo			; If one, its a factor
	jz	short win66
fail66a:mov	rax, fac3		; Load factor #3
	mov	rdx, fac3hi
sub66b:	sub	rem3lo, rax		; Subtract factor
	sbb	rem3hi, rdx
	jns	short sub66b		; Subtract until result is negative
	add	rem3lo, rax		; Make it positive
	adc	rem3hi, rdx
	jnz	short next66
	dec	rem3lo			; If one, its a factor
	jz	short win66

next66:	sub	rdi, rdi		; Clear queued factors count
	mov	r14, SAVED_REG6		; Restore sieve testing register
	mov	rbp, SAVED_REG5
	mov	rcx, SAVED_REG4
	mov	rbx, SAVED_REG3
	mov	rax, SAVED_REG2
	mov	rsi, SAVED_REG1
	jmp	bsf66			; Test next factor from sieve

win66:	mov	FACLSW, eax		; Factor found!
	shr	rax, 32
	mov	FACMSW, eax
	mov	FACHSW, edx
	mov	rax, 1			; Return TRUE
	jmp	done

;***************************************************************************
; For up to 74-bit factors using SSE2, may not be faster than the code above
;***************************************************************************

;
; Check all the bits in the sieve looking for a factor to test
;

tlp74:	finit				; Set for 64-bit precision
	mov	rsi, sieve		; Sieve address
	lea	r10, [rsi+sievesize]	; Sieve end address
	xor	rdi, rdi		; Count of queued factors to be tested
	mov	r8, rcx			; Save factor
	mov	rcx, shift_count
dec	rcx
	mov	r14, 1			; Form trial factor limit
	ror	r14, cl
	mov	rcx, r8			; Restore factor
svlp74:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf74:	bsf	rdx, rax		; Look for a set bit
	jnz	test74			; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	adc	rbx, 0
	cmp	rsi, r10		; End of sieve?
	jl	short svlp74		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem74a		; Yes, test them

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	FACHSW, ebx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle overflow

oflow74:and	rdi, rdi		; Are there untested factors?
	jnz	short rem74		; Yes, go do remaining trial factors
	mov	FACHSW, r14d		; Return end point
	mov	FACMSW, 0
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem74a:	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
	sbb	rbx, 0			; it will be re-added when we branch back to bsf74
rem74:	mov	r8, QWORD PTR XMM_F1	; Copy first trial factor
	mov	r9, QWORD PTR XMM_F2
	mov	rdx, QWORD PTR XMM_F3
	mov	QWORD PTR XMM_F1[rdi*8], r8
	mov	QWORD PTR XMM_F2[rdi*8], r9
	mov	QWORD PTR XMM_F3[rdi*8], rdx
	mov	r8, QWORD PTR XMM_INVFAC
	mov	QWORD PTR XMM_INVFAC[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 4			; Have enough been queued up?
	jne	short rem74		; No, go copy another
	jmp	do74			; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test74:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	mov	r9, rbx
	add	r8, facdists[rdx*8]	; Determine factor to test
	adc	r9, 0
	cmp	r9, r14			; Test for overflow
	jae	oflow74			; Jump if overflow

; Compute the factor to test and 63 most significant bits of 1 / factor
;; OPT...  old 32-bit code....  convert to divpd???

	mov	r12, r8			; Save the factor to test
	mov	r13, r9
	shld	r13, r12, 34
	mov	temp, r13
	fild	temp
	fmul	TWO_TO_30
	and	r12, 3FFFFFFFh
	mov	temp, r12
	fiadd	DWORD PTR temp		; We now have the factor to test
	fld	TWO_TO_FACSIZE_PLUS_62	; Constant to generate 63 bit inverse
	fdivrp	st(1), st
	fistp	QWORD PTR XMM_INVFAC[rdi*8]

; Compute the factor in 30 bit chunks
;  OPT - just save r8/r9 and use sse2/avx instr to convert to 30 bit chunks ??? (If better, do it in 32-bit code too)

	shld	r9, r8, 4
	mov	edx, 3FFFFFFFh
	and	rdx, r8
	mov	XMM_F3[rdi*8], edx
	shr	r8, 30
	and	r8, 3FFFFFFFh
	mov	XMM_F2[rdi*8], r8d
	mov	XMM_F1[rdi*8], r9d

; Check to see if we've accumulated enough factors to test

	inc	rdi			; One more factor queued up
	cmp	rdi, 4			; Have enough been queued up?
	jne	bsf74			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

; Work on initval.
; This is like the loop74 code except that we avoid the initial squaring.

do74:	sse2_fac_initval

; Square remainder and get new remainder.

	mov	edi, SSE2_LOOP_COUNTER	; Number of times to loop
loop74:	sse2_fac 74
	dec	edi			; Decrement loop counter
	jnz	loop74
;;	sub	rdi, rdi		; Clear queued factors count

; If result = factor + 1, then we found a divisor of 2**p - 1

	sse2_compare_part1		; Convert LSW of results
	pcmpeqd xmm2, XMMWORD PTR XMM_F3 ; See if any remainder is factor + 1
	pcmpeqd xmm8, XMMWORD PTR XMM_F3a ; See if any remainder is factor + 1
	movdqa	xmm0, xmm2
	por	xmm0, xmm8
	ptest	xmm0, XMMWORD PTR XMM_BITS30 ; See if any least significant words indicate a possible factor
	jz	bsf74			; No, test next factor from sieve

	sse2_compare_part2		; Convert high and middle words of the results
	pcmpeqd xmm1, XMMWORD PTR XMM_F2
	pcmpeqd xmm0, XMMWORD PTR XMM_F1
	pand	xmm2, XMMWORD PTR XMM_BITS30
	pand	xmm2, xmm1
	ptest	xmm2, xmm0
	jnz	short win74_1		; Jump if a factor found
	pcmpeqd xmm7, XMMWORD PTR XMM_F2a
	pcmpeqd xmm6, XMMWORD PTR XMM_F1a
	pand	xmm8, XMMWORD PTR XMM_BITS30
	pand	xmm8, xmm7
	ptest	xmm8, xmm6
	jnz	short win74_2		; Jump if a factor found
	jmp	bsf74			; Test next factor from sieve

win74_1:
	pand	xmm2, xmm0
	pmovmskb rcx, xmm2
	sub	rdx, rdx		; See if first SSE2 word was the factor
	shr	rcx, 1
	jc	short win74		; Yes! Factor found
	lea	rdx, [rdx+8]		; Last SSE2 word must have been the factor
	jmp	short win74		; Yes! Factor found

win74_2:
	pand	xmm8, xmm6
	pmovmskb rcx, xmm8
	mov	edx, 16			; See if first SSE2 word was the factor
	shr	rcx, 1
	jc	short win74		; Yes! Factor found
	lea	rdx, [rdx+8]		; Last SSE2 word must have been the factor
;;	jmp	short win74		; Yes! Factor found

win74:	mov	eax, XMM_F3[rdx]	; Factor found!!! Return it
	mov	ebx, XMM_F2[rdx]
	mov	ecx, XMM_F1[rdx]
	shl	eax, 2
	shrd	eax, ebx, 2
	shl	ebx, 2
	shrd	ebx, ecx, 4
	shr	ecx, 4
	mov	FACLSW, eax
	mov	FACMSW, ebx
	mov	FACHSW, ecx
	mov	rax, 1			; return TRUE
	jmp	done

;***************************************************************************
; For 75 to 86-bit factors using SSE2, may not be faster than the code above
;***************************************************************************

;
; Check all the bits in the sieve looking for a factor to test
;

tlp86:	finit				; Set for 64-bit precision
	mov	rsi, sieve		; Sieve address
	lea	r10, [rsi+sievesize]	; Sieve end address
	xor	rdi, rdi		; Count of queued factors to be tested
	mov	r8, rcx			; Save factor
	mov	rcx, shift_count
dec	rcx
	mov	r14, 1			; Form trial factor limit
	ror	r14, cl
	mov	rcx, r8			; Restore factor
svlp86:	mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
bsf86:	bsf	rdx, rax		; Look for a set bit
	jnz	test86			; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	adc	rbx, 0
	cmp	rsi, r10		; End of sieve?
	jl	short svlp86		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short rem86a		; Yes, test them

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	FACHSW, ebx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle overflow

oflow86:and	rdi, rdi		; Are there untested factors?
	jnz	short rem86		; Yes, go do remaining trial factors
	mov	FACHSW, r14d		; Return end point
	mov	FACMSW, 0
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

rem86a:	sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
	sbb	rbx, 0			; it will be re-added when we branch back to bsf86
rem86:	mov	r8, QWORD PTR XMM_F1	; Copy first trial factor
	mov	r9, QWORD PTR XMM_F2
	mov	rdx, QWORD PTR XMM_F3
	mov	QWORD PTR XMM_F1[rdi*8], r8
	mov	QWORD PTR XMM_F2[rdi*8], r9
	mov	QWORD PTR XMM_F3[rdi*8], rdx
	mov	r8, QWORD PTR XMM_INVFAC
	mov	QWORD PTR XMM_INVFAC[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 4			; Have enough been queued up?
	jne	short rem86		; No, go copy another
	jmp	do86			; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

test86:	btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	mov	r9, rbx
	add	r8, facdists[rdx*8]	; Determine factor to test
	adc	r9, 0
	cmp	r9, r14			; Test for overflow
	jae	oflow86			; Jump if overflow

; Compute the factor to test and 63 most significant bits of 1 / factor
;; OPT...  old 32-bit code....  convert to divpd???

	mov	r12, r8			; Save the factor to test
	mov	r13, r9
	shld	r13, r12, 34
	mov	temp, r13
	fild	temp
	fmul	TWO_TO_30
	and	r12, 3FFFFFFFh
	mov	temp, r12
	fiadd	DWORD PTR temp		; We now have the factor to test
	fld	TWO_TO_FACSIZE_PLUS_62	; Constant to generate 63 bit inverse
	fdivrp	st(1), st
	fistp	QWORD PTR XMM_INVFAC[rdi*8]

; Compute the factor in 30 bit chunks
;  OPT - just save r8/r9 and use sse2/avx instr to convert to 30 bit chunks ??? (If better, do it in 32-bit code too)

	shld	r9, r8, 4
	mov	edx, 3FFFFFFFh
	and	rdx, r8
	mov	XMM_F3[rdi*8], edx
	shr	r8, 30
	and	r8, 3FFFFFFFh
	mov	XMM_F2[rdi*8], r8d
	mov	XMM_F1[rdi*8], r9d

; Check to see if we've accumulated enough factors to test

	inc	rdi			; One more factor queued up
	cmp	rdi, 4			; Have enough been queued up?
	jne	bsf86			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

; Work on initval.
; This is like the loop86 code except that we avoid the initial squaring.

do86:	sse2_fac_initval

; Square remainder and get new remainder.

	mov	edi, SSE2_LOOP_COUNTER	; Number of times to loop
loop86:	sse2_fac 86
	dec	edi			; Decrement loop counter
	jnz	loop86
;;	sub	rdi, rdi		; Clear queued factors count

; If result = factor + 1, then we found a divisor of 2**p - 1

	sse2_compare_part1		; Convert LSW of results
	pcmpeqd xmm2, XMMWORD PTR XMM_F3 ; See if any remainder is factor + 1
	pcmpeqd xmm8, XMMWORD PTR XMM_F3a ; See if any remainder is factor + 1
	movdqa	xmm0, xmm2
	por	xmm0, xmm8
	ptest	xmm0, XMMWORD PTR XMM_BITS30 ; See if any least significant words indicate a possible factor
	jz	bsf86			; No, test next factor from sieve

	sse2_compare_part2		; Convert high and middle words of the results
	pcmpeqd xmm1, XMMWORD PTR XMM_F2
	pcmpeqd xmm0, XMMWORD PTR XMM_F1
	pand	xmm2, XMMWORD PTR XMM_BITS30
	pand	xmm2, xmm1
	ptest	xmm2, xmm0
	jnz	short win86_1		; Jump if a factor found
	pcmpeqd xmm7, XMMWORD PTR XMM_F2a
	pcmpeqd xmm6, XMMWORD PTR XMM_F1a
	pand	xmm8, XMMWORD PTR XMM_BITS30
	pand	xmm8, xmm7
	ptest	xmm8, xmm6
	jnz	short win86_2		; Jump if a factor found
	jmp	bsf86			; Test next factor from sieve

win86_1:
	pand	xmm2, xmm0
	pmovmskb rcx, xmm2
	sub	rdx, rdx		; See if first SSE2 word was the factor
	shr	rcx, 1
	jc	short win86		; Yes! Factor found
	lea	rdx, [rdx+8]		; Last SSE2 word must have been the factor
	jmp	short win86		; Yes! Factor found

win86_2:
	pand	xmm8, xmm6
	pmovmskb rcx, xmm8
	mov	edx, 16			; See if first SSE2 word was the factor
	shr	rcx, 1
	jc	short win86		; Yes! Factor found
	lea	rdx, [rdx+8]		; Last SSE2 word must have been the factor
;;	jmp	short win86		; Yes! Factor found

win86:	mov	eax, XMM_F3[rdx]	; Factor found!!! Return it
	mov	ebx, XMM_F2[rdx]
	mov	ecx, XMM_F1[rdx]
	shl	eax, 2
	shrd	eax, ebx, 2
	shl	ebx, 2
	shrd	ebx, ecx, 4
	shr	ecx, 4
	mov	FACLSW, eax
	mov	FACMSW, ebx
	mov	FACHSW, ecx
	mov	rax, 1			; return TRUE
	jmp	done


;************************************
; For up to 74-bit factors using AVX2
;************************************

;
; Check all the bits in the sieve looking for a factor to test
;

atlp74:	finit				; Set for 64-bit precision
	mov	rsi, sieve		; Sieve address
	lea	r10, [rsi+sievesize]	; Sieve end address
	xor	rdi, rdi		; Count of queued factors to be tested
	mov	r8, rcx			; Save factor
	mov	rcx, shift_count
	dec	rcx
	mov	r14, 1			; Form trial factor limit
	ror	r14, cl
	mov	rcx, r8			; Restore factor
asvlp74:mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
absf74:	bsf	rdx, rax		; Look for a set bit
	jnz	atest74			; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	adc	rbx, 0
	cmp	rsi, r10		; End of sieve?
	jl	short asvlp74		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short arem74a		; Yes, test them

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	FACHSW, ebx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle overflow

aoflow74:and	rdi, rdi		; Are there untested factors?
	jnz	short arem74		; Yes, go do remaining trial factors
	mov	FACHSW, r14d		; Return end point
	mov	FACMSW, 0
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

arem74a:sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
	sbb	rbx, 0			; it will be re-added when we branch back to bsf74
arem74:	mov	r8, QWORD PTR YMM_F1	; Copy first trial factor
	mov	r9, QWORD PTR YMM_F2
	mov	rdx, QWORD PTR YMM_F3
	mov	QWORD PTR YMM_F1[rdi*8], r8
	mov	QWORD PTR YMM_F2[rdi*8], r9
	mov	QWORD PTR YMM_F3[rdi*8], rdx
	mov	r8, QWORD PTR YMM_INVFAC
	mov	QWORD PTR YMM_INVFAC[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 8			; Have enough been queued up?
	jne	short arem74		; No, go copy another
	jmp	ado74			; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

atest74:btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	mov	r9, rbx
	add	r8, facdists[rdx*8]	; Determine factor to test
	adc	r9, 0
	cmp	r9, r14			; Test for overflow
	jae	aoflow74		; Jump if overflow

; Compute the factor to test and 63 most significant bits of 1 / factor
;; OPT...  old 32-bit code....  convert to divpd???

	mov	r12, r8			; Save the factor to test
	mov	r13, r9
	shld	r13, r12, 34
	mov	temp, r13
	fild	temp
	fmul	TWO_TO_30
	and	r12, 3FFFFFFFh
	mov	temp, r12
	fiadd	DWORD PTR temp		; We now have the factor to test
	fld	TWO_TO_FACSIZE_PLUS_62	; Constant to generate 63 bit inverse
	fdivrp	st(1), st
	fistp	QWORD PTR YMM_INVFAC[rdi*8]

; Compute the factor in 30 bit chunks
;  OPT - just save r8/r9 and use sse2/avx instr to convert to 30 bit chunks ??? (If better, do it in 32-bit code too)

	shld	r9, r8, 4
	mov	edx, 3FFFFFFFh
	and	rdx, r8
	mov	YMM_F3[rdi*8], edx
	shr	r8, 30
	and	r8, 3FFFFFFFh
	mov	YMM_F2[rdi*8], r8d
	mov	YMM_F1[rdi*8], r9d

; Check to see if we've accumulated enough factors to test

	inc	rdi			; One more factor queued up
	cmp	rdi, 8			; Have enough been queued up?
	jne	absf74			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

; Work on initval.
; This is like the aloop74 code except that we avoid the initial squaring.

ado74:	avx2_fac_initval

; Square remainder and get new remainder.

	mov	edi, SSE2_LOOP_COUNTER	; Number of times to loop
aloop74:avx2_fac 74
	dec	edi			; Decrement loop counter
	jnz	aloop74
;;	sub	rdi, rdi		; Clear queued factors count

; If result = factor + 1, then we found a divisor of 2**p - 1

	avx2_compare_part1
	vpcmpeqq ymm2, ymm2, YMMWORD PTR YMM_F3	; See if any remainder is factor + 1
	vpcmpeqq ymm8, ymm8, YMMWORD PTR YMM_F3a ; See if any remainder is factor + 1
	vpor	ymm0, ymm2, ymm8
	vptest	ymm0, ymm0		; See if any least significant words indicate a possible factor
	jz	absf74			; No, test next factor from sieve

	avx2_compare_part2
	vpcmpeqq ymm1, ymm1, YMMWORD PTR YMM_F2
	vpcmpeqq ymm0, ymm0, YMMWORD PTR YMM_F1
	vpand	ymm2, ymm2, ymm1
	vptest	ymm2, ymm0
	jnz	short awin74_1		; Jump if a factor found
	vpcmpeqq ymm7, ymm7, YMMWORD PTR YMM_F2a
	vpcmpeqq ymm6, ymm6, YMMWORD PTR YMM_F1a
	vpand	ymm8, ymm8, ymm7
	vptest	ymm8, ymm6
	jnz	short awin74_2		; Jump if a factor found
	jmp	absf74			; Test next factor from sieve

awin74_1:
	vpand	ymm2, ymm2, ymm0
	vpmovmskb rcx, ymm2
	sub	rdx, rdx		; See if first AVX word was the factor
	shr	rcx, 1
	jc	short awin74		; Yes! Factor found
	lea	rdx, [rdx+8]		; See if next AVX word was the factor
	shr	rcx, 8
	jc	short awin74		; Yes! Factor found
	lea	rdx, [rdx+8]		; See if next AVX word was the factor
	shr	rcx, 8
	jc	short awin74		; Yes! Factor found
	lea	rdx, [rdx+8]		; Last AVX word must have been the factor
	jmp	short awin74		; Yes! Factor found

awin74_2:
	vpand	ymm8, ymm8, ymm6
	vpmovmskb rcx, ymm8
	mov	edx, 32			; See if first AVX word was the factor
	shr	rcx, 1
	jc	short awin74		; Yes! Factor found
	lea	rdx, [rdx+8]		; See if next AVX word was the factor
	shr	rcx, 8
	jc	short awin74		; Yes! Factor found
	lea	rdx, [rdx+8]		; See if next AVX word was the factor
	shr	rcx, 8
	jc	short awin74		; Yes! Factor found
	lea	rdx, [rdx+8]		; Last AVX word must have been the factor
;;	jmp	short awin74		; Yes! Factor found

awin74:	mov	eax, YMM_F3[rdx]	; Factor found!!! Return it
	mov	ebx, YMM_F2[rdx]
	mov	ecx, YMM_F1[rdx]
	shl	eax, 2
	shrd	eax, ebx, 2
	shl	ebx, 2
	shrd	ebx, ecx, 4
	shr	ecx, 4
	mov	FACLSW, eax
	mov	FACMSW, ebx
	mov	FACHSW, ecx
	mov	rax, 1			; return TRUE
	jmp	done

;************************************
; For 75 to 86-bit factors using AVX2
;************************************

;
; Check all the bits in the sieve looking for a factor to test
;

atlp86:	finit				; Set for 64-bit precision
	mov	rsi, sieve		; Sieve address
	lea	r10, [rsi+sievesize]	; Sieve end address
	xor	rdi, rdi		; Count of queued factors to be tested
	mov	r8, rcx			; Save factor
	mov	rcx, shift_count
	dec	rcx
	mov	r14, 1			; Form trial factor limit
	ror	r14, cl
	mov	rcx, r8			; Restore factor
asvlp86:mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
absf86:	bsf	rdx, rax		; Look for a set bit
	jnz	atest86			; Found one, go test the factor
	add	rcx, facdist64      	; Add facdist * 64 to the factor
	adc	rbx, 0
	cmp	rsi, r10		; End of sieve?
	jl	short asvlp86		; Loop to test next sieve qword

; Check queued counter

	and	rdi, rdi		; Are there untested factors?
	jnz	short arem86a		; Yes, test them

; Return so caller can check for ESC

	shr	rcx, 32
	mov	FACMSW, ecx
	mov	FACHSW, ebx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Handle overflow

aoflow86:and	rdi, rdi		; Are there untested factors?
	jnz	short arem86		; Yes, go do remaining trial factors
	mov	FACHSW, r14d		; Return end point
	mov	FACMSW, 0
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

arem86a:sub	rcx, facdist64      	; Undo the add of facdist * 64 to rcx
	sbb	rbx, 0			; it will be re-added when we branch back to bsf86
arem86:	mov	r8, QWORD PTR YMM_F1	; Copy first trial factor
	mov	r9, QWORD PTR YMM_F2
	mov	rdx, QWORD PTR YMM_F3
	mov	QWORD PTR YMM_F1[rdi*8], r8
	mov	QWORD PTR YMM_F2[rdi*8], r9
	mov	QWORD PTR YMM_F3[rdi*8], rdx
	mov	r8, QWORD PTR YMM_INVFAC
	mov	QWORD PTR YMM_INVFAC[rdi*8], r8
	inc	rdi			; One more factor queued up
	cmp	rdi, 8			; Have enough been queued up?
	jne	short arem86		; No, go copy another
	jmp	ado86			; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

atest86:btr	rax, rdx		; Clear the sieve bit
	mov	r8, rcx			; Copy base factor
	mov	r9, rbx
	add	r8, facdists[rdx*8]	; Determine factor to test
	adc	r9, 0
	cmp	r9, r14			; Test for overflow
	jae	aoflow86		; Jump if overflow

; Compute the factor to test and 63 most significant bits of 1 / factor
;; OPT...  old 32-bit code....  convert to divpd???

	mov	r12, r8			; Copy the factor to test
	mov	r13, r9
	shld	r13, r12, 34
	mov	temp, r13
	fild	temp
	fmul	TWO_TO_30
	and	r12, 3FFFFFFFh
	mov	temp, r12
	fiadd	DWORD PTR temp		; We now have the factor to test
	fld	TWO_TO_FACSIZE_PLUS_62	; Constant to generate 63 bit inverse
	fdivrp	st(1), st
	fistp	QWORD PTR YMM_INVFAC[rdi*8]

; Compute the factor in 30 bit chunks
;  OPT - just save r8/r9 and use sse2/avx instr to convert to 30 bit chunks ??? (If better, do it in 32-bit code too)

	shld	r9, r8, 4
	mov	edx, 3FFFFFFFh
	and	rdx, r8
	mov	YMM_F3[rdi*8], edx
	shr	r8, 30
	and	r8, 3FFFFFFFh
	mov	YMM_F2[rdi*8], r8d
	mov	YMM_F1[rdi*8], r9d

; Check to see if we've accumulated enough factors to test

	inc	rdi			; One more factor queued up
	cmp	rdi, 8			; Have enough been queued up?
	jne	absf86			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

; Work on initval.
; This is like the aloop86 code except that we avoid the initial squaring.

ado86:	avx2_fac_initval

; Square remainder and get new remainder.

	mov	edi, SSE2_LOOP_COUNTER	; Number of times to loop
aloop86:avx2_fac 86
	dec	edi			; Decrement loop counter
	jnz	aloop86
;;	sub	rdi, rdi		; Clear queued factors count

; If result = factor + 1, then we found a divisor of 2**p - 1

	avx2_compare_part1
	vpcmpeqq ymm2, ymm2, YMMWORD PTR YMM_F3	; See if any remainder is factor + 1
	vpcmpeqq ymm8, ymm8, YMMWORD PTR YMM_F3a ; See if any remainder is factor + 1
	vpor	ymm0, ymm2, ymm8
	vptest	ymm0, ymm0		; See if any least significant words indicate a possible factor
	jz	absf86			; No, test next factor from sieve

	avx2_compare_part2
	vpcmpeqq ymm1, ymm1, YMMWORD PTR YMM_F2
	vpcmpeqq ymm0, ymm0, YMMWORD PTR YMM_F1
	vpand	ymm2, ymm2, ymm1
	vptest	ymm2, ymm0
	jnz	short awin86_1		; Jump if a factor found
	vpcmpeqq ymm7, ymm7, YMMWORD PTR YMM_F2a
	vpcmpeqq ymm6, ymm6, YMMWORD PTR YMM_F1a
	vpand	ymm8, ymm8, ymm7
	vptest	ymm8, ymm6
	jnz	short awin86_2		; Jump if a factor found
	jmp	absf86			; Test next factor from sieve

awin86_1:
	vpand	ymm2, ymm2, ymm0
	vpmovmskb rcx, ymm2
	sub	rdx, rdx		; See if first AVX word was the factor
	shr	rcx, 1
	jc	short awin86		; Yes! Factor found
	lea	rdx, [rdx+8]		; See if next AVX word was the factor
	shr	rcx, 8
	jc	short awin86		; Yes! Factor found
	lea	rdx, [rdx+8]		; See if next AVX word was the factor
	shr	rcx, 8
	jc	short awin86		; Yes! Factor found
	lea	rdx, [rdx+8]		; Last AVX word must have been the factor
	jmp	short awin86		; Yes! Factor found

awin86_2:
	vpand	ymm8, ymm8, ymm6
	vpmovmskb rcx, ymm8
	mov	edx, 32			; See if first AVX word was the factor
	shr	rcx, 1
	jc	short awin86		; Yes! Factor found
	lea	rdx, [rdx+8]		; See if next AVX word was the factor
	shr	rcx, 8
	jc	short awin86		; Yes! Factor found
	lea	rdx, [rdx+8]		; See if next AVX word was the factor
	shr	rcx, 8
	jc	short awin86		; Yes! Factor found
	lea	rdx, [rdx+8]		; Last AVX word must have been the factor
;;	jmp	short awin86		; Yes! Factor found

awin86:	mov	eax, YMM_F3[rdx]	; Factor found!!! Return it
	mov	ebx, YMM_F2[rdx]
	mov	ecx, YMM_F1[rdx]
	shl	eax, 2
	shrd	eax, ebx, 2
	shl	ebx, 2
	shrd	ebx, ecx, 4
	shr	ecx, 4
	mov	FACLSW, eax
	mov	FACMSW, ebx
	mov	FACHSW, ecx
	mov	rax, 1			; return TRUE
	jmp	done


;***************************************************************************
; For 48 to 96-bit factors using AVX and FMA (only tested to 2^92 in C code)
;***************************************************************************

;
; Check all the bits in the sieve looking for a factor to test
;

	QUEUE_COUNT = 12
ftlp96:	mov	rsi, sieve		; Sieve address
	lea	r10, [rsi+sievesize]	; Sieve end address
	mov	rdi, QUEUE_COUNT	; Count of queued factors to be tested
	mov	r12, rcx		; Move factor LSW to r12, we use rcx for shift counts within our loop
	mov	ecx, FMA_64_MINUS_LOW_BITS ; Load constant for shifting trial factors
fsvlp96:mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
fbsf96:	bsf	rdx, rax		; Look for a set bit
	jnz	ftest96			; Found one, go test the factor
	add	r12, facdist64      	; Add facdist * 64 to the factor
	adc	rbx, 0
	cmp	rsi, r10		; End of sieve?
	jl	short fsvlp96		; Loop to test next sieve qword

; Check queued counter

	cmp	rdi, QUEUE_COUNT	; Are there untested factors?
	jne	short frem96a		; Yes, test them

; Return so caller can check for ESC

	shr	r12, 32
	mov	FACMSW, r12d
	mov	FACHSW, ebx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

frem96a:sub	r12, facdist64      	; Undo the add of facdist * 64 to rcx
	sbb	rbx, 0			; it will be re-added when we branch back to bsf96
frem96:	mov	r8, YMM_FMA_FAC_HI_INT[QUEUE_COUNT*8-8]	; Copy first trial factor
	mov	r9, YMM_FMA_FAC_LO_INT[QUEUE_COUNT*8-8]
	mov	YMM_FMA_FAC_HI_INT[rdi*8-8], r8
	mov	YMM_FMA_FAC_LO_INT[rdi*8-8], r9
	dec	rdi			; One more factor queued up
	jnz	short frem96		; No, go copy another
	jmp	fdo96			; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

ftest96:btr	rax, rdx		; Clear the sieve bit
	mov	r8, r12			; Copy base factor
	mov	r9, rbx
	add	r8, facdists[rdx*8]	; Determine factor to test
	adc	r9, 0

; Save the factor as integers.  We start the normalize and convert to doubles here since we do not have AVX2 available.

;mov	r9, FOO1
;mov	r8, FOO2
	shld	r9, r8, cl
	mov	YMM_FMA_FAC_HI_INT[rdi*8-8], r9
	mov	YMM_FMA_FAC_LO_INT[rdi*8-8], r8

; Check to see if we've accumulated enough factors to test

	dec	rdi			; One more factor queued up, queue full?
	jnz	fbsf96			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

; Work on initval.
; This is like the floop96 code except that we avoid the initial squaring.

fdo96:	avx_fma_fac_initval

; Repeatedly square remainder and get new remainder.

	mov	edi, SSE2_LOOP_COUNTER	; Number of times to loop
floop96:avx_fma_fac
	dec	edi			; Decrement loop counter
	jnz	floop96
	mov	rdi, QUEUE_COUNT	; Count of queued factors to be tested

; If result = factor + 1, then we found a divisor of 2**p - 1

	avx_fma_compare
	vcmpeqpd ymm2, ymm2, YMMWORD PTR YMM_FMA_ONE	; See if any remainder is 1
	vcmpeqpd ymm6, ymm6, YMMWORD PTR YMM_FMA_ONE	; See if any remainder is 1
	vcmpeqpd ymm10, ymm10, YMMWORD PTR YMM_FMA_ONE	; See if any remainder is 1
	vorpd	ymm0, ymm2, ymm6
	vorpd	ymm0, ymm0, ymm10
	vptest	ymm0, ymm0
	jz	fbsf96			; If no factor found, test next factor from sieve

	vmovmskpd rcx, ymm10		; Merge masks into one register
	vmovmskpd rdx, ymm6
	shl	rcx, 4
	or	rcx, rdx
	vmovmskpd rdx, ymm2
	shl	rcx, 4
	or	rcx, rdx
	bsf	rdx, rcx		; Look for the set bit

	mov	rax, YMM_FMA_FAC_LO_INT[rdx*8]	; Factor found!!! Return it
	mov	rbx, YMM_FMA_FAC_HI_INT[rdx*8]
	mov	ecx, FMA_64_MINUS_LOW_BITS
	shl	rax, cl
	shrd	rax, rbx, cl
	mov	FACLSW, eax		; Save LSW (low 32 bits)
	shr	rax, 32
	mov	FACMSW, eax		; Save MSW
	shr	rbx, cl
	mov	FACHSW, ebx		; Save HSW
	mov	rax, 1			; return TRUE
	jmp	done


;******************************************************************************
; For 48 to 96-bit factors using AVX512 and FMA (only tested to 2^92 in C code)
;******************************************************************************

;
; Check all the bits in the sieve looking for a factor to test
;

	QUEUE_COUNT = 32
zftlp96:mov	rsi, sieve		; Sieve address
	lea	r10, [rsi+sievesize]	; Sieve end address
	mov	rdi, QUEUE_COUNT	; Count of queued factors to be tested
	mov	r12, rcx		; Move factor LSW to r12, we use rcx for shift counts within our loop
	mov	ecx, FMA_64_MINUS_LOW_BITS ; Load constant for shifting trial factors
zfsvlp96:mov	rax, [rsi]		; Load word from sieve
	lea	rsi, [rsi+8]		; Bump sieve address
zfbsf96:bsf	rdx, rax		; Look for a set bit
	jnz	zftest96		; Found one, go test the factor
	add	r12, facdist64      	; Add facdist * 64 to the factor
	adc	rbx, 0
	cmp	rsi, r10		; End of sieve?
	jl	short zfsvlp96		; Loop to test next sieve qword

; Check queued counter

	cmp	rdi, QUEUE_COUNT	; Are there untested factors?
	jne	short zfrem96a		; Yes, test them

; Return so caller can check for ESC

	shr	r12, 32
	mov	FACMSW, r12d
	mov	FACHSW, ebx
	mov	rax, 2			; Return for ESC check
	jmp	done

; Test less-than-full queue by copying first trial factor in the queue

zfrem96a:sub	r12, facdist64      	; Undo the add of facdist * 64 to rcx
	sbb	rbx, 0			; it will be re-added when we branch back to bsf96
zfrem96:mov	r8, ZMM_FMA_FAC_HI_INT[QUEUE_COUNT*8-8]	; Copy first trial factor
	mov	r9, ZMM_FMA_FAC_LO_INT[QUEUE_COUNT*8-8]
	mov	ZMM_FMA_FAC_HI_INT[rdi*8-8], r8
	mov	ZMM_FMA_FAC_LO_INT[rdi*8-8], r9
	dec	rdi			; One more factor queued up
	jnz	short zfrem96		; No, go copy another
	jmp	zfdo96			; Yes, go test them

;
; Gather trial factors to be tested all at once to minimize
; processor stalls.
;

zftest96:btr	rax, rdx		; Clear the sieve bit
	mov	r8, r12			; Copy base factor
	mov	r9, rbx
	add	r8, facdists[rdx*8]	; Determine factor to test
	adc	r9, 0

;;;BUG - can we do the adds above using AVX2?  carries are a bit of an issue

; Save the factor as integers.  We will normalize and convert to doubles in AVX code.

;mov	r9, FOO1
;mov	r8, FOO2
	shld	r9, r8, cl
	mov	ZMM_FMA_FAC_HI_INT[rdi*8-8], r9
	mov	ZMM_FMA_FAC_LO_INT[rdi*8-8], r8

; Check to see if we've accumulated enough factors to test

	dec	rdi			; One more factor queued up, queue full?
	jnz	zfbsf96			; No, go test more sieve bits

;
; Now test the accumulated trial factors
;

; Work on initval.
; This is like the zfloop96 code except that we avoid the initial squaring.

zfdo96:	avx512_fma_fac_initval

; Repeatedly square remainder and get new remainder.

	mov	edi, SSE2_LOOP_COUNTER	; Number of times to loop
zfloop96:avx512_fma_fac
	dec	edi			; Decrement loop counter
	jnz	zfloop96
	mov	rdi, QUEUE_COUNT	; Count of queued factors to be tested

; If result = factor + 1, then we found a divisor of 2**p - 1

	avx512_fma_compare
	vpcmpeqq k1, zmm2, ZMMWORD PTR ZMM_FMA_ONE	; See if any remainder is 1
	vpcmpeqq k2, zmm9, ZMMWORD PTR ZMM_FMA_ONE	; See if any remainder is 1
	vpcmpeqq k3, zmm16, ZMMWORD PTR ZMM_FMA_ONE	; See if any remainder is 1
	vpcmpeqq k4, zmm23, ZMMWORD PTR ZMM_FMA_ONE	; See if any remainder is 1
	korw	k5, k1, k2
	korw	k5, k5, k3
	kortestw k5, k4
	jz	zfbsf96			; If no factor found, test next factor from sieve

	kmovw	ecx, k4			; Merge the masks into one register
	kmovw	edx, k3
	shl	rcx, 8
	or	rcx, rdx
	kmovw	edx, k2
	shl	rcx, 8
	or	rcx, rdx
	kmovw	edx, k1
	shl	rcx, 8
	or	rcx, rdx
	bsf	rdx, rcx		; Look for the set bit

	mov	rax, ZMM_FMA_FAC_LO_INT[rdx*8]	; Factor found!!! Return it
	mov	rbx, ZMM_FMA_FAC_HI_INT[rdx*8]
	mov	ecx, FMA_64_MINUS_LOW_BITS
	shl	rax, cl
	shrd	rax, rbx, cl
	mov	FACLSW, eax		; Save LSW (low 32 bits)
	shr	rax, 32
	mov	FACMSW, eax		; Save MSW
	shr	rbx, cl
	mov	FACHSW, ebx		; Save HSW
	mov	rax, 1			; return TRUE
	jmp	done


; Pop registers and return

done:	ad_epilog 0,0,rbx,rbp,rsi,rdi,r12,r13,r14,r15,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

factor64_tf ENDP

_TEXT	ENDS
END
