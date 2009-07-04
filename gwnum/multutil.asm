; Copyright 2001-2009 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine used to be part of mult.asm.  However, linux x86-64 objcopy had
; difficulty with "mov rax, OFFSET variable" when variable was declared in the
; same file.
;

	TITLE   multutil

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xmult.mac

VERSION_NUMBER = 2511		;; Version 25.11

;; Make prctabs public.  Linux x86-64 objcopy handles loading offsets of
;; variables only if they are public variables in one file and defined
;; extrn wherever they are used.

IFNDEF X86_64
PUBLIC	prctab1, prctab2
ENDIF
PUBLIC	xprctab1, xprctab2, xprctab2a

;; Routines located in other asm files

IFNDEF X86_64
EXTRN	gwadd1:PROC
EXTRN	gwaddq1:PROC
EXTRN	gwsub1:PROC
EXTRN	gwsubq1:PROC
EXTRN	gwaddsub1:PROC
EXTRN	gwaddsubq1:PROC
EXTRN	gwcopyzero1:PROC
EXTRN	gwmuls1:PROC
EXTRN	gwadd2:PROC
EXTRN	gwaddq2:PROC
EXTRN	gwsub2:PROC
EXTRN	gwsubq2:PROC
EXTRN	gwaddsub2:PROC
EXTRN	gwaddsubq2:PROC
EXTRN	gwcopyzero2:PROC
EXTRN	gwmuls2:PROC
ENDIF

EXTRN	gwxadd1:PROC
EXTRN	gwxaddq1:PROC
EXTRN	gwxsub1:PROC
EXTRN	gwxsubq1:PROC
EXTRN	gwxaddsub1:PROC
EXTRN	gwxaddsubq1:PROC
EXTRN	gwxcopyzero1:PROC
EXTRN	gwxadds1:PROC
EXTRN	gwxmuls1:PROC
EXTRN	gwxadd2:PROC
EXTRN	gwxaddq2:PROC
EXTRN	gwxsub2:PROC
EXTRN	gwxsubq2:PROC
EXTRN	gwxaddsub2:PROC
EXTRN	gwxaddsubq2:PROC
EXTRN	gwxcopyzero2:PROC
EXTRN	gwxadds2:PROC
EXTRN	gwxmuls2:PROC

;; List of possible suffixes for the normalization routines

exnorm	MACRO num, prefix, suffix
	exnorm1	&prefix, &num&z&suffix
	exnorm1	&prefix, &num&ze&suffix

	exnorm1	&prefix, &num&suffix
	exnorm1	&prefix, &num&e&suffix
	exnorm1	&prefix, &num&c&suffix
	exnorm1	&prefix, &num&ec&suffix
	exnorm1	&prefix, &num&zp&suffix
	exnorm1	&prefix, &num&zpe&suffix
	exnorm1	&prefix, &num&zpc&suffix
	exnorm1	&prefix, &num&zpec&suffix
IFNB <prefix> 
	exnorm1	&prefix, &num&b&suffix
	exnorm1	&prefix, &num&eb&suffix
	exnorm1	&prefix, &num&cb&suffix
	exnorm1	&prefix, &num&ecb&suffix
	exnorm1	&prefix, &num&zpb&suffix
	exnorm1	&prefix, &num&zpeb&suffix
	exnorm1	&prefix, &num&zpcb&suffix
	exnorm1	&prefix, &num&zpecb&suffix

	exnorm1	&prefix, &num&s4&suffix
	exnorm1	&prefix, &num&es4&suffix
	exnorm1	&prefix, &num&cs4&suffix
	exnorm1	&prefix, &num&ecs4&suffix
	exnorm1	&prefix, &num&zps4&suffix
	exnorm1	&prefix, &num&zpes4&suffix
	exnorm1	&prefix, &num&zpcs4&suffix
	exnorm1	&prefix, &num&zpecs4&suffix

	exnorm1	&prefix, &num&bs4&suffix
	exnorm1	&prefix, &num&ebs4&suffix
	exnorm1	&prefix, &num&cbs4&suffix
	exnorm1	&prefix, &num&ecbs4&suffix
	exnorm1	&prefix, &num&zpbs4&suffix
	exnorm1	&prefix, &num&zpebs4&suffix
	exnorm1	&prefix, &num&zpcbs4&suffix
	exnorm1	&prefix, &num&zpecbs4&suffix
ENDIF
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

; Jmptable addresses for gwinfo1 to return

IFNDEF X86_64
EXTRN	jmptable:DWORD
EXTRN	jmptablep:DWORD
ENDIF
EXTRN	xjmptable:DWORD
EXTRN	xjmptablep:DWORD

;
; More global variables needed in FFT routines
;

IFNDEF X86_64
_GWDATA SEGMENT PAGE PUBLIC 'DATA'
ELSE
_GWDATA SEGMENT PAGE
ENDIF

IFNDEF X86_64
prctab1	DD	OFFSET gwadd1, OFFSET gwaddq1, OFFSET gwsub1, OFFSET gwsubq1
	DD	OFFSET gwaddsub1, OFFSET gwaddsubq1, OFFSET gwcopyzero1
	DD	OFFSET gwaddq1, OFFSET gwsubq1, OFFSET gwaddsubq1
	DD	0, OFFSET gwmuls1 ;; Buggy, won't propogate carry over multiple words. pad on zeropads
	DD	OFFSET r1z, OFFSET r1ze, OFFSET i1z, OFFSET i1ze
	DD	OFFSET r1, OFFSET r1e, OFFSET r1c, OFFSET r1ec
	DD	OFFSET r1zp, OFFSET r1zpe, OFFSET r1zpc, OFFSET r1zpec
	DD	OFFSET i1, OFFSET i1e, OFFSET i1c, OFFSET i1ec
	DD	OFFSET i1zp, OFFSET i1zpe, OFFSET i1zpc, OFFSET i1zpec
prctab2	DD	OFFSET gwadd2, OFFSET gwaddq2, OFFSET gwsub2, OFFSET gwsubq2
	DD	OFFSET gwaddsub2, OFFSET gwaddsubq2, OFFSET gwcopyzero2
	DD	OFFSET gwaddq2, OFFSET gwsubq2, OFFSET gwaddsubq2
	DD	0, OFFSET gwmuls2 ;; Buggy, won't propogate carry over multiple words. pad on zeropads
	DD	OFFSET r2z, OFFSET r2ze, OFFSET i2z, OFFSET i2ze
	DD	OFFSET r2, OFFSET r2e, OFFSET r2c, OFFSET r2ec
	DD	OFFSET r2zp, OFFSET r2zpe, OFFSET r2zpc, OFFSET r2zpec
	DD	OFFSET i2, OFFSET i2e, OFFSET i2c, OFFSET i2ec
	DD	OFFSET i2zp, OFFSET i2zpe, OFFSET i2zpc, OFFSET i2zpec
ENDIF

xprctab1 DP	OFFSET gwxadd1, OFFSET gwxaddq1, OFFSET gwxsub1
	DP	OFFSET gwxsubq1, OFFSET gwxaddsub1, OFFSET gwxaddsubq1
	DP	OFFSET gwxcopyzero1
	DP	OFFSET gwxaddq1, OFFSET gwxsubq1, OFFSET gwxaddsubq1
	DP	OFFSET gwxadds1, OFFSET gwxmuls1
	DP	OFFSET xr1z, OFFSET xr1ze, OFFSET xi1z, OFFSET xi1ze
	DP	OFFSET xr1, OFFSET xr1e, OFFSET xr1c, OFFSET xr1ec
	DP	OFFSET xr1zp, OFFSET xr1zpe, OFFSET xr1zpc, OFFSET xr1zpec
	DP	OFFSET xi1, OFFSET xi1e, OFFSET xi1c, OFFSET xi1ec
	DP	OFFSET xi1zp, OFFSET xi1zpe, OFFSET xi1zpc, OFFSET xi1zpec
	DP	OFFSET xr1b, OFFSET xr1eb, OFFSET xr1cb, OFFSET xr1ecb
	DP	OFFSET xr1zpb, OFFSET xr1zpeb, OFFSET xr1zpcb, OFFSET xr1zpecb
	DP	OFFSET xi1b, OFFSET xi1eb, OFFSET xi1cb, OFFSET xi1ecb
	DP	OFFSET xi1zpb, OFFSET xi1zpeb, OFFSET xi1zpcb, OFFSET xi1zpecb
	DP	OFFSET xr1s4, OFFSET xr1es4, OFFSET xr1cs4, OFFSET xr1ecs4
	DP	OFFSET xr1zps4, OFFSET xr1zpes4, OFFSET xr1zpcs4, OFFSET xr1zpecs4
	DP	OFFSET xi1s4, OFFSET xi1es4, OFFSET xi1cs4, OFFSET xi1ecs4
	DP	OFFSET xi1zps4, OFFSET xi1zpes4, OFFSET xi1zpcs4, OFFSET xi1zpecs4
	DP	OFFSET xr1bs4, OFFSET xr1ebs4, OFFSET xr1cbs4, OFFSET xr1ecbs4
	DP	OFFSET xr1zpbs4, OFFSET xr1zpebs4, OFFSET xr1zpcbs4, OFFSET xr1zpecbs4
	DP	OFFSET xi1bs4, OFFSET xi1ebs4, OFFSET xi1cbs4, OFFSET xi1ecbs4
	DP	OFFSET xi1zpbs4, OFFSET xi1zpebs4, OFFSET xi1zpcbs4, OFFSET xi1zpecbs4
xprctab2 DP	OFFSET gwxadd2, OFFSET gwxaddq2, OFFSET gwxsub2
	DP	OFFSET gwxsubq2, OFFSET gwxaddsub2, OFFSET gwxaddsubq2
	DP	OFFSET gwxcopyzero2
	DP	OFFSET gwxaddq2, OFFSET gwxsubq2, OFFSET gwxaddsubq2
	DP	OFFSET gwxadds2, OFFSET gwxmuls2
	DP	OFFSET xr2z, OFFSET xr2ze, OFFSET xi2z, OFFSET xi2ze
	DP	OFFSET xr2, OFFSET xr2e, OFFSET xr2c, OFFSET xr2ec
	DP	OFFSET xr2zp, OFFSET xr2zpe, OFFSET xr2zpc, OFFSET xr2zpec
	DP	OFFSET xi2, OFFSET xi2e, OFFSET xi2c, OFFSET xi2ec
	DP	OFFSET xi2zp, OFFSET xi2zpe, OFFSET xi2zpc, OFFSET xi2zpec
	DP	OFFSET xr2b, OFFSET xr2eb, OFFSET xr2cb, OFFSET xr2ecb
	DP	OFFSET xr2zpb, OFFSET xr2zpeb, OFFSET xr2zpcb, OFFSET xr2zpecb
	DP	OFFSET xi2b, OFFSET xi2eb, OFFSET xi2cb, OFFSET xi2ecb
	DP	OFFSET xi2zpb, OFFSET xi2zpeb, OFFSET xi2zpcb, OFFSET xi2zpecb
	DP	OFFSET xr2s4, OFFSET xr2es4, OFFSET xr2cs4, OFFSET xr2ecs4
	DP	OFFSET xr2zps4, OFFSET xr2zpes4, OFFSET xr2zpcs4, OFFSET xr2zpecs4
	DP	OFFSET xi2s4, OFFSET xi2es4, OFFSET xi2cs4, OFFSET xi2ecs4
	DP	OFFSET xi2zps4, OFFSET xi2zpes4, OFFSET xi2zpcs4, OFFSET xi2zpecs4
	DP	OFFSET xr2bs4, OFFSET xr2ebs4, OFFSET xr2cbs4, OFFSET xr2ecbs4
	DP	OFFSET xr2zpbs4, OFFSET xr2zpebs4, OFFSET xr2zpcbs4, OFFSET xr2zpecbs4
	DP	OFFSET xi2bs4, OFFSET xi2ebs4, OFFSET xi2cbs4, OFFSET xi2ecbs4
	DP	OFFSET xi2zpbs4, OFFSET xi2zpebs4, OFFSET xi2zpcbs4, OFFSET xi2zpecbs4
xprctab2a DP	OFFSET gwxadd2, OFFSET gwxaddq2, OFFSET gwxsub2
	DP	OFFSET gwxsubq2, OFFSET gwxaddsub2, OFFSET gwxaddsubq2
	DP	OFFSET gwxcopyzero2
	DP	OFFSET gwxaddq2, OFFSET gwxsubq2, OFFSET gwxaddsubq2
	DP	OFFSET gwxadds2, OFFSET gwxmuls2
	DP	OFFSET xr2zAMD, OFFSET xr2zeAMD, OFFSET xi2zAMD, OFFSET xi2zeAMD
	DP	OFFSET xr2AMD, OFFSET xr2eAMD, OFFSET xr2cAMD, OFFSET xr2ecAMD
	DP	OFFSET xr2zpAMD, OFFSET xr2zpeAMD, OFFSET xr2zpcAMD, OFFSET xr2zpecAMD
	DP	OFFSET xi2AMD, OFFSET xi2eAMD, OFFSET xi2cAMD, OFFSET xi2ecAMD
	DP	OFFSET xi2zpAMD, OFFSET xi2zpeAMD, OFFSET xi2zpcAMD, OFFSET xi2zpecAMD
	DP	OFFSET xr2bAMD, OFFSET xr2ebAMD, OFFSET xr2cbAMD, OFFSET xr2ecbAMD
	DP	OFFSET xr2zpbAMD, OFFSET xr2zpebAMD, OFFSET xr2zpcbAMD, OFFSET xr2zpecbAMD
	DP	OFFSET xi2bAMD, OFFSET xi2ebAMD, OFFSET xi2cbAMD, OFFSET xi2ecbAMD
	DP	OFFSET xi2zpbAMD, OFFSET xi2zpebAMD, OFFSET xi2zpcbAMD, OFFSET xi2zpecbAMD
	DP	OFFSET xr2s4AMD, OFFSET xr2es4AMD, OFFSET xr2cs4AMD, OFFSET xr2ecs4AMD
	DP	OFFSET xr2zps4AMD, OFFSET xr2zpes4AMD, OFFSET xr2zpcs4AMD, OFFSET xr2zpecs4AMD
	DP	OFFSET xi2s4AMD, OFFSET xi2es4AMD, OFFSET xi2cs4AMD, OFFSET xi2ecs4AMD
	DP	OFFSET xi2zps4AMD, OFFSET xi2zpes4AMD, OFFSET xi2zpcs4AMD, OFFSET xi2zpecs4AMD
	DP	OFFSET xr2bs4AMD, OFFSET xr2ebs4AMD, OFFSET xr2cbs4AMD, OFFSET xr2ecbs4AMD
	DP	OFFSET xr2zpbs4AMD, OFFSET xr2zpebs4AMD, OFFSET xr2zpcbs4AMD, OFFSET xr2zpecbs4AMD
	DP	OFFSET xi2bs4AMD, OFFSET xi2ebs4AMD, OFFSET xi2cbs4AMD, OFFSET xi2ecbs4AMD
	DP	OFFSET xi2zpbs4AMD, OFFSET xi2zpebs4AMD, OFFSET xi2zpcbs4AMD, OFFSET xi2zpecbs4AMD

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
