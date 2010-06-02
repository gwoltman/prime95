; Copyright 2001-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These routines implement the normalization part of a r4dwpn (radix-4 delayed with partial normalization) FFT.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE xnormal.mac

_TEXT SEGMENT

;;
;; Routines to do the normalization after a multiply
;;

;; When doing zero-padded FFTs, the multiplied 7 words around the halfway point
;; must be subtracted from the bottom of the FFT.  This must be done before
;; normalization multiplies the FFT data by k.  This macro does that.

xsub_7_words MACRO
	LOCAL	nozpad, zlp
	cmp	THIS_BLOCK, 8		;; Have we subtracted all 7 words?
	jge	short nozpad		;; Yes, skip this code
	mov	eax, cache_line_multiplier ;; Load loop counter
	mov	edx, THIS_BLOCK
	mov	rcx, rsi		;; Copy source ptr (we preserve rsi)
	mov	rdi, zpad_addr		;; Addr of first zpad element
zlp:	movsd	xmm0, Q [rcx]		;; Load FFT word
	movsd	xmm1, Q [rdi][rdx*8]	;; Load ZPAD data
	subsd	xmm0, xmm1
	movsd	Q [rcx], xmm0		;; Store FFT word
	bump	rcx, 64			;; Bump pointers
	inc	rdx
	dec	eax			;; Iterate 2*clm (up to 8) times
	jnz	short zlp		;; Loop if necessary
nozpad:
	ENDM

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

saved_rsi	EQU	PPTR [rsp+first_local]
loopcount1	EQU	DPTR [rsp+first_local+SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+8]

inorm	MACRO	lab, ttp, zero, echk, const, base2, sse4
	LOCAL	noadd, setlp, ilp0, ilp1, ilexit, done
	PROCFLP	lab
	int_prolog SZPTR+12,0,0
zero	mov	zero_fft, 1		;; Set flag saying zero upper half
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
no zero	mov	edx, ADDIN_ROW		;; Is this the time to do our addin?
no zero	cmp	edx, THIS_BLOCK
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	movsd	xmm0, Q [rsi][rdi]	;; Get the value
no zero	addsd	xmm0, ADDIN_VALUE	;; Add in the requested value
no zero	movsd	Q [rsi][rdi], xmm0	;; Save the new value
no zero	subsd	xmm7, ADDIN_VALUE	;; Do not include addin in sumout
noadd:	mov	saved_rsi, rsi		;; Save for xtop_carry_adjust

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	rax, rax		;; Clear big/lit flags
	sub	rcx, rcx
ttp	mov	al, [rdi+0]		;; Load big vs. little flags
ttp	mov	cl, [rdi+1]		;; Load big vs. little flags
IFDEF X86_64
	sub	r8, r8
	sub	r9, r9
ttp	mov	r8b, [rdi+2]		;; Load big vs. little flags
ttp	mov	r9b, [rdi+3]		;; Load big vs. little flags
ENDIF
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	xprefetcht1 [rdx+128]		;; Prefetch group multiplier
ilp1:	xprefetchw [rsi+64]
	xnorm_wpn ttp, zero, echk, const, base2, sse4 ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rdi, 4			;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	add	rsi, normblkdst		;; Skip gap in blkdst or clmblkdst
	lea	rbp, [rbp+64]		;; Next set of carries
ttp	lea	rdx, [rdx+128]		;; Next set of 8 group multipliers
	sub	loopcount2, 1		;; Test loop counter
	jz	ilexit			;; Jump when loop complete
	add	loopcount3, 80000000h/4 ;; 8 iterations
	jnc	ilp0
	add	rsi, normblkdst8	;; Add 128 every 8 clmblkdsts
	jmp	ilp0			;; Iterate
ilexit:	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, rdi		;; Save big/little flags array ptr

	; Handle adjusting the carry out of the topmost FFT word

	mov	eax, THIS_BLOCK		;; Check for processing last block
	cmp	eax, LAST_PASS1_BLOCK
	jne	done			;; Jump if not last block
	mov	rsi, saved_rsi		;; Restore FFT data ptr
	xnorm_top_carry			;; Adjust carry if k > 1

done:	int_epilog SZPTR+12,0,0
	ENDPP	lab
	ENDM

loopcount1z	EQU	DPTR [rsp+first_local]
loopcount2z	EQU	DPTR [rsp+first_local+4]
loopcount3z	EQU	DPTR [rsp+first_local+8]

zpnorm	MACRO	lab, ttp, echk, const, base2, sse4, khi, c1, cm1
	LOCAL	setlp, ilp0, ilp1, ilexit
	PROCFLP	lab
	int_prolog 12,0,0
const	mov	const_fft, 1		;; Set flag saying mul-by-const
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
	xsub_7_words

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2z, eax	;; Save loop counter
	mov	loopcount3z, 0		;; Clear outermost loop counter
	sub	rax, rax		;; Clear big/lit flags
ttp	mov	al, [rdi+0]		;; Load big vs. little flags
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1z, ebx	;; Save loop counter
	xprefetcht1 [rdx+128]		;; Prefetch group multiplier
ilp1:	xprefetchw [rsi+64]
	xnorm_wpn_zpad ttp, echk, const, base2, sse4, khi, c1, cm1 ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rdi, 4			;; Next big/little flags
	sub	loopcount1z, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	add	rsi, normblkdst		;; Skip gap in blkdst or clmblkdst
	bump	rbp, 64			;; Next set of carries
ttp	bump	rdx, 128		;; Next set of 8 group multipliers
	sub	loopcount2z, 1		;; Test loop counter
	jz	ilexit			;; Jump when loop complete
	add	loopcount3z, 80000000h/4 ;; 8 iterations
	jnc	ilp0
	add	rsi, normblkdst8	;; Add 128 every 8 clmblkdsts
	jmp	ilp0			;; Iterate
ilexit:	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, rdi		;; Save big/little flags array ptr
	int_epilog 12,0,0
	ENDPP	lab
	ENDM

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

PREFETCHING = 64

	inorm	xr3, noexec, noexec, noexec, noexec, exec, noexec
	inorm	xr3e, noexec, noexec, exec, noexec, exec, noexec
	inorm	xr3c, noexec, noexec, noexec, exec, exec, noexec
	inorm	xr3ec, noexec, noexec, exec, exec, exec, noexec
	inorm	xr3z, noexec, exec, noexec, noexec, exec, noexec
	inorm	xr3ze, noexec, exec, exec, noexec, exec, noexec
	inorm	xi3, exec, noexec, noexec, noexec, exec, noexec
	inorm	xi3e, exec, noexec, exec, noexec, exec, noexec
	inorm	xi3c, exec, noexec, noexec, exec, exec, noexec
	inorm	xi3ec, exec, noexec, exec, exec, exec, noexec
	inorm	xi3z, exec, exec, noexec, noexec, exec, noexec
	inorm	xi3ze, exec, exec, exec, noexec, exec, noexec
	zpnorm	xr3zp, noexec, noexec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xr3zpc1, noexec, noexec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xr3zpcm1, noexec, noexec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xr3zpe, noexec, exec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xr3zpec1, noexec, exec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xr3zpecm1, noexec, exec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xr3zpc, noexec, noexec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xr3zpec, noexec, exec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xi3zp, exec, noexec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xi3zpc1, exec, noexec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xi3zpcm1, exec, noexec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xi3zpe, exec, exec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xi3zpec1, exec, exec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xi3zpecm1, exec, exec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xi3zpc, exec, noexec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xi3zpec, exec, exec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xr3zpk, noexec, noexec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xr3zpkc1, noexec, noexec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xr3zpkcm1, noexec, noexec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xr3zpek, noexec, exec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xr3zpekc1, noexec, exec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xr3zpekcm1, noexec, exec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xr3zpck, noexec, noexec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	xr3zpeck, noexec, exec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi3zpk, exec, noexec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi3zpkc1, exec, noexec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xi3zpkcm1, exec, noexec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xi3zpek, exec, exec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi3zpekc1, exec, exec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xi3zpekcm1, exec, exec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xi3zpck, exec, noexec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi3zpeck, exec, exec, exec, exec, noexec, noexec, noexec, noexec

	inorm	xr3b, noexec, noexec, noexec, noexec, noexec, noexec
	inorm	xr3eb, noexec, noexec, exec, noexec, noexec, noexec
	inorm	xr3cb, noexec, noexec, noexec, exec, noexec, noexec
	inorm	xr3ecb, noexec, noexec, exec, exec, noexec, noexec
	inorm	xi3b, exec, noexec, noexec, noexec, noexec, noexec
	inorm	xi3eb, exec, noexec, exec, noexec, noexec, noexec
	inorm	xi3cb, exec, noexec, noexec, exec, noexec, noexec
	inorm	xi3ecb, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xr3zpb, noexec, noexec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr3zpbc1, noexec, noexec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xr3zpbcm1, noexec, noexec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xr3zpeb, noexec, exec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr3zpebc1, noexec, exec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xr3zpebcm1, noexec, exec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xr3zpcb, noexec, noexec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr3zpecb, noexec, exec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi3zpb, exec, noexec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi3zpbc1, exec, noexec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xi3zpbcm1, exec, noexec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xi3zpeb, exec, exec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi3zpebc1, exec, exec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xi3zpebcm1, exec, exec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xi3zpcb, exec, noexec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi3zpecb, exec, exec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr3zpbk, noexec, noexec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xr3zpbkc1, noexec, noexec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xr3zpbkcm1, noexec, noexec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xr3zpebk, noexec, exec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xr3zpebkc1, noexec, exec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xr3zpebkcm1, noexec, exec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xr3zpcbk, noexec, noexec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xr3zpecbk, noexec, exec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi3zpbk, exec, noexec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi3zpbkc1, exec, noexec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xi3zpbkcm1, exec, noexec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xi3zpebk, exec, exec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi3zpebkc1, exec, exec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xi3zpebkcm1, exec, exec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xi3zpcbk, exec, noexec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi3zpecbk, exec, exec, exec, noexec, noexec, noexec, noexec, noexec

	inorm	xr3s4, noexec, noexec, noexec, noexec, exec, exec
	inorm	xr3es4, noexec, noexec, exec, noexec, exec, exec
	inorm	xr3cs4, noexec, noexec, noexec, exec, exec, exec
	inorm	xr3ecs4, noexec, noexec, exec, exec, exec, exec
	inorm	xi3s4, exec, noexec, noexec, noexec, exec, exec
	inorm	xi3es4, exec, noexec, exec, noexec, exec, exec
	inorm	xi3cs4, exec, noexec, noexec, exec, exec, exec
	inorm	xi3ecs4, exec, noexec, exec, exec, exec, exec
	zpnorm	xr3zps4, noexec, noexec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xr3zps4c1, noexec, noexec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xr3zps4cm1, noexec, noexec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xr3zpes4, noexec, exec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xr3zpes4c1, noexec, exec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xr3zpes4cm1, noexec, exec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xr3zpcs4, noexec, noexec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xr3zpecs4, noexec, exec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xi3zps4, exec, noexec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xi3zps4c1, exec, noexec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xi3zps4cm1, exec, noexec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xi3zpes4, exec, exec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xi3zpes4c1, exec, exec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xi3zpes4cm1, exec, exec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xi3zpcs4, exec, noexec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xi3zpecs4, exec, exec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xr3zps4k, noexec, noexec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xr3zps4kc1, noexec, noexec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xr3zps4kcm1, noexec, noexec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xr3zpes4k, noexec, exec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xr3zpes4kc1, noexec, exec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xr3zpes4kcm1, noexec, exec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xr3zpcs4k, noexec, noexec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	xr3zpecs4k, noexec, exec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	xi3zps4k, exec, noexec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xi3zps4kc1, exec, noexec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xi3zps4kcm1, exec, noexec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xi3zpes4k, exec, exec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xi3zpes4kc1, exec, exec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xi3zpes4kcm1, exec, exec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xi3zpcs4k, exec, noexec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	xi3zpecs4k, exec, exec, exec, exec, exec, noexec, noexec, noexec

	inorm	xr3bs4, noexec, noexec, noexec, noexec, noexec, exec
	inorm	xr3ebs4, noexec, noexec, exec, noexec, noexec, exec
	inorm	xr3cbs4, noexec, noexec, noexec, exec, noexec, exec
	inorm	xr3ecbs4, noexec, noexec, exec, exec, noexec, exec
	inorm	xi3bs4, exec, noexec, noexec, noexec, noexec, exec
	inorm	xi3ebs4, exec, noexec, exec, noexec, noexec, exec
	inorm	xi3cbs4, exec, noexec, noexec, exec, noexec, exec
	inorm	xi3ecbs4, exec, noexec, exec, exec, noexec, exec
	zpnorm	xr3zpbs4, noexec, noexec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xr3zpbs4c1, noexec, noexec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xr3zpbs4cm1, noexec, noexec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xr3zpebs4, noexec, exec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xr3zpebs4c1, noexec, exec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xr3zpebs4cm1, noexec, exec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xr3zpcbs4, noexec, noexec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xr3zpecbs4, noexec, exec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xi3zpbs4, exec, noexec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xi3zpbs4c1, exec, noexec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xi3zpbs4cm1, exec, noexec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xi3zpebs4, exec, exec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xi3zpebs4c1, exec, exec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xi3zpebs4cm1, exec, exec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xi3zpcbs4, exec, noexec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xi3zpecbs4, exec, exec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xr3zpbs4k, noexec, noexec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xr3zpbs4kc1, noexec, noexec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xr3zpbs4kcm1, noexec, noexec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xr3zpebs4k, noexec, exec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xr3zpebs4kc1, noexec, exec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xr3zpebs4kcm1, noexec, exec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xr3zpcbs4k, noexec, noexec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	xr3zpecbs4k, noexec, exec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi3zpbs4k, exec, noexec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi3zpbs4kc1, exec, noexec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xi3zpbs4kcm1, exec, noexec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xi3zpebs4k, exec, exec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi3zpebs4kc1, exec, exec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xi3zpebs4kcm1, exec, exec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xi3zpcbs4k, exec, noexec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi3zpecbs4k, exec, exec, exec, noexec, exec, noexec, noexec, noexec

_TEXT	ENDS
END
