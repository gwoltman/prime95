; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 8 levels of two pass FFTs that use the
; SSE2 instructions.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xmult.mac
INCLUDE	xlucas.mac
INCLUDE xpass2.mac
INCLUDE xnormal.mac

;; Routines to do the last 8 levels in a two-pass FFT

PUBLICP	xpass2_8_levels
PUBLICP	xpass2_8_levels_p

_TEXT SEGMENT

PREFETCHING = 1

;; Distance between two pass 2 data blocks.  Pass 2 does 8 FFT levels
;; 2 sets of data (2 * 2^8 complex values = 2^10 doubles = 8KB).

blkdst	EQU	(8192+128)


;; Do the last 8 levels of a two pass FFT

PROCP	xpass2_8_levels
	xpass2_8_levels_real blkdst
LABELP	xpass2_8_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p2lp:	xpass2_8_levels_complex blkdst
	sub	ecx, 1
	jnz	p2lp
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_8_levels

;;
;; Routines to do the normalization after a multiply
;;

PROCP _xnorm2

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	noadd, ilp0, ilp1
	PUBLICP	lab
	LABELP	lab
zero	mov	zero_fft, 1		;; Set flag saying zero upper half
	movapd	xmm7, XMM_SUMOUT	;; Load SUMOUT
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
no zero	cmp	edx, _ADDIN_ROW		;; Is this the time to do our addin?
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, _ADDIN_OFFSET	;; Get address to add value into
no zero	_movsd	xmm0, [rsi][rdi]	;; Get the value
no zero	addsd	xmm0, _ADDIN_VALUE	;; Add in the requested value
no zero	_movsd	[rsi][rdi], xmm0	;; Save the new value
no zero	subsd	xmm7, _ADDIN_VALUE	;; Do not include addin in sumout
noadd:	push	rdx
	push	rsi
	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	sub	rax, rax		;; Clear big/lit flags
	sub	rcx, rcx
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	mov	rbx, norm_ptr2		;; Load column multipliers ptr
	push	rsi
ilp1:	xprefetchw [rsi+blkdst]
	xnorm_2d ttp, zero, echk, const	;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbx, [rbx+32]		;; Next column multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	pop	rsi
	lea	rsi, [rsi+blkdst]	;; Next source pointer
	lea	rbp, [rbp+64]		;; Next set of carries
ttp	lea	rdx, [rdx+128]		;; Next set of 8 group multipliers
	sub	loopcount2, 1		;; Test loop counter
	jnz	ilp0
	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, rdi		;; Save big/little flags array ptr
ttp	mov	norm_ptr2, rbx		;; Save column multipliers ptr
	pop	rsi
	pop	rdx
	sub	rbx, rbx
	cmp	edx, 1			;; Check for last iteration
	je	top_carry_adjust_2	;; Top carry may require adjusting
	ret
	ENDM

zpnorm	MACRO	lab, ttp, echk, const
	LOCAL	ilp0, ilp1
	PUBLICP	lab
	LABELP	lab
const	mov	const_fft, 1		;; Set flag saying mul-by-const
	movapd	xmm7, XMM_SUMOUT	;; Load SUMOUT
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
	push	rdx
	push	rsi
	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	sub	rax, rax		;; Clear big/lit flags
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	mov	rbx, norm_ptr2		;; Load column multipliers ptr
	push	rsi
ilp1:	xprefetchw [rsi+blkdst]
	xnorm_2d_zpad ttp, echk, const	;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbx, [rbx+32]		;; Next column multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	pop	rsi
	lea	rsi, [rsi+blkdst]	;; Next source pointer
	lea	rbp, [rbp+64]		;; Next set of carries
ttp	lea	rdx, [rdx+128]		;; Next set of 8 group multipliers
	sub	loopcount2, 1		;; Test loop counter
	jnz	ilp0
	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, rdi		;; Save big/little flags array ptr
ttp	mov	norm_ptr2, rbx		;; Save column multipliers ptr
	pop	rsi
	pop	rdx
	sub	rbx, rbx
	ret
	ENDM

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

	inorm	xr2, noexec, noexec, noexec, noexec
	inorm	xr2e, noexec, noexec, exec, noexec
	inorm	xr2c, noexec, noexec, noexec, exec
	inorm	xr2ec, noexec, noexec, exec, exec
	inorm	xr2z, noexec, exec, noexec, noexec
	inorm	xr2ze, noexec, exec, exec, noexec
	inorm	xi2, exec, noexec, noexec, noexec
	inorm	xi2e, exec, noexec, exec, noexec
	inorm	xi2c, exec, noexec, noexec, exec
	inorm	xi2ec, exec, noexec, exec, exec
	inorm	xi2z, exec, exec, noexec, noexec
	inorm	xi2ze, exec, exec, exec, noexec

	zpnorm	xr2zp, noexec, noexec, noexec
	zpnorm	xr2zpe, noexec, exec, noexec
	zpnorm	xr2zpc, noexec, noexec, exec
	zpnorm	xr2zpec, noexec, exec, exec
	zpnorm	xi2zp, exec, noexec, noexec
	zpnorm	xi2zpe, exec, exec, noexec
	zpnorm	xi2zpc, exec, noexec, exec
	zpnorm	xi2zpec, exec, exec, exec

; Special code to handle adjusting the carry out of the topmost FFT word

top_carry_adjust_2:
	xnorm_top_carry
	ret

ENDPP _xnorm2

_TEXT	ENDS
END
