; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 11 levels of two pass FFTs that use the
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

IFDEF X86_64
INCLUDE	xlucas64.mac
ENDIF

;; Routines to do the last 11 levels in a two-pass FFT

PUBLICP	xpass2_11_levels
PUBLICP	xpass2_11_levels_p

_TEXT SEGMENT

PREFETCHING = 1

;; Do the last 11 levels of a two pass FFT

PROCP	xpass2_11_levels
	start_timer 2
	xpass2_11_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_11_levels_p
	start_timer 3
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p2lp:	xpass2_11_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p2lp
	mov	rsi, _DESTARG		; Restore source pointer
	end_timer 3
	ret
ENDPP	xpass2_11_levels

;;
;; Routines to do the normalization after a multiply
;;

PROCP _xnorm3

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	noadd, ilp0, ilp1, ilexit
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
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	rax, rax		;; Clear big/lit flags
	sub	rcx, rcx
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	mov	rbx, norm_ptr2		;; Load column multipliers ptr
ilp1:	xprefetchw [rsi+64]
	xnorm_2d ttp, zero, echk, const	;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbx, [rbx+32]		;; Next column multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
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
ilexit:	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, rdi		;; Save big/little flags array ptr
ttp	mov	norm_ptr2, rbx		;; Save column multipliers ptr
	pop	rsi
	pop	rdx
	sub	rbx, rbx
	cmp	edx, 256		;; Check for last iteration
	je	top_carry_adjust_3	;; Top carry may require adjusting
	ret
	ENDM

zpnorm	MACRO	lab, ttp, echk, const
	LOCAL	ilp0, ilp1, ilexit
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
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	rax, rax		;; Clear big/lit flags
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	mov	rbx, norm_ptr2		;; Load column multipliers ptr
ilp1:	xprefetchw [rsi+64]
	xnorm_2d_zpad ttp, echk, const	;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbx, [rbx+32]		;; Next column multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
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
ilexit:	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
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

	inorm	xr3, noexec, noexec, noexec, noexec
	inorm	xr3e, noexec, noexec, exec, noexec
	inorm	xr3c, noexec, noexec, noexec, exec
	inorm	xr3ec, noexec, noexec, exec, exec
	inorm	xr3z, noexec, exec, noexec, noexec
	inorm	xr3ze, noexec, exec, exec, noexec
	inorm	xi3, exec, noexec, noexec, noexec
	inorm	xi3e, exec, noexec, exec, noexec
	inorm	xi3c, exec, noexec, noexec, exec
	inorm	xi3ec, exec, noexec, exec, exec
	inorm	xi3z, exec, exec, noexec, noexec
	inorm	xi3ze, exec, exec, exec, noexec

	zpnorm	xr3zp, noexec, noexec, noexec
	zpnorm	xr3zpe, noexec, exec, noexec
	zpnorm	xr3zpc, noexec, noexec, exec
	zpnorm	xr3zpec, noexec, exec, exec
	zpnorm	xi3zp, exec, noexec, noexec
	zpnorm	xi3zpe, exec, exec, noexec
	zpnorm	xi3zpc, exec, noexec, exec
	zpnorm	xi3zpec, exec, exec, exec

; Special code to handle adjusting the carry out of the topmost FFT word

top_carry_adjust_3:
	xnorm_top_carry
	ret

ENDPP _xnorm3

_TEXT	ENDS
END
