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
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p11lp:	xpass2_11_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p11lp
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_11_levels



;; Routines to do the last 12 levels in a two-pass FFT

PUBLICP	xpass2_12_levels
PUBLICP	xpass2_12_levels_p

;; Do the last 12 levels of a two pass FFT

PROCP	xpass2_12_levels
	start_timer 2
	xpass2_12_levels_real 0
	add	rsi, pass1blkdst
	end_timer 2
LABELP	xpass2_12_levels_p
	mov	ecx, count1		; Number of complex iterations
	mov	rdx, pass2_premults	; Address of the group multipliers
p12lp:	xpass2_12_levels_complex 0
	add	rsi, pass1blkdst
	sub	ecx, 1
	jnz	p12lp
	mov	rsi, _DESTARG		; Restore source pointer
	ret
ENDPP	xpass2_12_levels


_TEXT	ENDS
END
