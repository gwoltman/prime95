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

_TEXT SEGMENT

INCLUDE extrn.mac
INCLUDE	unravel.mac
INCLUDE xmult.mac
INCLUDE	xlucas.mac
INCLUDE xpass2.mac

PREFETCHING = 1

;; Routines to do the last 8 levels in a two-pass FFT

PUBLIC	xpass2_8_levels
PUBLIC	xpass2_8_levels_p

;; Distance between two pass 2 data blocks.  Pass 2 does 8 FFT levels
;; 2 sets of data (2 * 2^8 complex values = 2^10 doubles = 8KB).

blkdst	EQU	(8192+128)


;; Do the last 8 levels of a two pass FFT

xpass2_8_levels PROC
	xpass2_8_levels_real blkdst
xpass2_8_levels_p:
	mov	ecx, count1		; Number of complex iterations
	mov	edx, pass2_premults	; Address of the group multipliers
p2lp:	xpass2_8_levels_complex blkdst
	sub	ecx, 1
	jnz	p2lp
	mov	esi, _DESTARG		; Restore source pointer
	ret
xpass2_8_levels ENDP

_TEXT	ENDS
END
