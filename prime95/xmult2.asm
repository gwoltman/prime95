; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes between than 5K and 128K doubles are supported.
; This code does two passes, 8 levels on the second pass.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
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
INCLUDE xfft2.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac
INCLUDE xnormal.mac

EXTRN	xpass2_8_levels:PROC
EXTRN	xpass2_8_levels_p:PROC

;; Distance between two pass 2 data blocks.  Pass 2 does 8 FFT levels
;; 2 sets of data (2 * 2^8 complex values = 2^10 doubles = 8KB).

blkdst	EQU	(8192+128)

;; All the FFT routines for each FFT length

_xmm_gw_ffts2 PROC
	EXPANDING = 2
;	xfft	5K
;	xfft	6K
;	xfft	7K
;	xfft	8K
	xfft	10K
	xfft	12K
	xfft	12Kp
	xfft	14K
	xfft	16K
	xfft	16Kp
	xfft	20K
	xfft	24K
	xfft	24Kp
	xfft	28K
	xfft	32K
	xfft	32Kp
;	xfft	40K
;	xfft	48K
;	xfft	56K
;	xfft	64K
;	xfft	80K
;	xfft	96K
;	xfft	112K
;	xfft	128K

; Common code to finish the FFT by restoring state and returning.

gw_finish_fft_2:
	mov	DWORD PTR [esi-28], 3	; Set has-been-FFTed flags
	xfft_1_ret

; Common code to finish up multiplies

; Split the accumulated carries into two carries - a high carry and a
; low carry.  Handle both the with and without two-to-phi array cases.
; Add these carries back into the FFT data.

	PUBLIC	gw_carries_2
gw_carries_2:
	cmp	_ZERO_PADDED_FFT, 0	; Special case the zero padded FFT case
	jne	gw_carries_2_zpad
	mov	esi, carries		; Addr of the carries
	mov	ebx, addcount1		; Compute end of carries addr
	shl	ebx, 6
	add	ebx, esi
	xnorm012_2d_part1
	mov	ebp, _DESTARG		; Addr of the FFT data
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	edx, norm_grp_mults	; Addr of the group multipliers
	sub	ecx, ecx		; Clear big/little flag
	mov	eax, count3		; Load 3 section counts
	mov	loopcount1, eax		; Save for later
ilp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	ebx, 03FFh
	jz	spldn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 10			; Move counts list along
	mov	loopcount1, eax
	shl	ebx, 6			; Compute addr of the last carries row
	add	ebx, esi
	xnorm012_2d_part2
	sub	eax, eax		; Clear big/little flag
ilp1:	mov	ebx, norm_col_mults	; Addr of the column multipliers
	xnorm012_2d			; Split carries for one cache line
	mov	ebx, cache_line_multiplier ; Cache lines in each pass1 loop
	lea	esi, [esi+64]		; Next carries pointer
	lea	ebp, [ebp+blkdst]	; Next FFT data pointer
	cmp	_RATIONAL_FFT, 0	; Don't bump these two pointers
	jne	short iskip		; for rational FFTs
	lea	edx, [edx+128]		; Next group multiplier
	lea	edi, [edi+ebx*4]	; Next big/little flags pointer
iskip:	sub	loopcount2, 1		; Test loop counter
	jnz	ilp1			; Next carry row in section
	jmp	ilp0			; Next section
spldn:	mov	zero_fft, 0		; Clear zero-high-words-fft flag
	ret

gw_carries_2_zpad:
	mov	esi, carries		; Addr of the carries
	mov	ebx, addcount1		; Compute end of carries addr
	shl	ebx, 6
	add	ebx, esi
	mov	esi, _DESTARG		; Addr of the FFT data
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	ebp, norm_col_mults	; Addr of the group multipliers
	cmp	const_fft, 0		; Call correct part1 macro
	je	c2a			; Jump if not const
	xnorm012_2d_zpad_part1 exec
	jmp	c2b
c2a:	xnorm012_2d_zpad_part1 noexec
c2b:	mov	esi, carries		; Addr of the carries
	mov	ebp, _DESTARG		; Addr of the FFT data
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	edx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load 3 section counts
	mov	loopcount1, eax		; Save for later
zlp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	ebx, 03FFh
	jz	zpldn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 10			; Move counts list along
	mov	loopcount1, eax
	shl	ebx, 6			; Compute addr of the last carries row
	add	ebx, esi
	xnorm012_2d_zpad_part2
	sub	eax, eax		; Clear big/little flag
zlp1:	mov	ebx, norm_col_mults	; Addr of the column multipliers
	cmp	const_fft, 0		; Call correct zpad macro
	je	c2c			; Jump if not const
	xnorm012_2d_zpad exec		; Split carries for one cache line
	jmp	c2d
c2c:	xnorm012_2d_zpad noexec		; Split carries for one cache line
c2d:	mov	ebx, cache_line_multiplier ; Cache lines in each pass1 loop
	lea	esi, [esi+64]		; Next carries pointer
	lea	ebp, [ebp+blkdst]	; Next FFT data pointer
	cmp	_RATIONAL_FFT, 0	; Don't bump these two pointers
	jne	short zskip		; for rational FFTs
	lea	edx, [edx+128]		; Next group multiplier
	lea	edi, [edi+ebx*4]	; Next big/little flags pointer
zskip:	sub	loopcount2, 1		; Test loop counter
	jnz	zlp1			; Next carry row in section
	jmp	zlp0			; Next section
zpldn:	mov	const_fft, 0		; Clear mul-by-const-fft flag
	ret

; Finish the multiply

gw_finish_mult_2:

; Set FFT-started flag

cmnend:	mov	esi, _DESTARG		; Addr of FFT data
	mov	eax, _POSTFFT		; Set FFT started flag
	mov	DWORD PTR [esi-28], eax

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

	_movsd	xmm7, XMM_SUMOUT	; Add together the two partial sumouts
	addsd	xmm7, XMM_SUMOUT+8
	mulsd	xmm7, ttmp_ff_inv
	_movsd	[esi-24], xmm7		; Save sum of FFT outputs
	_movsd	xmm6, XMM_MAXERR	; Compute new maximum error
	maxsd	xmm6, XMM_MAXERR+8
	_movsd	_MAXERR, xmm6

; Return

exit:	pop	esi
	pop	edi
	pop	ebx
	pop	ebp
	ret

_xmm_gw_ffts2 ENDP

_TEXT	ENDS
END
