; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes of 40K and above are supported.
; This code does two passes, 11 levels on the second pass.
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

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xfft3.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac
INCLUDE xnormal.mac

EXTRN	xpass2_11_levels:PROC
EXTRN	xpass2_11_levels_p:PROC

_TEXT SEGMENT

;; All the FFT routines for each FFT length

_xmm_gw_ffts3 PROC
	EXPANDING = 2
	xfft	40K
	xfft	48K
	xfft	48Kp
	xfft	56K
	xfft	64K
	xfft	64Kp
	xfft	80K
	xfft	96K
	xfft	96Kp
	xfft	112K
	xfft	128K
	xfft	128Kp
	xfft	160K
	xfft	192K
	xfft	192Kp
	xfft	224K
	xfft	256K
	xfft	256Kp
	xfft	320K
	xfft	384K
	xfft	384Kp
	xfft	448K
	xfft	512K
	xfft	512Kp
	EXPANDING = 3
	xfftclm	640K, 4
	xfftclm	640K, 2
	xfftclm	768K, 4
	xfftclm	768K, 2
	xfftclm	768Kp, 4
	xfftclm	768Kp, 2
	xfftclm	896K, 4
	xfftclm	896K, 2
	xfftclm	896K, 1
	xfftclm	1024K, 4
	xfftclm	1024K, 2
	xfftclm	1024K, 1
	xfftclm	1024Kp, 4
	xfftclm	1024Kp, 2
	xfftclm	1024Kp, 1
	xfftclm	1280K, 4
	xfftclm	1280K, 2
	xfftclm	1280K, 1
	xfftclm	1536K, 4
	xfftclm	1536K, 2
	xfftclm	1536K, 1
	xfftclm	1536Kp, 4
	xfftclm	1536Kp, 2
	xfftclm	1536Kp, 1
	xfftclm	1792K, 4
	xfftclm	1792K, 2
	xfftclm	1792K, 1
	xfftclm	1792K, 0
	xfftclm	2048K, 4
	xfftclm	2048K, 2
	xfftclm	2048K, 1
	xfftclm	2048K, 0
	xfftclm	2048Kp, 4
	xfftclm	2048Kp, 2
	xfftclm	2048Kp, 1
	xfftclm	2048Kp, 0
	xfftclm	2560K, 2
	xfftclm	2560K, 1
	xfftclm	2560K, 0
	xfftclm	3072K, 2
	xfftclm	3072K, 1
	xfftclm	3072K, 0
	xfftclm	3072Kp, 2
	xfftclm	3072Kp, 1
	xfftclm	3072Kp, 0
	xfftclm	3584K, 2
	xfftclm	3584K, 1
	xfftclm	3584K, 0
	xfftclm	4096K, 2
	xfftclm	4096K, 1
	xfftclm	4096K, 0
	xfftclm	4096Kp, 2
	xfftclm	4096Kp, 1
	xfftclm	4096Kp, 0

; Common code to finish the FFT by restoring state and returning.

gw_finish_fft_3:
	mov	DWORD PTR [rsi-28], 3	; Set has-been-FFTed flags
	xfft_1_ret

; Common code to finish up multiplies

; Split the accumulated carries into two carries - a high carry and a
; low carry.  Handle both the with and without two-to-phi array cases.
; Add these carries back into the FFT data.

	PUBLIC	gw_carries_3
gw_carries_3:
	cmp	_ZERO_PADDED_FFT, 0	; Special case the zero padded FFT case
	jne	gw_carries_3_zpad
	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Compute addr after carries
	shl	rbx, 6
	add	rbx, rsi
	xnorm012_2d_part1
	mov	rbp, _DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	sub	rcx, rcx		; Clear big/little flag
	mov	eax, count3		; Load 3 section counts
	mov	loopcount1, eax		; Save for later
ilp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	rbx, 03FFh
	jz	spldn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 10			; Move counts list along
	mov	loopcount1, eax
	shl	rbx, 6			; Compute addr of the last carries row
	add	rbx, rsi
	xnorm012_2d_part2
	sub	rax, rax		; Clear big/little flag
ilp1:	mov	rbx, norm_col_mults	; Addr of the column multipliers
	xnorm012_2d			; Split carries for one cache line
	mov	ebx, cache_line_multiplier; Cache lines in each pass1 loop
	lea	rsi, [rsi+64]		; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	_RATIONAL_FFT, 0	; Don't bump these two pointers
	jne	short iskip		; for rational FFTs
	lea	rdx, [rdx+128]		; Next group multiplier
	lea	rdi, [rdi+rbx*4]	; Next big/little flags pointer
iskip:	sub	loopcount2, 1		; Test loop counter
	jnz	ilp1			; Next carry row in section
	jmp	ilp0			; Next section
spldn:	mov	zero_fft, 0		; Clear zero-high-words-fft flag
	ret

gw_carries_3_zpad:
	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Compute addr after carries
	shl	rbx, 6
	add	rbx, rsi
	mov	rsi, _DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rbp, norm_col_mults	; Addr of the group multipliers
	cmp	const_fft, 0		; Call correct part1 macro
	je	c2a			; Jump if not const
	xnorm012_2d_zpad_part1 exec
	jmp	c2b
c2a:	xnorm012_2d_zpad_part1 noexec
c2b:	mov	rsi, carries		; Addr of the carries
	mov	rbp, _DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load 3 section counts
	mov	loopcount1, eax		; Save for later
zlp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	rbx, 03FFh
	jz	zpldn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 10			; Move counts list along
	mov	loopcount1, eax
	shl	rbx, 6			; Compute addr of the last carries row
	add	rbx, rsi
	xnorm012_2d_zpad_part2
	sub	rax, rax		; Clear big/little flag
zlp1:	mov	rbx, norm_col_mults	; Addr of the column multipliers
	cmp	const_fft, 0		; Call correct zpad macro
	je	c2c			; Jump if not const
	xnorm012_2d_zpad exec		; Split carries for one cache line
	jmp	c2d
c2c:	xnorm012_2d_zpad noexec		; Split carries for one cache line
c2d:	mov	ebx, cache_line_multiplier; Cache lines in each pass1 loop
	lea	rsi, [rsi+64]		; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	_RATIONAL_FFT, 0	; Don't bump these two pointers
	jne	short zskip		; for rational FFTs
	lea	rdx, [rdx+128]		; Next group multiplier
	lea	rdi, [rdi+rbx*4]	; Next big/little flags pointer
zskip:	sub	loopcount2, 1		; Test loop counter
	jnz	zlp1			; Next carry row in section
	jmp	zlp0			; Next section
zpldn:	mov	const_fft, 0		; Clear mul-by-const-fft flag
	ret

; Finish the multiply

gw_finish_mult_3:

; Set FFT-started flag

	mov	rsi, _DESTARG		; Addr of FFT data
	mov	eax, _POSTFFT		; Set FFT started flag
	mov	DWORD PTR [rsi-28], eax

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

	_movsd	xmm7, XMM_SUMOUT	; Add together the two partial sumouts
	addsd	xmm7, XMM_SUMOUT+8
	mulsd	xmm7, ttmp_ff_inv
	_movsd	[rsi-24], xmm7		; Save sum of FFT outputs
	_movsd	xmm6, XMM_MAXERR	; Compute new maximum error
	maxsd	xmm6, XMM_MAXERR+8
	_movsd	_MAXERR, xmm6

; Return

exit:	pop	rsi
	pop	rdi
	pop	rbx
	pop	rbp
	ret

_xmm_gw_ffts3 ENDP

_TEXT	ENDS
END
