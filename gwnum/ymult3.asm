; Copyright 2011 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These routines implement some common cleanup code for yr4dwpn (r4delay with partial normalization) FFTs
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE yarch.mac
INCLUDE ybasics.mac
INCLUDE ymult.mac
INCLUDE ynormal.mac

PUBLIC	ygw_finish_fft
PUBLIC	ygw_finish_mult_wpn

_TEXT SEGMENT

;;*****************************************
;; Routine for finishing off a r4dwpn FFT
;;*****************************************

;; Multiple versions of the carry propagation routines

biglit_incr	EQU	PPTR [rsp+first_local]
grp_incr	EQU	PPTR [rsp+first_local+SZPTR]
saveptr1	EQU	PPTR [rsp+first_local+2*SZPTR]
saveptr2	EQU	PPTR [rsp+first_local+3*SZPTR]
saveptr3	EQU	PPTR [rsp+first_local+4*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+5*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+5*SZPTR+4]

;; Base 2, irrational, no zero padding
PROCFL	ygw_carries_wpn3
	int_prolog 5*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	biglit_incr, rax	; Save biglit increment

	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Load block count
	shl	rbx, 6
	add	rbx, rsi
	ynorm012_wpn_part1 exec

;;BUG	cmp	zero_fft, 1		;; Are we zeroing high words?				;; what to do here?
;;BUG	je	done			;; Yes, skip high words

	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1:	ynorm012_wpn exec, exec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1
	bump	rdx, 2*YMM_GMD		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpn3 ENDP


;; Base 2, rational, no zero padding
PROCFL	ygw_carries_wpnr3
	int_prolog 5*SZPTR+8,0,0

	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Load block count
	shl	rbx, 6
	add	rbx, rsi
	ynorm012_wpn_part1 exec

;;BUG	cmp	zero_fft, 1		;; Are we zeroing high words?				;; what to do here?
;;BUG	je	done			;; Yes, skip high words

	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array			;; Needed???
	mov	rdx, norm_grp_mults	; Addr of the group multipliers				;; Needed???

	mov	eax, addcount1		; Load count of pass 1 blocks
	mov	loopcount2, eax		; Save count
rlp1:	ynorm012_wpn noexec, exec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	sub	loopcount2, 1		; Test loop counter
	jnz	rlp1			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnr3 ENDP


;; Base 2, irrational, zero padding
PROCFL	ygw_carries_wpnzp3
	int_prolog 5*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	biglit_incr, rax	; Save biglit increment

	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Load block count
	shl	rbx, 6
	add	rbx, rsi
	ynorm012_wpn_zpad_part1 exec, exec

	mov	rsi, carries		; Addr of the carries
	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0zp:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1zp:	ynorm012_wpn_zpad exec, exec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1zp
	bump	rdx, YMM_GMD		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0zp			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnzp3 ENDP


;; Base 2, rational, zero padding
PROCFL	ygw_carries_wpnrzp3
	int_prolog 5*SZPTR+8,0,0

	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Load block count
	shl	rbx, 6
	add	rbx, rsi
	ynorm012_wpn_zpad_part1 noexec, exec

	mov	rsi, carries		; Addr of the carries
	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array			;; Needed???
	mov	rdx, norm_grp_mults	; Addr of the group multipliers				;; Needed???

	mov	eax, addcount1		; Load count of pass 1 blocks
	mov	loopcount2, eax		; Save count
rlp1zp:	ynorm012_wpn_zpad noexec, exec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	sub	loopcount2, 1		; Test loop counter
	jnz	rlp1zp			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnrzp3 ENDP


;; Not base 2, irrational, no zero padding
PROCFL	ygw_carries_wpnn3
	int_prolog 5*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	biglit_incr, rax	; Save biglit increment

	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Load block count
	shl	rbx, 6
	add	rbx, rsi
	ynorm012_wpn_part1 noexec

	mov	rsi, carries		; Addr of the carries
	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0n:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1n:	ynorm012_wpn exec, noexec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1n
	bump	rdx, 2*YMM_GMD		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0n			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnn3 ENDP


;; Not base 2, rational, no zero padding
PROCFL	ygw_carries_wpnnr3
	int_prolog 5*SZPTR+8,0,0

	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Load block count
	shl	rbx, 6
	add	rbx, rsi
	ynorm012_wpn_part1 noexec

	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array			;; Needed???
	mov	rdx, norm_grp_mults	; Addr of the group multipliers				;; Needed???

	mov	eax, addcount1		; Load count of pass 1 blocks
	mov	loopcount2, eax		; Save count
rlp1n:	ynorm012_wpn noexec, noexec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	sub	loopcount2, 1		; Test loop counter
	jnz	rlp1n			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnnr3 ENDP


;; Not base 2, irrational, zero padding
PROCFL	ygw_carries_wpnnzp3
	int_prolog 5*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	biglit_incr, rax	; Save biglit increment

	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Load block count
	shl	rbx, 6
	add	rbx, rsi
	ynorm012_wpn_zpad_part1 exec, noexec

	mov	rsi, carries		; Addr of the carries
	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0nzp: mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1nzp: ynorm012_wpn_zpad exec, noexec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1nzp
	bump	rdx, YMM_GMD		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0nzp			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnnzp3 ENDP


;; Not base 2, rational, zero padding
PROCFL	ygw_carries_wpnnrzp3
	int_prolog 5*SZPTR+8,0,0

	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Load block count
	shl	rbx, 6
	add	rbx, rsi
	ynorm012_wpn_zpad_part1 noexec, noexec

	mov	rsi, carries		; Addr of the carries
	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array			;; Needed???
	mov	rdx, norm_grp_mults	; Addr of the group multipliers				;; Needed???

	mov	eax, addcount1		; Load count of pass 1 blocks
	mov	loopcount2, eax		; Save count
rlp1nzp: ynorm012_wpn_zpad noexec, noexec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	sub	loopcount2, 1		; Test loop counter
	jnz	rlp1nzp			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnnrzp3 ENDP

; Common code to finish off r4dwpn FFTs.  The Windows 64-bit ABI
; frowns on us jumping from one procedure into another.
; However, my reading of the spec is that as long as the two procedures have
; identical prologs then stack unwinding for exception handling will work OK.
; Of course, this code won't be linked into a 64-bit Windows executable,
; but we include the dummy prolog to be consistent.

PROCF	__common_wpn_yfft_exit_code

	;; Create a dummy prolog
	ad_prolog 0,1,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

; Common code to finish up ffts

ygw_finish_fft:
	mov	rsi, DESTARG		; Restore source pointer
	mov	DWORD PTR [rsi-28], 3	; Set has-been-FFTed flags
	yfft_1_ret

; Common code to finish up multiplies

; Finish the multiply

ygw_finish_mult_wpn:

; Set FFT-started flag

	mov	rsi, DESTARG		; Addr of FFT data
	mov	al, POSTFFT		; Set FFT started flag
	mov	BYTE PTR [rsi-28], al

; Clear other flags

	mov	zero_fft, 0		; Clear zero-high-words-fft flag for next FFT

; Calculate SUMOUT and MAXERR values

	vmovsd	xmm7, YMM_SUMOUT	; Add together the two partial sumouts
	vaddsd	xmm7, xmm7, YMM_SUMOUT+8
	vaddsd	xmm7, xmm7, YMM_SUMOUT+16
	vaddsd	xmm7, xmm7, YMM_SUMOUT+24
	vmovsd	Q [rsi-24], xmm7	; Save sum of FFT outputs
	vmovsd	xmm6, YMM_MAXERR	; Compute new maximum error
	vmaxsd	xmm6, xmm6, YMM_MAXERR+8
	vmaxsd	xmm6, xmm6, YMM_MAXERR+16
	vmaxsd	xmm6, xmm6, YMM_MAXERR+24
	vmovsd	MAXERR, xmm6

; Return

	ad_epilog 0,1,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

__common_wpn_yfft_exit_code ENDP

_TEXT	ENDS
END
