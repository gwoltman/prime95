; Copyright 2001-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
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
INCLUDE memory.mac
INCLUDE xnormal.mac

_TEXT SEGMENT

blkrows	EQU	1024

gwx3procs PROC

;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwxaddf3
gwxaddf3:
	PUBLIC	gwxaddq3
gwxaddq3:
	push	rsi			; Save registers
	push	rdi			; Save registers
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	eax, addcount1		; Load loop counter
uadd1:	mov	edi, blkrows		; Cache lines in a block
uaddlp:	movapd	xmm0, [rdx]		; Load second number
	addpd	xmm0, [rcx]		; Add in first number
	movapd	xmm1, [rdx+16]		; Load second number
	addpd	xmm1, [rcx+16]		; Add in first number
	movapd	xmm2, [rdx+32]		; Load second number
	addpd	xmm2, [rcx+32]		; Add in first number
	movapd	xmm3, [rdx+48]		; Load second number
	addpd	xmm3, [rcx+48]		; Add in first number
	movapd	[rsi], xmm0		; Save result
	movapd	[rsi+16], xmm1		; Save result
	movapd	[rsi+32], xmm2		; Save result
	movapd	[rsi+48], xmm3		; Save result
	lea	rcx, [rcx+64]		; Next source
	lea	rdx, [rdx+64]		; Next source
	lea	rsi, [rsi+64]		; Next dest
	sub	edi, 1			; Test inner counter
	jnz	short uaddlp		; Loop if necessary
	lea	rcx, [rcx-blkrows*64]	; Next source
	lea	rdx, [rdx-blkrows*64]	; Next source
	lea	rsi, [rsi-blkrows*64]	; Next dest
	add	rcx, pass1blkdst	; Next source
	add	rdx, pass1blkdst	; Next source
	add	rsi, pass1blkdst	; Next dest
	sub	eax, 1			; Check loop counter
	jnz	short uadd1		; Loop if necessary
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret

;;
;; Add two numbers with carry propogation
;;

	PUBLIC	gwxadd3
gwxadd3:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbp			; Save registers
	push	rbx			; Save registers
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	movapd	xmm7, XMM_BIGVAL	; Init 4 carries
	movapd	XMM_TMP1, xmm7
	movapd	XMM_TMP2, xmm7
	movapd	XMM_TMP3, xmm7
	movapd	XMM_TMP4, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

asec:	mov	loopcount1, eax		; Save section counts
	and	eax, 03FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address

	;; Do a block

ablk:	cmp	_RATIONAL_FFT, 0	; Test for irrational FFTs
	je	iadd0   		; Yes, use two-to-phi multipliers

	mov	loopcount3, blkrows	; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
radd1:	xnorm_op_2d addpd, noexec	; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement inner loop counter
	jnz	radd1 			; Loop til done
	jmp	ablkdn			; Jump to block done code

iadd0:	mov	eax, normval1		; Load middle count (clms in a blk)
	mov	loopcount3, eax		; Save middle count
	mov	rbx, norm_col_mults	; Addr of the column multipliers
iadd1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount4, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
iadd2:	xnorm_op_2d addpd, exec		; Add and normalize 8 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	iadd2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	iadd1			; Loop til done

	;; Block done

ablkdn:	lea	rsi, [rsi-blkrows*64]	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	lea	rcx, [rcx-blkrows*64]
	lea	rdx, [rdx-blkrows*64]
	add	rcx, pass1blkdst	; Next source
	add	rdx, pass1blkdst	; Next source
	add	rsi, pass1blkdst	; Next dest
	lea	rbp, [rbp+128]		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	ablk 			; Loop til done

	;; Section done

	mov	rax, norm_ptr1		; Reload section start ptr
	mov	rbx, norm_ptr2		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 10			; Move counts list along
	jnz	asec			; Do next section if necessary

	;; All sections done

	mov	rsi, _DESTARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	xnorm_op_2d_fft			; Add 2 carries to start of fft

	pop	rbx			; Restore registers
	pop	rbp			; Restore registers
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret



;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

	PUBLIC	gwxsubf3
gwxsubf3:
	PUBLIC	gwxsubq3
gwxsubq3:
	push	rsi			; Save registers
	push	rdi			; Save registers
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	eax, addcount1		; Load loop counter
usub1:	mov	edi, blkrows		; Cache lines in a block
usublp:	movapd	xmm0, [rdx]		; Load second number
	subpd	xmm0, [rcx]		; Subtract first number
	movapd	xmm1, [rdx+16]		; Load second number
	subpd	xmm1, [rcx+16]		; Subtract first number
	movapd	xmm2, [rdx+32]		; Load second number
	subpd	xmm2, [rcx+32]		; Subtract first number
	movapd	xmm3, [rdx+48]		; Load second number
	subpd	xmm3, [rcx+48]		; Subtract first number
	movapd	[rsi], xmm0		; Save result
	movapd	[rsi+16], xmm1		; Save result
	movapd	[rsi+32], xmm2		; Save result
	movapd	[rsi+48], xmm3		; Save result
	lea	rcx, [rcx+64]		; Next source
	lea	rdx, [rdx+64]		; Next source
	lea	rsi, [rsi+64]		; Next dest
	sub	edi, 1			; Test inner counter
	jnz	short usublp		; Loop if necessary
	lea	rcx, [rcx-blkrows*64]	; Next source
	lea	rdx, [rdx-blkrows*64]	; Next source
	lea	rsi, [rsi-blkrows*64]	; Next dest
	add	rcx, pass1blkdst	; Next source
	add	rdx, pass1blkdst	; Next source
	add	rsi, pass1blkdst	; Next dest
	sub	eax, 1			; Check loop counter
	jnz	short usub1		; Loop if necessary
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret

;;
;; Subtract two numbers with carry propogation
;;

	PUBLIC	gwxsub3
gwxsub3:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbp			; Save registers
	push	rbx			; Save registers
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	movapd	xmm7, XMM_BIGVAL	; Init 4 carries
	movapd	XMM_TMP1, xmm7
	movapd	XMM_TMP2, xmm7
	movapd	XMM_TMP3, xmm7
	movapd	XMM_TMP4, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

ssec:	mov	loopcount1, eax		; Save section counts
	and	eax, 03FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address

	;; Do a block

sblk:	cmp	_RATIONAL_FFT, 0	; Test for irrational FFTs
	je	isub0   		; Yes, use two-to-phi multipliers

	mov	loopcount3, blkrows	; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
rsub1:	xnorm_op_2d subpd, noexec	; Subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement inner loop counter
	jnz	rsub1 			; Loop til done
	jmp	sblkdn			; Jump to block done code

isub0:	mov	eax, normval1		; Load middle count (clms in a blk)
	mov	loopcount3, eax		; Save middle count
	mov	rbx, norm_col_mults	; Addr of the column multipliers
isub1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount4, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
isub2:	xnorm_op_2d subpd, exec		; Subtract and normalize 8 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	isub2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	isub1			; Loop til done

	;; Block done

sblkdn:	lea	rsi, [rsi-blkrows*64]	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	lea	rcx, [rcx-blkrows*64]
	lea	rdx, [rdx-blkrows*64]
	add	rcx, pass1blkdst	; Next source
	add	rdx, pass1blkdst	; Next source
	add	rsi, pass1blkdst	; Next dest
	lea	rbp, [rbp+128]		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	sblk 			; Loop til done

	;; Section done

	mov	rax, norm_ptr1		; Reload section start ptr
	mov	rbx, norm_ptr2		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 10			; Move counts list along
	jnz	ssec			; Do next section if necessary

	;; All sections done

	mov	rsi, _DESTARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	xnorm_op_2d_fft			; Add 2 carries to start of fft

	pop	rbx			; Restore registers
	pop	rbp			; Restore registers
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret


;;
;; Add and subtract two numbers without carry propogation.
;;

	PUBLIC	gwxaddsubf3
gwxaddsubf3:
	PUBLIC	gwxaddsubq3
gwxaddsubq3:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbp			; Save registers
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination #1
	mov	rbp, _DEST2ARG  	; Address of destination #2
	mov	eax, addcount1		; Load loop counter
uaddsub1:mov	edi, blkrows		; Cache lines in a block
uaddsublp:
	movapd	xmm0, [rcx]		; Load first number
	movapd	xmm1, xmm0		; Dup first number
	addpd	xmm0, [rdx]		; Add in second number
	subpd	xmm1, [rdx]		; Subtract out second number
	movapd	xmm2, [rcx+16]		; Load first number
	movapd	xmm3, xmm2		; Dup first number
	addpd	xmm2, [rdx+16]		; Add in second number
	subpd	xmm3, [rdx+16]		; Subtract out second number
	movapd	xmm4, [rcx+32]		; Load first number
	movapd	xmm5, xmm4		; Dup first number
	addpd	xmm4, [rdx+32]		; Add in second number
	subpd	xmm5, [rdx+32]		; Subtract out second number
	movapd	xmm6, [rcx+48]		; Load first number
	movapd	xmm7, xmm6		; Dup first number
	addpd	xmm6, [rdx+48]		; Add in second number
	subpd	xmm7, [rdx+48]		; Subtract out second number
	movapd	[rsi], xmm0		; Save result
	movapd	[rbp], xmm1		; Save result
	movapd	[rsi+16], xmm2		; Save result
	movapd	[rbp+16], xmm3		; Save result
	movapd	[rsi+32], xmm4		; Save result
	movapd	[rbp+32], xmm5		; Save result
	movapd	[rsi+48], xmm6		; Save result
	movapd	[rbp+48], xmm7		; Save result
	lea	rcx, [rcx+64]		; Next source
	lea	rdx, [rdx+64]		; Next source
	lea	rsi, [rsi+64]		; Next dest
	lea	rbp, [rbp+64]		; Next dest
	sub	edi, 1			; Test inner counter
	jnz	uaddsublp		; Loop if necessary
	lea	rcx, [rcx-blkrows*64]	; Next source
	lea	rdx, [rdx-blkrows*64]	; Next source
	lea	rsi, [rsi-blkrows*64]	; Next dest
	lea	rbp, [rbp-blkrows*64]	; Next dest
	add	rcx, pass1blkdst	; Next source
	add	rdx, pass1blkdst	; Next source
	add	rsi, pass1blkdst	; Next dest
	add	rbp, pass1blkdst	; Next dest
	sub	eax, 1			; Check loop counter
	jnz	uaddsub1		; Loop if necessary
	pop	rbp			; Restore registers
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret

;;
;; Add and subtract two numbers with carry propogation
;;

	PUBLIC	gwxaddsub3
gwxaddsub3:
	push	rsi			; Save registers
	push	rdi			; Save registers
	push	rbp			; Save registers
	push	rbx			; Save registers
	mov	rcx, _SRCARG		; Address of first number
	mov	rdx, _SRC2ARG		; Address of second number
	mov	rsi, _DESTARG		; Address of destination
	mov	rbp, _DEST2ARG  	; Address of destination #2
	mov	rbx, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	movapd	xmm7, XMM_BIGVAL	; Init 8 carries
	movapd	XMM_TMP1, xmm7
	movapd	XMM_TMP2, xmm7
	movapd	XMM_TMP3, xmm7
	movapd	XMM_TMP4, xmm7
	movapd	XMM_TMP5, xmm7
	movapd	XMM_TMP6, xmm7
	movapd	XMM_TMP7, xmm7
	movapd	XMM_TMP8, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

assec:	mov	loopcount1, eax		; Save section counts
	and	eax, 03FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section start address
	mov	norm_ptr3, rbx		; Save section group ttp address

	;; Do a block

asblk:	cmp	_RATIONAL_FFT, 0	; Test for irrational FFTs
	je	ias0	   		; Yes, use two-to-phi multipliers

	mov	loopcount3, blkrows	; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
ras1:	xnorm_addsub_2d noexec		; Add & subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement inner loop counter
	jnz	ras1 			; Loop til done
	jmp	asblkdn			; Jump to block done code

ias0:	mov	eax, normval1		; Load middle count (clms in a blk)
	mov	loopcount3, eax		; Save middle count
	mov	rax, norm_col_mults	; Addr of the column multipliers
ias1:	push	rax
	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount4, eax		; Save inner loop count
	pop	rax
ias2:	xnorm_addsub_2d exec		; Add & subtract and normalize 8 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	ias2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	ias1			; Loop til done

	;; Block done

asblkdn:lea	rsi, [rsi-blkrows*64]	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rsi, rbx, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	lea	rbp, [rbp-blkrows*64]	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rbp, rbx, XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8
	lea	rcx, [rcx-blkrows*64]
	lea	rdx, [rdx-blkrows*64]
	add	rcx, pass1blkdst	; Next source
	add	rdx, pass1blkdst	; Next source
	add	rsi, pass1blkdst	; Next dest
	add	rbp, pass1blkdst	; Next dest
	lea	rbx, [rbx+128]		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	asblk 			; Loop til done

	;; Section done

	push	rbx
	mov	rax, norm_ptr1		; Reload section start ptr for dest #1
	mov	rbx, norm_ptr3		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	rax, norm_ptr2		; Reload section start ptr for dest #2
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8
	pop	rbx
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 10			; Move counts list along
	jnz	assec			; Do next section if necessary

	;; All sections done

	mov	rsi, _DESTARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	xnorm_op_2d_fft			; Add 2 carries to start of fft

	mov	rsi, _DEST2ARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP7		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	movapd	xmm6, XMM_TMP5		; Load non-wraparound carry
	xnorm_op_2d_fft			; Add 2 carries to start of fft

	pop	rbx			; Restore registers
	pop	rbp			; Restore registers
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret



gwx3procs ENDP

;;
;; Copy one number and zero some low order words.
;;

	PUBLIC	gwxcopyzero3
gwxcopyzero3 PROC
	push	rsi			; Save registers
	push	rdi			; Save registers
	mov	rsi, _SRCARG		; Address of first number
	mov	rdi, _DESTARG		; Address of destination
	sub	ecx, ecx		; Offset to compare to COPYZERO
	mov	eax, addcount1		; Get number of blocks
cz1:	mov	edx, blkrows		; Number of cache lines in one block
cz2:	xcopyzero
	lea	rsi, [rsi+64]		; Next source
	lea	rdi, [rdi+64]		; Next dest
	lea	rcx, [rcx+64]		; Next compare offset
	sub	edx, 1			; Test loop counter
	jnz	cz2			; Loop if necessary
	lea	rsi, [rsi-blkrows*64]
	lea	rdi, [rdi-blkrows*64]
	lea	rcx, [rcx-blkrows*64]
	add	rsi, pass1blkdst
	add	rdi, pass1blkdst
	add	rcx, pass1blkdst
	sub	eax, 1			; Test loop counter
	jnz	cz1			; Loop if necessary
	pop	rdi			; Restore registers
	pop	rsi			; Restore registers
	ret
gwxcopyzero3 ENDP

;;
;; Routines to do the normalization after a multiply
;;

_xnorm3 PROC

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	noadd, ilp0, ilp1, ilexit
	PUBLIC	lab
lab:
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
ilp1:	xnorm_2d ttp, zero, echk, const	;; Normalize 8 values
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
	PUBLIC	lab
lab:
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
ilp1:	xnorm_2d_zpad ttp, echk, const	;; Normalize 8 values
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

_xnorm3 ENDP

_TEXT	ENDS
END
