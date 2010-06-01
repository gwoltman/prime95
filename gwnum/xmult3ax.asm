; Copyright 2001-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu

;
; Additional routines used with r4dwpn (r4delay with partial normalization) FFTs,
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

;;
;; Add two numbers with carry propogation
;;

saved_reg	EQU	PPTR [rsp+first_local]
loopcount1	EQU	DPTR [rsp+first_local+SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+SZPTR+16]

PROCFL	gwxadd3
	ad_prolog SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
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
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address

	;; Do a block

ablk:	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax

add0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2add			; yes, do simpler rounding
nb2add0:cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2iadd0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
	sub	rax, rax		; Clear big/lit flag
nb2radd1:
	xnorm_op_wpn addpd, noexec, noexec, saved_reg ; Add and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2radd1		; Loop til done
	jmp	achunkdn		; Jump to chunk done code
nb2iadd0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2iadd1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax	; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
nb2iadd2:xnorm_op_wpn addpd, exec, noexec, saved_reg ; Add and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2iadd2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2iadd1		; Loop til done
	jmp	achunkdn		; Jump to chunk done code
b2add:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	iadd0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
	sub	rax, rax		; Clear big/lit flag
radd1:	xnorm_op_wpn addpd, noexec, exec, saved_reg ; Add and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	radd1 			; Loop til done
	jmp	achunkdn		; Jump to chunk done code
iadd0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
iadd1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax	; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
iadd2:	xnorm_op_wpn addpd, exec, exec, saved_reg ; Add and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	iadd2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	iadd1			; Loop til done

	;; Chunk done

achunkdn:bump	rcx, 128		; Skip 128 bytes every 8KB
	bump	rdx, 128		; Skip 128 bytes every 8KB
	bump	rsi, 128		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	add0
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest

	;; Block done

ablkdn:	sub	rsi, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_wpn_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	add	rsi, pass1blkdst
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	askip2			; Yes, skip bumping ttp/biglit ptrs
	bump	rbp, 128		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
askip2:	dec	loopcount2		; Decrement outer loop counter
	jnz	ablk 			; Loop til done

	;; Section done

	mov	rax, norm_ptr1		; Reload section start ptr
	mov	rbx, norm_ptr2		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_wpn_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	asec			; Do next section if necessary

	;; All sections done

	mov	rsi, DESTARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	xnorm_op_wpn_fft		; Add 2 carries to start of fft

	ad_epilog SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxadd3 ENDP

;;
;; Subtract two numbers with carry propogation
;;

saved_reg	EQU	PPTR [rsp+first_local]
loopcount1	EQU	DPTR [rsp+first_local+SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+SZPTR+16]

PROCFL	gwxsub3
	ad_prolog SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
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
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address

	;; Do a block

sblk:	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax

sub0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2sub			; yes, do simpler rounding
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2isub0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
	sub	rax, rax		; Clear big/lit flag
nb2rsub1:
	xnorm_op_wpn subpd, noexec, noexec, saved_reg ; Subtract and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2rsub1		; Loop til done
	jmp	schunkdn		; Jump to chunk done code
nb2isub0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2isub1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
nb2isub2:xnorm_op_wpn subpd, exec, noexec, saved_reg ; Subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2isub2		; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2isub1		; Loop til done
	jmp	schunkdn		; Jump to chunk done code
b2sub:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	isub0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
	sub	rax, rax		; Clear big/lit flag
rsub1:	xnorm_op_wpn subpd, noexec, exec, saved_reg ; Subtract and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	rsub1 			; Loop til done
	jmp	schunkdn		; Jump to chunk done code
isub0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
isub1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
isub2:	xnorm_op_wpn subpd, exec, exec, saved_reg ; Subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	isub2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	isub1			; Loop til done

	;; Chunk done

schunkdn:bump	rcx, 128		; Skip 128 bytes every 8KB
	bump	rdx, 128		; Skip 128 bytes every 8KB
	bump	rsi, 128		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	sub0
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest

	;; Block done

sblkdn:	sub	rsi, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_wpn_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	add	rsi, pass1blkdst
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	sskip2			; Yes, skip bumping ttp/biglit ptrs
	bump	rbp, 128		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
sskip2:	dec	loopcount2		; Decrement outer loop counter
	jnz	sblk 			; Loop til done

	;; Section done

	mov	rax, norm_ptr1		; Reload section start ptr
	mov	rbx, norm_ptr2		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_wpn_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	ssec			; Do next section if necessary

	;; All sections done

	mov	rsi, DESTARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	xnorm_op_wpn_fft		; Add 2 carries to start of fft

	ad_epilog SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxsub3 ENDP

;;
;; Add and subtract two numbers with carry propogation
;;

saved_dest1_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_dest2_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
saved_grp_ptr	EQU	PPTR [rsp+first_local+2*SZPTR]
saved_reg	EQU	PPTR [rsp+first_local+3*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+4*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+4*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+4*SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+4*SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+4*SZPTR+16]

PROCFL	gwxaddsub3
	ad_prolog 4*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
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
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	saved_dest1_ptr, rsi	; Save section start address
	mov	saved_dest2_ptr, rbp	; Save section start address
	mov	saved_grp_ptr, rbx	; Save section group ttp address

	;; Do a block

asblk:	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax

as0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2as			; yes, do simpler rounding
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2ias0	   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in an 8KB chunk
	sub	rax, rax		; Clear big/lit flag
nb2ras1: xnorm_addsub_wpn noexec, noexec ; Add & sub and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2ras1			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code
nb2ias0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2ias1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
nb2ias2:xnorm_addsub_wpn exec, noexec	; Add & subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2ias2			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2ias1			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code
b2as:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	ias0	   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in an 8KB chunk
	sub	rax, rax		; Clear big/lit flag
ras1:	xnorm_addsub_wpn noexec, exec	; Add & sub and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	ras1 			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code
ias0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
ias1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flag
ias2:	xnorm_addsub_wpn exec, exec	; Add & subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	ias2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	ias1			; Loop til done

	;; Chunk done

aschunkdn:bump	rcx, 128		; Skip 128 bytes every 8KB
	bump	rdx, 128		; Skip 128 bytes every 8KB
	bump	rsi, 128		; Skip 128 bytes every 8KB
	bump	rbp, 128		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	as0
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	add	rbp, pass2gapsize	; Next dest

	;; Block done

asblkdn:sub	rsi, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_wpn_blk rsi, rbx, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	sub	rbp, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_wpn_blk rbp, rbx, XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8
	add	rsi, pass1blkdst
	add	rbp, pass1blkdst
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	asskip2			; Yes, skip bumping ttp/biglit ptrs
	bump	rbx, 128		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
asskip2:dec	loopcount2		; Decrement outer loop counter
	jnz	asblk 			; Loop til done

	;; Section done

	mov	saved_reg, rbx
	mov	rax, saved_dest1_ptr	; Reload section start ptr for dest #1
	mov	rbx, saved_grp_ptr	; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_wpn_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	rax, saved_dest2_ptr	; Reload section start ptr for dest #2
					; Add 2 carries to start of section
	xnorm_op_wpn_sec XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8
	mov	rbx, saved_reg
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	assec			; Do next section if necessary

	;; All sections done

	mov	rsi, DESTARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	xnorm_op_wpn_fft		; Add 2 carries to start of fft

	mov	rsi, DEST2ARG		; Addr of FFT data
	movapd	xmm7, XMM_TMP7		; Load wraparound carry
	xnorm_top_carry_cmn rsi, xmm7, 2
	movapd	xmm6, XMM_TMP5		; Load non-wraparound carry
	xnorm_op_wpn_fft		; Add 2 carries to start of fft
	ad_epilog 4*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxaddsub3 ENDP

;;
;; Add in a small number with carry propogation
;;

PROCFL	gwxadds3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	movsd	xmm7, DBLARG		; Small addin value

	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	sub	rax, rax		; Clear biglit flag
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2addsm			; yes, do simpler rounding
	xnorm_smalladd_wpn noexec	; Similar to add last carry code
	jmp	addsmdn
b2addsm:xnorm_smalladd_wpn exec		; Similar to add last carry code
addsmdn:
	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxadds3 ENDP

;;
;; Multiply a number by a small value with carry propogation
;;

saved_sec_biglit EQU	PPTR [rsp+first_local+0*SZPTR]
saved_blk_start EQU	PPTR [rsp+first_local+1*SZPTR]
saved_blk_biglit EQU	PPTR [rsp+first_local+2*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+3*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+3*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+3*SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+3*SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+3*SZPTR+16]

PROCFL	gwxmuls3
	ad_prolog 3*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	movlpd	xmm0, DBLARG		; Small multiplier value
	movhpd	xmm0, DBLARG		; Small multiplier value
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	jne	skip			; No, skip mul by two-to-phi fudge factor
skip:	movapd	XMM_TMP5, xmm0		; Save small value
	movapd	xmm7, XMM_BIGVAL	; Init 4 carries
	movapd	XMM_TMP1, xmm7
	movapd	XMM_TMP2, xmm7
	movapd	XMM_TMP3, xmm7
	movapd	XMM_TMP4, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

msec:	mov	loopcount1, eax		; Save section counts
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address
	mov	saved_sec_biglit, rdi	; Save section biglit ptr

	;; Do a block

mblk:	mov	saved_blk_start, rsi
	mov	saved_blk_biglit, rdi
	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax

mul0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2mul			; yes, do simpler rounding
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2imul0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
	sub	rax, rax		; Clear big/lit flags
	sub	rcx, rcx
nb2rmul1:xnorm_smallmul_wpn noexec, noexec ; Mul and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2rmul1		; Loop til done
	jmp	mchunkdn		; Jump to chunk done code
nb2imul0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2imul1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flags
	sub	rcx, rcx
nb2imul2:xnorm_smallmul_wpn exec, noexec ; Mul and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2imul2		; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2imul1		; Loop til done
	jmp	mchunkdn		; Jump to chunk done code

b2mul:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	imul0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
	sub	rax, rax		; Clear big/lit flags
	sub	rcx, rcx
rmul1:	xnorm_smallmul_wpn noexec, exec	; Mul and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	rmul1 			; Loop til done
	jmp	mchunkdn		; Jump to chunk done code
imul0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
imul1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	sub	rax, rax		; Clear big/lit flags
	sub	rcx, rcx
imul2:	xnorm_smallmul_wpn exec, exec	; Mul and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	imul2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	imul1			; Loop til done

	;; Chunk done

mchunkdn:bump	rsi, 128		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	mul0

	;; Block done

	xchg	rsi, saved_blk_start	; Restore/save block start ptr
	xchg	rdi, saved_blk_biglit	; Restore/save block biglit ptr
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2mblk			; yes, do simpler rounding
	xnorm_smallmul_wpn_blk noexec	; Add 4 carries to start of block
	jmp	mblkdn
b2mblk:	xnorm_smallmul_wpn_blk exec	; Add 4 carries to start of block
mblkdn:	mov	rsi, saved_blk_start	; Restore start ptr
	mov	rdi, saved_blk_biglit	; Restore biglit ptr
	add	rsi, pass2gapsize	; Next dest
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	mskip2			; Yes, skip bumping ttp/biglit ptrs
	add	rdi, normval3		; Adjust little/big flags ptr
	bump	rbp, 128		; Next set of group multipliers
mskip2:	dec	loopcount2		; Decrement outer loop counter
	jnz	mblk 			; Loop til section done

	;; Section done

	xchg	rsi, norm_ptr1		; Restore/save section start ptr
	xchg	rbp, norm_ptr2		; Restore/save section group ttp ptr
	xchg	rdi, saved_sec_biglit	; Restore/save section biglit ptr
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2msec			; yes, do simpler rounding
	xnorm_smallmul_wpn_sec noexec	; Add 2 carries to start of section
	jmp	msecdn
b2msec:	xnorm_smallmul_wpn_sec exec	; Add 2 carries to start of section
msecdn:	mov	rsi, norm_ptr1		; Restore next section start ptr
	mov	rbp, norm_ptr2		; Restore next section group ttp ptr
	mov	rdi, saved_sec_biglit	; Restore/save section biglit ptr
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	msec			; Do next section if necessary

	;; All sections done

	mov	rsi, DESTARG		; Addr of FFT data
	movapd	xmm6, XMM_TMP1		; Load non-wraparound carry
	movapd	xmm7, XMM_TMP3		; Load wraparound carry
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array

	cmp	ZERO_PADDED_FFT, 0	;; Zero-padded FFT?
	jne	mulzp			;; Yes, do special zpad carry

	xnorm_top_carry_cmn rsi, xmm7, 2
	sub	rax, rax		; Trashed by xnorm_top_carry_cmn
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2mfft			; yes, do simpler rounding
	xnorm_smallmul_wpn_fft noexec	; Add 2 carries to start of fft
	jmp	muldn
b2mfft:	xnorm_smallmul_wpn_fft exec	; Add 2 carries to start of fft
	jmp	muldn

mulzp:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2mzpfft		; yes, do simpler rounding
	xnorm_smallmul_wpn_fft_zpad noexec ; Do the special zpad carry not base 2
	jmp	muldn
b2mzpfft:
	xnorm_smallmul_wpn_fft_zpad exec ; Do the special zpad carry base 2
muldn:

	ad_epilog 3*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxmuls3 ENDP

_TEXT	ENDS
END
