; Copyright 1995-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;

INCLUDE pfa.mac
INCLUDE lucas.mac
INCLUDE lucasp.mac
INCLUDE xmult.mac
INCLUDE xlucas.mac
INCLUDE xnormal.mac
INCLUDE xpass2.mac

IFDEF X86_64
X87_CASES	EQU	0
ELSE
X87_CASES	EQU	13
ENDIF
SSE2_CASES	EQU	120

loopent	MACRO	y,z		; Create a entry in the loop entry table
	DPTR	&y&z
	ENDM
looptab	MACRO	y, cnt		; Create the loop entry table
	x = 0
	REPT	cnt
	loopent	y, %x
	x = x + 1
	ENDM
	ENDM

x87table: looptab case, X87_CASES
sse2table: looptab sscase, SSE2_CASES

;; Macros from p4notes.doc

read1	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
	cnt = mem/64		        ;; Read 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:	mov     eax, [rsi]		;; Read one cache line
	mov     eax, [rsi+4]		;; 4 bytes at a time
	mov     eax, [rsi+8]
	mov     eax, [rsi+12]
	mov     eax, [rsi+16]
	mov     eax, [rsi+20]
	mov     eax, [rsi+24]
	mov     eax, [rsi+28]
	mov     eax, [rsi+32]
	mov     eax, [rsi+36]
	mov     eax, [rsi+40]
	mov     eax, [rsi+44]
	mov     eax, [rsi+48]
	mov     eax, [rsi+52]
	mov     eax, [rsi+56]
	mov     eax, [rsi+60]
	lea     rsi, [rsi+64]			; Next cache line
	sub     ecx, 1
	jnz     loop1
	lea     rsi, [rsi-mem]		; Restore esi
	dec	edx
	jnz	loop2
	ENDM

read2	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
        cnt = mem/64                    ; 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  xmm1, [rsi]			; Read one cache line
        movapd  xmm1, [rsi+16]		; 16 bytes at a time
        movapd  xmm1, [rsi+32]
        movapd  xmm1, [rsi+48]
        lea     rsi, [rsi+64]			; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

write1	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/64                    ; 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  [rsi], xmm1			; Write one cache line
        movapd  [rsi+16], xmm1		; 16 bytes at a time
        movapd  [rsi+32], xmm1
        movapd  [rsi+48], xmm1
        lea     rsi, [rsi+64]			; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

write2	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/(4*128)               ; 128 bytes per iteration
        dist = 64
	mov	edx, c
	sub     ebx, ebx
loop2:	mov     ecx, cnt
loop1:  movapd  [rsi+0*dist], xmm1      ; Write 8 cache lines
        movapd  [rsi+1*dist], xmm1
        movapd  [rsi+2*dist], xmm1
        movapd  [rsi+3*dist], xmm1
        movapd  [rsi+4*dist], xmm1
        movapd  [rsi+5*dist], xmm1
        movapd  [rsi+6*dist], xmm1
        movapd  [rsi+7*dist], xmm1
        lea     rsi, [rsi+16]           ; Same cache lines
        add     bl, 256/4			; 4 inner loop iterations
        jnc     loop1
        lea     rsi, [rsi-4*16+8*dist]  ; Next set of 8 cache lines
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

readwrite1 MACRO mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/64
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  xmm0, [rsi]             ; Read one cache line
        movapd  xmm1, [rsi+16]
        movapd  xmm2, [rsi+32]
        movapd  xmm3, [rsi+48]
        subpd   xmm0, xmm0              ; Operate on the data
        pxor    xmm1, xmm1
        subpd   xmm2, xmm2
        pxor    xmm3, xmm3
        movapd  [rsi], xmm0             ; Write the cache line
        movapd  [rsi+16], xmm1
        movapd  [rsi+32], xmm2
        movapd  [rsi+48], xmm3
        lea     rsi, [rsi+64]           ; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

;; Time one of the basic FFT building blocks

x87mac	MACRO	memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / memused
	outer_iters = 10000 / inner_iters
	mov	eax, outer_iters
ss0a:	mov	rdi, _SRCARG
	lea	rsi, [rdi+4096]
	mov	ecx, inner_iters
ss0b:	disp	&ops
	lea	rsi, [rsi+memused]	;; Next source pointer
	lea	rdi, [rdi+SCD]		;; Next sine/cosine pointer
	dec	ecx
	jnz	ss0b
	dec	eax
	jnz	ss0a
	ENDM

sse2mac	MACRO	memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / memused
	outer_iters = 10000 / inner_iters
	mov	rbx, 262144		;; Offset for some sse2 macros
	mov	rbp, 524288		;; Offset for mulf sse2 macros
	mov	eax, outer_iters
ss0a:	mov	rdi, _SRCARG		;; Sincos data
	lea	rsi, [rdi+4096]		;; Source & dest ptr
	lea	rdx, [rsi+524288+256]	;; Destination for "g" macros
	mov	ecx, inner_iters
ss0b:	&ops
	lea	rdi, [rdi+2*XMM_SCD]	;; Next sine/cosine pointer
	dec	ecx
	jnz	ss0b
	dec	eax
	jnz	ss0a
	ENDM


PUBLIC _timeit
_timeit	PROC
	push	rsi
	push	rdi
	push	rbp
	push	rbx

	clear_timers

	sub	ebx, ebx
	sub	ecx, ecx
	sub	ebp, ebp
	mov	rdi, _SRCARG
	mov	rsi, _SRCARG

	mov	edx, _NUMARG		; Get test number

	mov	eax, X87_CASES
	cmp	edx, -1			; -1 = get num x87 cases
	je	exit
	mov	eax, SSE2_CASES
	cmp	edx, -2			; -2 = get num sse2 cases
	je	exit

	cmp	edx, 1000		; Tests above 1000 are SSE2 code
	jl	short x87
	subpd	xmm0, xmm0
	subpd	xmm1, xmm1
	subpd	xmm2, xmm2
	subpd	xmm3, xmm3
	subpd	xmm4, xmm4
	subpd	xmm5, xmm5
	subpd	xmm6, xmm6
	subpd	xmm7, xmm7

; Jump to desired test case

	sub	rdx, 1000
	mov	rax, sse2table[rdx*SZPTR]; Get address of test to execute
	jmp	rax
x87:	mov	rax, x87table[rdx*SZPTR]; Get address of test to execute
	jmp	rax

; Time the "do-nothing" case.

IFNDEF X86_64
case0:
	jmp	exit

; Time the loop of "do-nothing" case.  1000 iterations.

case1:
	mov	ecx, 1000
c1a:	dec	ecx
	jnz	c1a
	jmp	exit

; This code reads a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will read from the L1 cache
; only, 96KB will read from the L2 cache only, and 2MB will test reading
; from main memory.

case2:	read1	4096, 1000	; Read 4KB
	jmp	exit

case3:	read1	96*1024, 100	; Read 96KB
	jmp	exit

case4:	read1	2048*1024, 2	; Read 2MB
	jmp	exit

case5:	x87mac	64, 4096, eight_reals_fft, 8, 16, 32
	jmp	exit

case6:	x87mac	64, 100000, eight_reals_fft, 8, 16, 32
	jmp	exit

case7:	x87mac	64, 4096, eight_reals_unfft, 8, 16, 32
	jmp	exit

case8:	x87mac	64, 100000, eight_reals_unfft, 8, 16, 32
	jmp	exit

case9:	x87mac	64, 4096, four_complex_fft, 8, 16, 32
	jmp	exit

case10:	x87mac	64, 100000, four_complex_fft, 8, 16, 32
	jmp	exit

case11:	x87mac	64, 4096, four_complex_unfft, 8, 16, 32
	jmp	exit

case12:	x87mac	64, 100000, four_complex_unfft, 8, 16, 32
	jmp	exit
ENDIF

; This code reads a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will read from the L1 cache
; only, 96KB will read from the L2 cache only, and 2MB will test reading
; from main memory.

sscase0:
	read2	4096, 1000	; Read 4KB
	jmp	exit

sscase1:
	read2	96*1024, 100	; Read 96KB
	jmp	exit

sscase2:
	read2	2048*1024, 2	; Read 2MB
	jmp	exit

; This code writes a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 2MB will test writing
; to main memory.

sscase3:
	write1	4096, 1000	; Read 4KB
	jmp	exit

sscase4:
	write1	96*1024, 100	; Read 96KB
	jmp	exit

sscase5:
	write1	2048*1024, 2	; Read 2MB
	jmp	exit

; This code writes a block of memory non-contiguously.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 2MB will test writing
; to main memory.

sscase6:
	write2	4096, 1000	; Read 4KB
	jmp	exit

sscase7:
	write2	96*1024, 100	; Read 96KB
	jmp	exit

sscase8:
	write2	2048*1024, 2	; Read 2MB
	jmp	exit

; This code reads & writes a block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 2MB will test writing
; to main memory.

sscase9:
	readwrite1	4096, 1000	; Read 4KB
	jmp	exit

sscase10:
	readwrite1	96*1024, 100	; Read 96KB
	jmp	exit

sscase11:
	readwrite1	2048*1024, 2	; Read 2MB
	jmp	exit

; Time ~10000 iterations of the SSE2 macros in L1 and L2 caches

sscase12:
	sse2mac	128, 4096, x2cl_eight_reals_fft rsi, 2*64, 64
	jmp	exit
sscase13:
	sse2mac	128, 100000, x2cl_eight_reals_fft rsi, 2*64, 64
	jmp	exit
sscase14:
	sse2mac	128, 4096, x2cl_eight_reals_first_fft rsi, 2*64, 64
	jmp	exit
sscase15:
	sse2mac	128, 100000, x2cl_eight_reals_first_fft rsi, 2*64, 64
	jmp	exit
sscase16:
	sse2mac	128, 4096, x2cl_eight_reals_fft_2 rsi, 2*64, 64
	jmp	exit
sscase17:
	sse2mac	128, 100000, x2cl_eight_reals_fft_2 rsi, 2*64, 64
	jmp	exit
sscase18:
	sse2mac	128, 4096, x2cl_eight_reals_fft_1 rsi, 2*64, 64
	jmp	exit
sscase19:
	sse2mac	128, 100000, x2cl_eight_reals_fft_1 rsi, 2*64, 64
	jmp	exit
sscase20:
	sse2mac	128, 4096, s2cl_eight_reals_first_fft rsi, 2*64, 64
	jmp	exit
sscase21:
	sse2mac	128, 100000, s2cl_eight_reals_first_fft rsi, 2*64, 64
	jmp	exit
sscase22:
	sse2mac	128, 4096, s2cl_eight_reals_fft_1 rsi, 2*64, 64
	jmp	exit
sscase23:
	sse2mac	128, 100000, s2cl_eight_reals_fft_1 rsi, 2*64, 64
	jmp	exit
sscase24:
	sse2mac	128, 4096, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	jmp	exit
sscase25:
	sse2mac	128, 100000, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	jmp	exit
sscase26:
	sse2mac	128, 4096, s2cl_eight_reals_fft_2_final rsi, 2*64, 64
	jmp	exit
sscase27:
	sse2mac	128, 100000, s2cl_eight_reals_fft_2_final rsi, 2*64, 64
	jmp	exit
sscase28:
	sse2mac	128, 4096, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	jmp	exit
sscase29:
	sse2mac	128, 100000, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	jmp	exit
sscase30:
	sse2mac	128, 4096, s2cl_eight_reals_with_mult_2 rsi, 2*64, 64
	jmp	exit
sscase31:
	sse2mac	128, 100000, s2cl_eight_reals_with_mult_2 rsi, 2*64, 64
	jmp	exit
sscase32:
	sse2mac	128, 4096, s2cl_eight_reals_with_mulf_2 rsi, 2*64, 64
	jmp	exit
sscase33:
	sse2mac	128, 100000, s2cl_eight_reals_with_mulf_2 rsi, 2*64, 64
	jmp	exit
sscase34:
	sse2mac	256, 4096, x4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64
	jmp	exit
sscase35:
	sse2mac	256, 100000, x4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64
	jmp	exit
sscase36:
	sse2mac	256, 4096, g4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	jmp	exit
sscase37:
	sse2mac	256, 100000, g4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	jmp	exit

sscase38:
	sse2mac	128, 4096, x2cl_eight_reals_unfft_2 rsi, 2*64, 64
	jmp	exit
sscase39:
	sse2mac	128, 100000, x2cl_eight_reals_unfft_2 rsi, 2*64, 64
	jmp	exit
sscase40:
	sse2mac	128, 4096, x2cl_half_eight_reals_unfft_2 rsi, 2*64, 64
	jmp	exit
sscase41:
	sse2mac	128, 100000, x2cl_half_eight_reals_unfft_2 rsi, 2*64, 64
	jmp	exit
sscase42:
	sse2mac	128, 4096, g2cl_eight_reals_unfft_2 rsi, 2*64, 64, rdx, 2*64, 64
	jmp	exit
sscase43:
	sse2mac	128, 100000, g2cl_eight_reals_unfft_2 rsi, 2*64, 64, rdx, 2*64, 64
	jmp	exit
sscase44:
	sse2mac	256, 4096, x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase45:
	sse2mac	256, 100000, x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase46:
	sse2mac	256, 4096, x4cl_eight_reals_unfft_2 rsi, 4*64, 64, 2*64
	jmp	exit
sscase47:
	sse2mac	256, 100000, x4cl_eight_reals_unfft_2 rsi, 4*64, 64, 2*64
	jmp	exit
sscase48:
	sse2mac	256, 4096, s4cl_eight_reals_unfft_1 rsi, 4*64, 64, 2*64
	jmp	exit
sscase49:
	sse2mac	256, 100000, s4cl_eight_reals_unfft_1 rsi, 4*64, 64, 2*64
	jmp	exit

sscase50:
	sse2mac	128, 4096, x2cl_two_complex_fft rsi, 2*64, 64
	jmp	exit
sscase51:
	sse2mac	128, 100000, x2cl_two_complex_fft rsi, 2*64, 64
	jmp	exit
sscase52:
	sse2mac	128, 4096, x2cl_two_complex_fft_in_place rsi, 2*64, 64
	jmp	exit
sscase53:
	sse2mac	128, 100000, x2cl_two_complex_fft_in_place rsi, 2*64, 64
	jmp	exit

sscase54:
	sse2mac	128, 4096, x2cl_two_complex_unfft rsi, 2*64, 64
	jmp	exit
sscase55:
	sse2mac	128, 100000, x2cl_two_complex_unfft rsi, 2*64, 64
	jmp	exit

sscase56:
	sse2mac	128, 4096, x2cl_four_complex_fft rsi, 2*64, 64
	jmp	exit
sscase57:
	sse2mac	128, 100000, x2cl_four_complex_fft rsi, 2*64, 64
	jmp	exit
sscase58:
	sse2mac	128, 4096, x2cl_four_complex_first_fft rsi, 2*64, 64
	jmp	exit
sscase59:
	sse2mac	128, 100000, x2cl_four_complex_first_fft rsi, 2*64, 64
	jmp	exit
sscase60:
	sse2mac	128, 4096, s2cl_four_complex_gpm_fft rsi, 2*64, 64
	jmp	exit
sscase61:
	sse2mac	128, 100000, s2cl_four_complex_gpm_fft rsi, 2*64, 64
	jmp	exit
sscase62:
	sse2mac	128, 4096, s2cl_four_complex_first_fft rsi, 2*64, 64
	jmp	exit
sscase63:
	sse2mac	128, 100000, s2cl_four_complex_first_fft rsi, 2*64, 64
	jmp	exit
sscase64:
	sse2mac	128, 4096, s2cl_four_complex_fft_final rsi, 2*64, 64
	jmp	exit
sscase65:
	sse2mac	128, 100000, s2cl_four_complex_fft_final rsi, 2*64, 64
	jmp	exit
sscase66:
	sse2mac	128, 4096, s2cl_four_complex_with_square rsi, 2*64, 64
	jmp	exit
sscase67:
	sse2mac	128, 100000, s2cl_four_complex_with_square rsi, 2*64, 64
	jmp	exit
sscase68:
	sse2mac	128, 4096, s2cl_four_complex_with_mult rsi, 2*64, 64
	jmp	exit
sscase69:
	sse2mac	128, 100000, s2cl_four_complex_with_mult rsi, 2*64, 64
	jmp	exit
sscase70:
	sse2mac	128, 4096, s2cl_four_complex_with_mulf rsi, 2*64, 64
	jmp	exit
sscase71:
	sse2mac	128, 100000, s2cl_four_complex_with_mulf rsi, 2*64, 64
	jmp	exit
sscase72:
	sse2mac	256, 4096, x4cl_four_complex_fft rsi, 4*64, 64, 2*64
	jmp	exit
sscase73:
	sse2mac	256, 100000, x4cl_four_complex_fft rsi, 4*64, 64, 2*64
	jmp	exit
sscase74:
	sse2mac	256, 4096, x4cl_four_complex_cpm01_fft rsi, 4*64, 64, 2*64
	jmp	exit
sscase75:
	sse2mac	256, 100000, x4cl_four_complex_cpm01_fft rsi, 4*64, 64, 2*64
	jmp	exit
sscase76:
	sse2mac	256, 4096, x4cl_four_complex_cpm23_fft rsi, 4*64, 64, 2*64
	jmp	exit
sscase77:
	sse2mac	256, 100000, x4cl_four_complex_cpm23_fft rsi, 4*64, 64, 2*64
	jmp	exit
sscase78:
	sse2mac 256, 4096, g4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	jmp	exit
sscase79:
	sse2mac 256, 100000, g4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	jmp	exit
sscase80:
	sse2mac 256, 4096, x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	jmp	exit
sscase81:
	sse2mac 256, 100000, x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	jmp	exit
sscase82:
	sse2mac 256, 4096, x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	jmp	exit
sscase83:
	sse2mac 256, 100000, x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	jmp	exit
sscase84:
	sse2mac 256, 4096, x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64
	jmp	exit
sscase85:
	sse2mac 256, 100000, x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64
	jmp	exit

sscase86:
	sse2mac	128, 4096, x2cl_four_complex_unfft rsi, 2*64, 64
	jmp	exit
sscase87:
	sse2mac	128, 100000, x2cl_four_complex_unfft rsi, 2*64, 64
	jmp	exit
sscase88:
	sse2mac	128, 4096, g2cl_four_complex_unfft rsi, 2*64, 64, rdx, 2*64, 64
	jmp	exit
sscase89:
	sse2mac	128, 100000, g2cl_four_complex_unfft rsi, 2*64, 64, rdx, 2*64, 64
	jmp	exit
sscase90:
	sse2mac	256, 4096, x4cl_four_complex_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase91:
	sse2mac	256, 100000, x4cl_four_complex_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase92:
	sse2mac	256, 4096, x4cl_four_complex_last_unfft rsi, 4*64, 64, 2*64, 0
	jmp	exit
sscase93:
	sse2mac	256, 100000, x4cl_four_complex_last_unfft rsi, 4*64, 64, 2*64, 0
	jmp	exit
sscase94:
	sse2mac	256, 4096, s4cl_four_complex_gpm_unfft rsi, 4*64, 64, 2*64, 0
	jmp	exit
sscase95:
	sse2mac	256, 100000, s4cl_four_complex_gpm_unfft rsi, 4*64, 64, 2*64, 0
	jmp	exit
sscase96:
	sse2mac	256, 4096, x4cl_four_complex_cpm0_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase97:
	sse2mac	256, 100000, x4cl_four_complex_cpm0_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase98:
	sse2mac	256, 4096, x4cl_four_complex_cpm1_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase99:
	sse2mac	256, 100000, x4cl_four_complex_cpm1_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase100:
	sse2mac	256, 4096, x4cl_four_complex_cpm2_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase101:
	sse2mac	256, 100000, x4cl_four_complex_cpm2_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase102:
	sse2mac	256, 4096, x4cl_four_complex_cpm3_unfft rsi, 4*64, 64, 2*64
	jmp	exit
sscase103:
	sse2mac	256, 100000, x4cl_four_complex_cpm3_unfft rsi, 4*64, 64, 2*64
	jmp	exit

sscase104:
	sse2mac	192, 4096, x3cl_six_reals_first_fft rsi, 3*64, 64
	jmp	exit
sscase105:
	sse2mac	192, 100000, x3cl_six_reals_first_fft rsi, 3*64, 64
	jmp	exit
sscase106:
	sse2mac	192, 4096, x3cl_six_reals_last_unfft rsi, 3*64, 64
	jmp	exit
sscase107:
	sse2mac	192, 100000, x3cl_six_reals_last_unfft rsi, 3*64, 64
	jmp	exit

sscase108:
	sse2mac	192, 4096, x3cl_three_complex_first_fft rsi, 3*64, 64
	jmp	exit
sscase109:
	sse2mac	192, 100000, x3cl_three_complex_first_fft rsi, 3*64, 64
	jmp	exit
sscase110:
	sse2mac	192, 4096, x3cl_three_complex_last_unfft rsi, 3*64, 64
	jmp	exit
sscase111:
	sse2mac	192, 100000, x3cl_three_complex_last_unfft rsi, 3*64, 64
	jmp	exit

sscase112:
	sse2mac	320, 4096, x5cl_five_reals_first_fft rsi, 5*64, 64
	jmp	exit
sscase113:
	sse2mac	320, 100000, x5cl_five_reals_first_fft rsi, 5*64, 64
	jmp	exit
sscase114:
	sse2mac	320, 4096, x5cl_five_reals_last_unfft rsi, 5*64, 64
	jmp	exit
sscase115:
	sse2mac	320, 100000, x5cl_five_reals_last_unfft rsi, 5*64, 64
	jmp	exit

sscase116:
	sse2mac	448, 4096, x7cl_seven_reals_first_fft rsi, 7*64, 64
	jmp	exit
sscase117:
	sse2mac	448, 100000, x7cl_seven_reals_first_fft rsi, 7*64, 64
	jmp	exit
sscase118:
	sse2mac	448, 4096, x7cl_seven_reals_last_unfft rsi, 7*64, 64
	jmp	exit
sscase119:
	sse2mac	448, 100000, x7cl_seven_reals_last_unfft rsi, 7*64, 64
	jmp	exit




; Exit the timing code

exit:	pop	rbx
	pop	rbp
	pop	rdi
	pop	rsi
	ret

_timeit	ENDP
