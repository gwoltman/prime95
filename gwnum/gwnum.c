/*----------------------------------------------------------------------
| gwnum.c
|
| This file contains the C routines and global variables that are used
| in the multi-precision arithmetic routines.  That is, all routines
| that deal with the gwnum data type.
| 
|  Copyright 2002-2005 Just For Fun Software, Inc.
|  All Rights Reserved.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include "cpuid.h"
#include "gwnum.h"
#include "gwutil.h"
#include "gwdbldbl.h"

/* global variables */

double	GWSAFETY_MARGIN = 0.0;	/* Reduce maximum allowable bits per */
				/* FFT data word by this amount. */
long	GWMAXMULBYCONST = 3;	/* Gwsetup needs to know the maximum value */
				/* the caller will use in gwsetmulbyconst. */
				/* The default value is 3, commonly used */
				/* in a base-3 Fermat PRP test. */
double	KARG = 1.0;		/* K in K*B^N+C */
unsigned long BARG = 0;		/* B in K*B^N+C */
unsigned long PARG = 0;		/* N in K*B^N+C */
signed long CARG = -1;		/* C in K*B^N+C */
unsigned long FFTLEN = 0;	/* The FFT size we are using */
unsigned long RATIONAL_FFT = 0;	/* TRUE if bits per FFT word is integer */
unsigned long BITS_PER_WORD = 0;/* Bits in a little word */
int	GWERROR = 0;		/* Set if an error is detected */
double	MAXERR = 0.0;	/* Convolution error in a multiplication */
double	MAXDIFF = 0.0;		/* Maximum allowable difference between */
				/* sum of inputs and outputs */
unsigned long COPYZERO[8] = {0};/* Ptrs to help in gwcopyzero */
unsigned long ASM_TIMERS[32] = {0};/* Assembly language timers */
void	(*GWPROCPTRS[24])()={NULL}; /* Ptrs to assembly routines */
unsigned int NORMNUM = 0;	/* The post-multiply normalize routine index */
void (*NORMRTN)() = NULL;	/* The post-multiply normalization routine */
unsigned long POSTFFT = 0;	/* True if assembly code can start the */
				/* FFT process on the result of a multiply */
unsigned long ADDIN_ROW = 0;	/* For adding a constant after multiply */
unsigned long ADDIN_OFFSET = 0;
double ADDIN_VALUE = 0.0;

unsigned long EXTRA_BITS=0;	/* Number of unnormalized adds that can */
				/* be safely performed. */
int	SPREAD_CARRY_OVER_4_WORDS=0;/* True when carry out of top word */
				/* must be spread over more than 2 words */
int	TOP_CARRY_NEEDS_ADJUSTING=0;/* True when carry out of top word */
				/* needs adjusting */
double	INVERSE_KARG = 0.0;	/* 1/K */
double	KARG_HI = 0.0;		/* Upper bits of K */
double	KARG_LO = 0.0;		/* Lower bits of K */
double	CARRY_ADJUST1 = 0.0;	/* Adjustment constant #1 in wrapping carry */
double	CARRY_ADJUST2 = 0.0;	/* Adjustment constant #2 in wrapping carry */
double	CARRY_ADJUST3 = 0.0;	/* Adjustment constant #3 in wrapping carry */
double	CARRY_ADJUST4 = 0.0;	/* Adjustment constant #4 in wrapping carry */
double	CARRY_ADJUST5 = 0.0;	/* Adjustment constant #5 in wrapping carry */
double	CARRY_ADJUST6 = 0.0;	/* Adjustment constant #6 in wrapping carry */
unsigned long HIGH_WORD1_OFFSET=0;/* Offset of top FFT word */
unsigned long HIGH_WORD2_OFFSET=0;/* Offset of second high FFT word */
unsigned long HIGH_WORD3_OFFSET=0;/* Offset of third high FFT word */
unsigned long HIGH_SCRATCH1_OFFSET=0;
				/* Offset of top FFT word in scratch area */
unsigned long HIGH_SCRATCH2_OFFSET=0;
				/* Offset of second highest FFT word */
unsigned long HIGH_SCRATCH3_OFFSET=0;
				/* Offset of third highest FFT word */

unsigned long ZPAD_TYPE = 0;	/* 1,2,or 3 words in k (used by zero pad) */
double	ZPAD_INVERSE_K6 = 0.0;	/* Zero padded FFT constants */
double	ZPAD_INVERSE_K5 = 0.0;
double	ZPAD_INVERSE_K4 = 0.0;
double	ZPAD_INVERSE_K3 = 0.0;
double	ZPAD_INVERSE_K2 = 0.0;
double	ZPAD_INVERSE_K1 = 0.0;
double	ZPAD_K6_HI = 0.0;
double	ZPAD_K5_HI = 0.0;
double	ZPAD_K4_HI = 0.0;
double	ZPAD_K3_HI = 0.0;
double	ZPAD_K2_HI = 0.0;
double	ZPAD_K1_HI = 0.0;
double	ZPAD_K6_MID = 0.0;
double	ZPAD_K5_MID = 0.0;
double	ZPAD_K4_MID = 0.0;
double	ZPAD_K3_MID = 0.0;
double	ZPAD_K2_MID = 0.0;
double	ZPAD_K1_MID = 0.0;
double	ZPAD_K6_LO = 0.0;
double	ZPAD_K5_LO = 0.0;
double	ZPAD_K4_LO = 0.0;
double	ZPAD_K3_LO = 0.0;
double	ZPAD_K2_LO = 0.0;
double	ZPAD_K1_LO = 0.0;
double	ZPAD_SHIFT6 = 0.0;
double	ZPAD_SHIFT5 = 0.0;
double	ZPAD_SHIFT4 = 0.0;
double	ZPAD_SHIFT3 = 0.0;
double	ZPAD_SHIFT2 = 0.0;
double	ZPAD_SHIFT1 = 0.0;

unsigned long BIGLIT_INCR2 = 0;	/* Offset to step in big/lit array */
unsigned long BIGLIT_INCR4 = 0;	/* Offset to step in big/lit array */

unsigned long* INFT[4] = {0};	/* For assembly language arg passing */
unsigned long GWVERNUM = 0;	/* For assembly language arg passing */
void	*SRCARG = NULL;		/* For assembly language arg passing */
void	*SRC2ARG = NULL;	/* For assembly language arg passing */
void	*DESTARG = NULL;	/* For assembly language arg passing */
void	*DEST2ARG = NULL;	/* For assembly language arg passing */
long	NUMARG = 0;		/* For assembly language arg passing */
long	NUM2ARG = 0;		/* For assembly language arg passing */
double	fft_count = 0;		/* Count of forward and inverse FFTs */
double	*gwnum_memory;		/* Allocated memory */
unsigned long GW_ALIGNMENT = 0;	/* How to align allocated gwnums */
unsigned long GW_ALIGNMENT_MOD = 0; /* How to align allocated gwnums */
unsigned long PASS1_CACHE_LINES = 0; /* Cache lines grouped together in */
				/* first pass of an FFT. */
unsigned long PASS2_LEVELS = 0; /* FFT levels done in pass 2. */
unsigned long SCRATCH_SIZE = 0;	/* Size of the pass 1 scratch area */
double	bit_length;		/* Bit length of k*b^n */
int	ZERO_PADDED_FFT=0;	/* True if doing a zero pad FFT */
double	fft_bits_per_word;	/* Num bits in each fft word */
double	fft_max_bits_per_word;	/* Maximum bits per data word that */
				/* this FFT size can support */

int	GENERAL_MOD = 0;	/* True if doing general-purpose mod */
				/* as defined in gwsetup_general_mod. */
giant	GW_MODULUS = NULL;	/* In the general purpose mod case, this is */
				/* the number operations are modulo. */
gwnum	GW_MODULUS_FFT = NULL;	/* In the general purpose mod case, this is */
				/* the FFT of GW_MODULUS. */
gwnum	GW_RECIP = NULL;	/* Shifted reciprocal of of GW_MODULUS */
unsigned long GW_ZEROWORDSLOW=0;/* Count of words to zero during copy step */
				/* of a general purpose mod. */

gwnum	GW_RANDOM = NULL;	/* A random number used in */
				/* gwsquare_carefully. */

char	GWSTRING_REP[40];	/* The gwsetup modulo number as a string. */

char	*GW_BIGBUF = NULL;	/* Optional buffer to allocate gwnums in */
unsigned long GW_BIGBUF_SIZE = 0;/* Size of the optional buffer */

gwnum	*gwnum_alloc = NULL;	/* Array of allocated gwnums */
unsigned int gwnum_alloc_count = 0; /* Count of allocated gwnums */
unsigned int gwnum_alloc_array_size = 0; /* Size of gwnum_alloc array */
gwnum	*gwnum_free = NULL;	/* Array of available gwnums */
unsigned int gwnum_free_count = 0; /* Count of available gwnums */

/* When debugging gwnum and giants, I sometimes write code that "cheats" */
/* by calling a routine that is part of prime95 rather than the gwnum and */
/* giants library.  Prime95 will set this routine pointer so that gwnum */
/* code can cheat while keeping the gwnum library interface clean. */

void (*OutputBothRoutine)(char *) = NULL;
#define OutputBoth(x)	(*OutputBothRoutine)(x)

/* Assembly helper routines */

void gwinfo1 (void);
void gwsetup2 (void);

/* gwnum assembly routine pointers */

#define gw_fft()	(*GWPROCPTRS[0])()
#define gw_square()	(*GWPROCPTRS[1])()
#define gw_mul()	(*GWPROCPTRS[2])()
#define gw_mulf()	(*GWPROCPTRS[3])()
#define gw_add()	(*GWPROCPTRS[4])()
#define gw_addq()	(*GWPROCPTRS[5])()
#define gw_sub()	(*GWPROCPTRS[6])()
#define gw_subq()	(*GWPROCPTRS[7])()
#define gw_addsub()	(*GWPROCPTRS[8])()
#define gw_addsubq()	(*GWPROCPTRS[9])()
#define gw_copyzero()	(*GWPROCPTRS[10])()
#define gw_addf()	(*GWPROCPTRS[11])()
#define gw_subf()	(*GWPROCPTRS[12])()
#define gw_addsubf()	(*GWPROCPTRS[13])()
#define norm_routines	14

/* Helper macros */

#define fftinc(x)		(fft_count += x)
#define isinf(x)		((x) != 0.0 && (x) == 2.0*(x))
#define isnan(x)		((x) != (x))
#define is_valid_double(x)	(! isnan (x) && ! isinf (x))

/* Error codes returned */

#define GWERROR_BAD_FFT_DATA	-1

/* Forward declarations */

int internal_gwsetup (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. Must be two. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c,		/* C in K*B^N+C. Must be rel. prime to K. */
	unsigned long fftlen);	/* Specific FFT size to use (or zero) */
double virtual_bits_per_word ();
void raw_gwsetaddin (unsigned long word, long val);

/* Find the power of two greater than or equal to N. */

unsigned long pow_two_above_or_equal (
	unsigned long n)
{
static	unsigned long save_n = 0;
static	unsigned long save_result = 0;

	if (n != save_n) {
		save_n = n;
		save_result = 1;
		for (n = n - 1; n; n = n >> 1) save_result = save_result << 1;
	}
	return (save_result);
}


/* This routine builds a sin/cos table - used by gwsetup */

double *build_sin_cos_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long N,	/* Number of DATA values processed by this */
				/* FFT level.  This explains the divide by 2 */
				/* for complex FFTs later in this routine */
	int	hermetian_skip,	/* True if some sin/cos values are skipped */
	int	type)		/* 0 = old style - a plain old array */
				/* 1 = SSE2 - data is duplicated */
				/* 2 = SSE2 - data is interleaved */
{
	unsigned long i;

/* Handle hermetian skip when interleaving.  First data slot is left */
/* undefined. */

	if (type == 2 && hermetian_skip) type = 3;

/* Special case the really small sin/cos tables.  If N is between 9 and 16 */
/* or between 33 and 64, then the assembly code is only doing one FFT level. */
/* In this case, the code just uses the middle sin/cos values of a 2N sized */
/* table.  We could optimize this inefficient memory usage at a later date. */

	if (N <= 8) return (table);
	if (N >= 9 && N <= 16) N = N * 2;
	if (N >= 33 && N <= 64 && type == 1 && hermetian_skip) N = N * 2;

/* In the all-complex case. build the same size table as the hermetian */
/* case which skips half the i values. */

	if (!hermetian_skip) N = N / 2;

/* Loop to build table. */

	for (i = hermetian_skip ? ((N & 4) ? 4 : 8) : 0; i < N; i += 4) {
		unsigned long shifted_i, shifted_N, flipped_i;
		double	sincos[6];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (hermetian_skip) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);
				if (j == 3) continue;
			}
		}

/* Compute the 3 sin/cos values */

		NUMARG = (long) flipped_i;
		NUM2ARG = (long) N;
		DESTARG = (void *) &sincos;
		gwsincos3 (NUMARG, NUM2ARG, (double *) &sincos);

/* Copy the sin/cos values in the appropriate way */

		if (type == 0) {
			memcpy (table, sincos, sizeof (sincos));
			table += 6;
		} else if (type == 1) {
			table[0] = table[1] = sincos[0];
			table[2] = table[3] = sincos[1];
			table[4] = table[5] = sincos[2];
			table[6] = table[7] = sincos[3];
			table[8] = table[9] = sincos[4];
			table[10] = table[11] = sincos[5];
			table += 12;
		} else if (type == 2) {
			table[0] = sincos[0];
			table[2] = sincos[1];
			table[4] = sincos[2];
			table[6] = sincos[3];
			table[8] = sincos[4];
			table[10] = sincos[5];
			type++;
		} else {
			table[1] = sincos[0];
			table[3] = sincos[1];
			table[5] = sincos[2];
			table[7] = sincos[3];
			table[9] = sincos[4];
			table[11] = sincos[5];
			type--;
			table += 12;
		}
	}
	return (table);
}

/* This routine builds a pass 2 premultiplier table - used by gwsetup */

double *build_premult_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned long i, N, incr, type;

/* Build a premultiplier table for the second pass incrementing by */
/* the pre-calculated pass2_size. */

	N = FFTLEN;
	incr = pass2_size;
	if (CARG > 0 && !ZERO_PADDED_FFT) N = N / 2;

/* Mod 2^N+1 arithmetic starts at first data set, */
/* mod 2^N-1 skips some data sets */

 	if (CARG > 0 && !ZERO_PADDED_FFT) i = 0;
	else i = incr * 4;

/* Loop to build table. */

	type = 0;
	for ( ; i < N; i += incr) {
		unsigned long shifted_i, shifted_N, flipped_i, k, l;
		double	sincos[2];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (CARG < 0 || ZERO_PADDED_FFT) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);	
				if (j == 3) continue;
			}
		}

/* Generate the group multipliers */

		for (k = 0; k < incr / 4; k += 4) {

/* There are 4 multipliers in a XMM_PMD set */

			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value (root of unity) */

				if (CARG < 0 || ZERO_PADDED_FFT) {
					NUMARG = (long) (((l * incr/4 + k) * flipped_i) % N);
					NUM2ARG = (long) N;
					DESTARG = (void *) &sincos;
					gwsincos (NUMARG, NUM2ARG, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					NUMARG = (long) (((l * incr/4 + k) * flipped_i * 4 + l*incr/4+k) % (N*4));
					NUM2ARG = (long) (N*4);
					DESTARG = (void *) &sincos;
					gwsincos (NUMARG, NUM2ARG, (double *) &sincos);
				}

/* Save the premultiplier value */

				table[l*4+type] = sincos[0];
				table[l*4+2+type] = sincos[1];
			}
			table += 16;
		}
	
/* Generate the 4 column multipliers * 4 sin/cos values */

		for (k = 0; k < 4; k++) {
			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value (root of unity) */

				if (CARG < 0 || ZERO_PADDED_FFT) {
					NUMARG = (long) ((k * flipped_i + l * N/16) % N);
					NUM2ARG = (long) N;
					DESTARG = (void *) &sincos;
					gwsincos (NUMARG, NUM2ARG, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					NUMARG = (long) (((k * flipped_i * 2 + l * N/8) *2 + k) % (N*4));
					NUM2ARG = (long) (N*4);
					DESTARG = (void *) &sincos;
					gwsincos (NUMARG, NUM2ARG, (double *) &sincos);
				}

/* Save the premultiplier value */

				table[l*4+type] = sincos[0];
				table[l*4+2+type] = sincos[1];
			}
			table += 16;
		}

		if (type == 0) table -= (incr / 4 + 16) * 4;
		type = 1 - type;
 	}

	return (table);
}

/* This routine builds a plus 1 premultiplier table - used by gwsetup */
/* when c is positive. */

double *build_plus1_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned long i, j, k, l, N;
	int	pfa;

/* Set flag if this is a 3*2^n FFT */

	pfa = (FFTLEN != pow_two_above_or_equal (FFTLEN));

/* Adjust for two-pass FFTs */

	if (pass2_size == 1) N = FFTLEN;
	else N = FFTLEN / (pass2_size / 2);

/* Loop to build premultiplier table in the same order as the underlying */
/* assembly macro needs them.  The pfa macro operates on 3 cache lines */
/* while the power-of-two macro operates on 2 cache lines. */
/* A 64 length FFT needs 0,8,16,24 for the macro then 3 more iterations */
/* for the cache lines beginning with 2,4,6. */
/* A 48 length FFT needs 0,8,16 and 4,12,20 for the first macro then */
/* one more iteration for the cache lines beginning with 2. */

	for (i = 0; i < N / (pfa ? 24 : 32); i++) {
	for (l = 0; l < 2; l++) {
		double	sincos[2];

/* Generate the pre multipliers (roots of -1). */

		for (k = 0; k < (unsigned long) (pfa ? 3 : 4); k++) {
		for (j = 0; j < 2; j++) {

/* Compute the sin/cos value */

			if (pfa)
				NUMARG = (long) ((i * 2 + l * N/12 + j + k * N/6) % N);
			else
				NUMARG = (long) ((i * 4 + l * 2 + j + k * N/8) % N);
			NUM2ARG = (long) (N*2);
			DESTARG = (void *) &sincos;
			gwsincos (NUMARG, NUM2ARG, (double *) &sincos);

/* Save the premultiplier value */

			table[0+j] = sincos[0];
			table[2+j] = sincos[1];

/* For two-pass FFTs we could apply the root of -1 for the upper SSE2 */
/* double here or in the pass 2 premultipliers.  We've arbitrarily chosen */
/* to do it in the pass 2 premults. */

			if (pass2_size > 1) {
				j = 1;
				table[0+j] = sincos[0];
				table[2+j] = sincos[1];
			}
		}
		table += 4;
		}
	}
 	}

	return (table);
}

/* This routine builds a normalization table - used by SSE2 normalizaion */
/* routines */

double *build_norm_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size, /* Size of pass2 */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (pass2_size == 1) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < FFTLEN; i++) {
			unsigned long j, table_entry;
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (i, &ttp, NULL, &ttmp);

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (FFTLEN, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the entry for the MSW or LSW in an SSE2 pair */

			table[table_entry*4+(j&1)] = ttmp;
			table[table_entry*4+2+(j&1)] = ttp;
		}
		return (table + FFTLEN + FFTLEN);
	}

/* Two pass FFTs are handled here */

	num_cols = pass2_size / 2;
	if (col) {

/* Loop to build table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (i, &ttp, NULL, &ttmp);

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

			table[i*4] = ttmp;
			table[i*4+1] = ttmp;
			table[i*4+2] = ttp;
			table[i*4+3] = ttp;
		}
		return (table + num_cols * 4);
	}

/* Build the group multipliers table */

	else {
		unsigned long pfa, h, hlimit, haddin, m, mmult, u, umult;

/* Determine if this is a PFA 5, 6, 7, or 8 */

		for (pfa = FFTLEN; pfa > 8; pfa >>= 1);

/* Loop to build table */

		umult = FFTLEN / 2;
		hlimit = FFTLEN / 4 / (2*num_cols);
		for (h = 0; h < hlimit; h++) {
			if (pfa == 5) {
				if (h < hlimit / 5) {
					haddin = h * 2 * num_cols;
					mmult = FFTLEN / 20;
				} else {
					haddin = FFTLEN/10 + (h - hlimit/5) * 2 * num_cols;
					mmult = FFTLEN / 5;
				}
			} else if (pfa == 7) {
				if (h < hlimit / 7) {
					haddin = h * 2 * num_cols;
					mmult = FFTLEN / 28;
				} else if (h < 3 * hlimit / 7) {
					haddin = FFTLEN/14 + (h - hlimit/7) * 2 * num_cols;
					mmult = FFTLEN / 14;
				} else {
					haddin = 3*FFTLEN/14 + (h - 3*hlimit/7) * 2 * num_cols;
					mmult = FFTLEN / 7;
				}
			} else {
				haddin = h * 2 * num_cols;
				mmult = FFTLEN / 4;
			}
			for (u = 0; u < 2; u++) {
			for (m = 0; m < 2; m++) {
			for (k = 0; k < 2; k++) {
				double	ttp, ttmp;
				long	n;

/* Call asm routines to compute the two multipliers */

				n = haddin + u * umult + m * mmult + k * num_cols;
				gwfft_weights3 (n, &ttp, &ttmp, NULL);

/* Now set the entry for BOTH the MSW and LSW in an SSE2 pair */

				table[k] = ttmp;
				table[2+k] = ttp;
			}
			table += 4;
			}
			}
		}
		return (table);
	}
}

/* This routine builds a big/little flags table - used by SSE2 normalizaion */
/* routines */

double *build_biglit_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned char *p;
	unsigned long h, i, j, k, m, u, gap;
	unsigned long pfa, hlimit, haddin, mmult, umult;

/* Handle one pass FFTs differently */

	if (pass2_size == 1) {

/* Loop to build table */

		p = (unsigned char *) table;
		for (i = 0; i < FFTLEN; i++) {
			unsigned long table_entry;

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (FFTLEN, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the biglit table entry for a LSW in an SSE2 pair */

			if ((j & 1) == 0) {
				p[table_entry] = is_big_word (i) * 16;
			}

/* Otherwise, set the biglit table entry for a MSW in an SSE2 pair */

			else {
				if (is_big_word (i)) p[table_entry] += 32;
			}
		}
		return ((double *) (p + FFTLEN / 2));
	}

/* Determine if this is a PFA 5, 6, 7, or 8 */

	for (pfa = FFTLEN; pfa > 8; pfa >>= 1);

/* Determine the gap between XMM high and low words */

	gap = pass2_size / 2;

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code.  This is especially ugly in the PFA cases */

	p = (unsigned char *) table;
	umult = FFTLEN / 2;
	hlimit = FFTLEN / 4 / (2*gap);
	for (i = 0; i < gap; i += PASS1_CACHE_LINES) {
	for (h = 0; h < hlimit; h++) {
		if (pfa == 5) {
			if (h < hlimit / 5) {
				haddin = h * 2 * gap;
				mmult = FFTLEN / 20;
			} else {
				haddin = FFTLEN/10 + (h - hlimit/5) * 2 * gap;
				mmult = FFTLEN / 5;
			}
		} else if (pfa == 7) {
			if (h < hlimit / 7) {
				haddin = h * 2 * gap;
				mmult = FFTLEN / 28;
			} else if (h < 3 * hlimit / 7) {
				haddin = FFTLEN/14 + (h - hlimit/7) * 2 * gap;
				mmult = FFTLEN / 14;
			} else {
				haddin = 3*FFTLEN/14 + (h - 3*hlimit/7) * 2 * gap;
				mmult = FFTLEN / 7;
			}
		} else {
			haddin = h * 2 * gap;
			mmult = FFTLEN / 4;
		}
	for (j = 0; j < PASS1_CACHE_LINES; j++) {
	for (u = 0; u < 2; u++) {
	for (m = 0; m < 2; m++) {
	for (k = 0; k < 2 * gap; k += gap) {
		unsigned long word;

/* Now set the big/little flag for a LSW in an SSE2 pair */
/* Otherwise, set the big/little flag for a MSW in an SSE2 pair */

		word = haddin + i + j + u * umult + m * mmult + k;
		if (k == 0) *p = is_big_word (word) * 16;
		else if (is_big_word (word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs.  The fudge flag is */
/* set if the col mult * the grp mult is twice the correct fft_weight, */
/* meaning a mul by 0.5 is required to generate the correct multiplier. */
/* Since we can't do equality compares on floats, this test is a little bit */
/* cryptic. */

		if (gwfft_weight_exponent (word&(gap-1)) +
		    gwfft_weight_exponent (word&~(gap-1)) > 1.0) {
			if (k == 0) *p += 64;
			else *p += 128;
		}

/* Set some offsets that help the assembly code step through the big/lit */
/* array in a non-traditional order.  Two pass-FFTs step through the array */
/* in chunks of PASS1_CACHE_LINES, but the add, sub, and carry propagation */
/* code need to access the big/lit array linearly.  Set two global variables */
/* that tell the assembly code the big/lit array distance between words */
/* 0 and 2, and words 0 and 4. */

		if (word == 2) 
			BIGLIT_INCR2 = (char *) p - (char *) table;
		if (word == 4) 
			BIGLIT_INCR4 = (char *) p - (char *) table;
	}
	p++;
	}
	}
	}
	}
	}
	return ((double *) p);
}


/* This routine builds an x87 sin/cos table - used by gwsetup */

double *build_x87_sin_cos_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long N,
	int	hermetian_skip)	/* True if some sin/cos values are skipped */
{
	unsigned long i;

/* Special case the really small sin/cos tables.  If N is between 9 and 16 */
/* then the assembly code is only doing one FFT level. */
/* In this case, the code just uses the middle sin/cos values of a 2N sized */
/* table.  We could optimize this inefficient memory usage at a later date. */

	if (N <= 8) return (table);
	if (N >= 9 && N <= 16) N = N * 2;

/* The N value passed in represents the number of real numbers that are */
/* processed in a section.  If heremetian_skip is not set, then we are */
/* instead dealing with complex numbers and there are half as many complex */
/* numbers in a section.  For example, when doing 8 levels in pass 2, this */
/* routine is called with N=512.  The first real section has 512 values, */
/* while the remaining pass 2 sections have 256 complex values. */

	if (!hermetian_skip) N = N / 2;

/* Loop to build table */

	for (i = hermetian_skip ? ((N & 4) ? 4 : 8) : 0; i < N; i += 4) {
		unsigned long shifted_i, shifted_N, flipped_i;
		double	sincos[6];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (hermetian_skip) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);
				if (j == 3) continue;
			}
		}

/* Compute the 3 sin/cos values */

		NUMARG = (long) flipped_i;
		NUM2ARG = (long) N;
		DESTARG = (void *) &sincos;
		gwsincos3 (NUMARG, NUM2ARG, (double *) &sincos);

/* Copy the sin/cos values to the table */

		memcpy (table, sincos, sizeof (sincos));
		table += 6;
	}
	return (table);
}

/* This routine builds a pass 2 premultiplier table - used by gwsetup */

double *build_x87_premult_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned long i, N, incr;

/* Build a premultiplier table for the second pass incrementing by */
/* the pre-calculated pass2_size. */

	N = FFTLEN;
	incr = pass2_size;
	if (CARG > 0 && !ZERO_PADDED_FFT) N = N / 2;

/* Mod 2^N+1 arithmetic starts at first data set, */
/* mod 2^N-1 skips some data sets */

	if (CARG > 0 && !ZERO_PADDED_FFT) i = 0;
	else i = incr * 2;

/* Loop to build table */

	for ( ; i < N; i += incr) {
		unsigned long shifted_i, shifted_N, flipped_i, k, l;
		double	sincos[2];

/* Flip the bits in i.  Our prime-factor-FFT makes this a little complex. */
/* The algorithm below works, but I've long since forgotten why. */

		shifted_i = i; shifted_N = N; flipped_i = 0;
		while ((shifted_N & 1) == 0) {
			flipped_i <<= 1;
			if (shifted_i & 1) flipped_i++;
			shifted_i >>= 1;
			shifted_N >>= 1;
		}
		flipped_i = (flipped_i * shifted_N) + shifted_i;

/* When the FFT is working on real data Hermetian symettry allows us to */
/* eliminate half of the FFT data and consequently half of the sin/cos data */
/* Case 1:  If shifted source is > shifted N/2, then we */
/* do not need these sin/cos values. */
/* Case 2:  If shifted source is zero, loop to find the top */
/* two bits.  Skip the number if the top two bits equal 3. */

		if (CARG < 0 || ZERO_PADDED_FFT) {
			if (shifted_i > shifted_N / 2) continue;
			if (shifted_i == 0) {
				unsigned long j;
				for (j = i; j > 3; j >>= 1);	
				if (j == 3) continue;
			}
		}

/* Generate the group multipliers */

		for (k = 0; k < incr / 4; k += 4) {

/* There are 4 multipliers in a PMD set */

			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value (root of unity) */

				if (CARG < 0 || ZERO_PADDED_FFT) {
					NUMARG = (long) (((l * incr/4 + k) * flipped_i) % N);
					NUM2ARG = (long) N;
					DESTARG = (void *) &sincos;
					gwsincos (NUMARG, NUM2ARG, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					NUMARG = (long) (((l * incr/4 + k) * flipped_i * 4 + l*incr/4+k) % (N*4));
					NUM2ARG = (long) (N*4);
					DESTARG = (void *) &sincos;
					gwsincos (NUMARG, NUM2ARG, (double *) &sincos);
				}

/* Save the premultiplier values */

				table[l*2] = sincos[0];
				table[l*2+1] = sincos[1];
			}
			table += 8;
		}
	
/* Generate the 4 column multipliers * 4 sin/cos values */

		for (k = 0; k < 4; k++) {
			for (l = 0; l < 4; l++) {

/* Compute the sin/cos value (root of unity) */

				if (CARG < 0 || ZERO_PADDED_FFT) {
					NUMARG = (long) ((k * flipped_i + l * N/16) % N);
					NUM2ARG = (long) N;
					DESTARG = (void *) &sincos;
					gwsincos (NUMARG, NUM2ARG, (double *) &sincos);
				}

/* If C > 0, then also multiply by the proper root of -1.  This is done */
/* by changing the value we are taking the sin/cos of */

				else {
					NUMARG = (long) (((k * flipped_i * 2 + l * N/8) *2 + k) % (N*4));
					NUM2ARG = (long) (N*4);
					DESTARG = (void *) &sincos;
					gwsincos (NUMARG, NUM2ARG, (double *) &sincos);
				}

/* Save the premultiplier value */

				table[l*2] = sincos[0];
				table[l*2+1] = sincos[1];
			}
			table += 8;
		}
 	}

	return (table);
}

/* This routine builds a plus 1 premultiplier table - used by gwsetup */
/* when c is positive. */

double *build_x87_plus1_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned long i, k, N;
	int	pfa;

/* Set flag if this is a 3*2^n FFT */

	pfa = (FFTLEN != pow_two_above_or_equal (FFTLEN));

/* Adjust for two-pass FFTs */

	if (pass2_size == 1) N = FFTLEN;
	else N = FFTLEN / pass2_size;

/* Loop to build premultiplier table in the same order as the underlying */
/* assembly macro needs them. */

	for (i = 0; i < N / (pfa ? 6 : 8); i++) {
		double	sincos[2];

/* Generate the pre multipliers (roots of -1) used in one three_complex */
/* or four complex macro. */

		for (k = 0; k < (unsigned long) (pfa ? 3 : 4); k++) {

/* Compute the sin/cos value */

			if (pfa)
				NUMARG = (long) ((i + k * N/6) % N);
			else
				NUMARG = (long) ((i + k * N/8) % N);
			NUM2ARG = (long) (N*2);
			DESTARG = (void *) &sincos;
			gwsincos (NUMARG, NUM2ARG, (double *) &sincos);

/* Save the premultiplier value */

			table[0] = sincos[0];
			table[1] = sincos[1];
			table += 2;
		}
	}

	return (table);
}

/* This routine builds a normalization table - used by x87 normalizaion */
/* routines */

double *build_x87_norm_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size, /* Size of pass2 */
	int	col)		/* TRUE if building column, not group, table */
{
	unsigned long i, k, num_cols;

/* Handle one-pass FFTs first, there are no group multipliers */

	if (pass2_size == 1) {
		if (!col) return (table);

/* Loop to build table */

		for (i = 0; i < FFTLEN; i++) {
			unsigned long j;
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (i, &ttp, NULL, &ttmp);

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (FFTLEN, i) / sizeof (double);

/* Now set the appropriate table entry.  These are put into the array */
/* in the same order that the normalization code needs them. */

			table[j*2] = ttmp;
			table[j*2+1] = ttp;
		}
		return (table + FFTLEN + FFTLEN);
	}

/* Two pass FFTs are handled here */

	num_cols = pass2_size;
	if (col) {

/* Loop to build columns table */

		for (i = 0; i < num_cols; i++) {
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (i, &ttp, NULL, &ttmp);

/* Now set the appropriate table entry.  These are put into the array */
/* in the same order that the normalization code needs them. */

			table[i+i] = ttmp;
			table[i+i+1] = ttp;
		}
		return (table + num_cols * 2);
	}

/* Build the group multipliers table */

	else {
		unsigned long num_grps;
		
/* Loop to build group table */

		num_grps = FFTLEN / num_cols;
		for (i = 0; i < num_grps; i++) {
			double	ttp, ttmp;

/* Call asm routines to compute the two multipliers */

			gwfft_weights3 (i * num_cols, &ttp, &ttmp, NULL);

/* Now set the appropriate table entry.  These are put into the array */
/* in the same order that the normalization code needs them. */

			if (i < num_grps / 2) k = i * 2;
			else k = (i - num_grps / 2) * 2 + 1;
			table[k+k] = ttmp;
			table[k+k+1] = ttp;
		}
		return (table + num_grps * 2);
	}
}

/* This routine builds a big/little flags table - used by x87 normalizaion */
/* routines */

double *build_x87_biglit_table (
	double	*table,		/* Pointer to the table to fill in */
	unsigned long pass2_size)
{
	unsigned char *p;
	unsigned long i, j, k, m;

/* Handle one pass FFTs differently */

	if (pass2_size == 1) {

/* Loop to build table */

		p = (unsigned char *) table;
		for (i = 0; i < FFTLEN; i++) {
			unsigned long table_entry;

/* Find where this data appears in the FFT array and in the table we are building. */

			j = addr_offset (FFTLEN, i) / sizeof (double);
			table_entry = j >> 1;

/* Now set the biglit table entry for a LSW in a pair */

			if ((j & 1) == 0) {
				p[table_entry] = is_big_word (i) * 16;
			}

/* Otherwise, set the biglit table entry for a MSW in a pair */

			else {
				if (is_big_word (i)) p[table_entry] += 32;
			}
		}
		return ((double *) (p + FFTLEN / 2));
	}

/* Loop to build table in exactly the same order that it will be */
/* used by the assembly code. */

	p = (unsigned char *) table;
	for (i = 0; i < pass2_size; i += PASS1_CACHE_LINES * 2) {
	for (j = 0; j < FFTLEN / 2; j += pass2_size) {
	for (k = 0; k < PASS1_CACHE_LINES * 2; k++) {
	for (m = 0; m < FFTLEN; m += FFTLEN / 2) {
		unsigned long word;

/* Now set the big/little flag for a LSW in a pair */
/* Otherwise, set the big/little flag for a MSW in a pair */

		word = i + j + k + m;
		if (m == 0) *p = is_big_word (word) * 16;
		else if (is_big_word (word)) *p += 32;

/* Set the ttp and ttmp fudge flags for two pass FFTs */
/* The fudge flag is set if col mult * grp mult will be greater than 2 */

		if (gwfft_weight_exponent (word & (pass2_size-1)) +
		    gwfft_weight_exponent (word & ~(pass2_size-1)) > 1.0) {
			if (m == 0) *p += 64;
			else *p += 128;
		}

/* Set some offsets that help the assembly code step through the big/lit */
/* array in a non-traditional order.  Two pass-FFTs step through the array */
/* in chunks of PASS1_CACHE_LINES, but the add, sub, and carry propagation */
/* code need to access the big/lit array linearly.  Set two global variables */
/* that tell the assembly code the big/lit array distance between words */
/* 0 and 2, and words 0 and 4. */

		if (word == 2) 
			BIGLIT_INCR2 = (char *) p - (char *) table;
		if (word == 4) 
			BIGLIT_INCR4 = (char *) p - (char *) table;
	}
	p++;
	}
	}
	}
	return ((double *) p);
}


/* This routine used to be in assembly language.  It scans the assembly */
/* code arrays looking for the best FFT size to implement our k*b^n+c FFT. */
/* Returns 0 for IBDWT FFTs, 1 for zero padded FFTs, or a gwsetup error */
/* code. */

int gwinfo (			/* Return zero-padded fft flag or error code */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* N in K*B^N+C. Base must be two. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c,		/* C in K*B^N+C. Must be rel. prime to K. */
	unsigned long fftlen)	/* Specific FFT size to use (or zero) */
{
	unsigned long *jmptab, *zpad_jmptab;
	double	log2k, log2c, log2maxmulbyconst;
	double	max_bits_per_word, bits_per_word;
	unsigned long l2_cache_size, max_exp;
	char	buf[20];

/* If L2 cache size is unknown, assume it is 128KB */

	if (CPU_L2_CACHE_SIZE >= 0)
		l2_cache_size = CPU_L2_CACHE_SIZE;
	else
		l2_cache_size = 128;

/* Get pointer to 4 assembly jmptables and the version number */

	gwinfo1 ();

/* Make sure that the assembly code version number matches the C version */
/* number.  If they do not match, then the user linked in the wrong gwnum */
/* object files! */

	sprintf (buf, "%d.%d", GWVERNUM / 100, GWVERNUM % 100);
	if (strcmp (buf, GWNUM_VERSION)) return (GWERROR_VERSION);

/* Precalculate some needed values */

	log2k = log (k) / log ((double) 2.0);
	log2c = log ((double) abs (c)) / log ((double) 2.0);
	log2maxmulbyconst = log ((double) abs (GWMAXMULBYCONST)) / log ((double) 2.0);

/* First, see what FFT length we would get if we emulate the k*b^n+c modulo */
/* with a zero padded FFT.  If k is 1 and abs (c) is 1 then we can skip this */
/* loop as we're sure to find an IBDWT that will do the job. */

	zpad_jmptab = NULL;
	if (fftlen == 0 && (k > 1.0 || n < 500 || abs (c) > 1)) {

/* Use the proper 2^N-1 jmptable */

		if (CPU_FLAGS & CPU_SSE2) zpad_jmptab = INFT[0];
		else zpad_jmptab = INFT[2];

/* Find the table entry for the FFT that can do a mod 2^2n FFT, handling */
/* k and c in the normalization routines.  We will compare this to the */
/* non-zero-padded FFT length later.  The zeroes in the upper half of FFT */
/* input data let us get about another 0.3 bits per input word. */

		while ((max_exp = zpad_jmptab[0]) != 0) {

/* Check L2 cache size constraints */
		
			if (l2_cache_size < ((zpad_jmptab[8] >> 16) & 0x7FFF))
				goto next1;

/* Check FFT requires prefetch capability */

			if (zpad_jmptab[8] & 0x80000000 &&
			    ! (CPU_FLAGS & CPU_PREFETCH))
				goto next1;

/* Compare the maximum number of bits allowed in the FFT input word */
/* with the number of bits we would use.  Break when we find an acceptable */
/* FFT length. */

			max_bits_per_word = (double) max_exp / zpad_jmptab[1];
			max_bits_per_word -= GWSAFETY_MARGIN;
			bits_per_word = (double) (n + n) / zpad_jmptab[1];
			if (bits_per_word < max_bits_per_word + 0.3) break;

/* Move to next jmptable entry */

next1:			for (zpad_jmptab += 14; *zpad_jmptab; zpad_jmptab++);
			zpad_jmptab++;
		}
	}

/* Now see what FFT length we would use a DWT does the k*b^n+c modulo. */
/* with a zero padded FFT. */

/* Use the proper 2^N+1 or 2^N-1 jmptable */

	if (c < 0) {
		if (CPU_FLAGS & CPU_SSE2) jmptab = INFT[0];
		else jmptab = INFT[2];
	} else {
		if (CPU_FLAGS & CPU_SSE2) jmptab = INFT[1];
		else jmptab = INFT[3];
	}

/* Find the table entry using either the specified fft length or */
/* the that can handle the k,b,n,c being tested. */

	while ((max_exp = jmptab[0]) != 0) {

/* Check L2 cache size constraints */
		
		if (l2_cache_size < ((jmptab[8] >> 16) & 0x7FFF))
			goto next2;

/* Check FFT requires prefetch capability */

		if (jmptab[8] & 0x80000000 && ! (CPU_FLAGS & CPU_PREFETCH))
			goto next2;

/* Check if this table entry matches the specified FFT length. */

		if (fftlen) {
			if (fftlen == jmptab[1]) break;
		}

/* Or check that this FFT length will work with this k,n,c pair */

		else {
			double max_bits_per_word;
			double bits_per_word;

/* Compute the maximum number of bits allowed in the FFT input word */

			max_bits_per_word = (double) max_exp / jmptab[1];
			max_bits_per_word -= GWSAFETY_MARGIN;

/* For historical reasons, the jmptable computes maximum exponent based on */
/* a Mersenne-mod FFT (i.e k=1.0, c=-1).  Handle more complex cases here. */
/* A Mersenne-mod FFT produces 2 * bits_per_word in each FFT result word. */
/* The more general case yields 2 * bits_per_word + log2(k) + 1.5 * log2(c) */
/* in each FFT result word. */

			bits_per_word = (log2k + n) / jmptab[1];
			if (2.0 * bits_per_word + log2k + 1.5 * log2c <
					2.0 * max_bits_per_word) {
				double total_bits, loglen;

/* Because carries are spread over 4 words, there is a minimum limit on */
/* the bits per word.  An FFT result word cannot be more than 5 times */
/* bits-per-word (bits-per-word are stored in the current word and the */
/* 4 words we propogate carries to.  How many bits are in an FFT result */
/* word?  Well, because of balanced representation the abs(input word) is */
/* (bits_per_word-1) bits long. An FFT result word contains multiplied data */
/* words, that's (bits_per_word-1)*2 bits.  Adding up many multiplied data */
/* words adds some bits proportional to the size of the FFT.  Experience */
/* has shown this to be 0.6 * log (FFTLEN).  This entire result is */
/* multiplied by k in the normalization code, so add another log2(k) bits. */
/* Finally, the mulbyconst adds to the size of the carry. */

				loglen = log ((double) jmptab[1]) / log (2.0);
				total_bits = (bits_per_word - 1.0) * 2.0 +
					     1.5 * log2c + loglen * 0.6 +
					     log2k + log2maxmulbyconst;
				if (total_bits > 5.0 * bits_per_word) {
					ASSERTG (zpad_jmptab == NULL ||
						 jmptab[1] >= zpad_jmptab[1]);
					goto next2;
				}

/* Because of limitations in the top_carry_adjust code, there is a limit */
/* on the size of k that can be handled.  This isn't a big deal since the */
/* zero-padded implementation will use the same FFT length.  Check to see */
/* if this is this k can be handled.  K must fit in the top three words */
/* for one-pass FFTs and within the top two words of two-pass FFTs. */

				if (jmptab[10] == 0 &&
				    log2k > floor (3.0 * bits_per_word)) {
					ASSERTG (zpad_jmptab == NULL ||
						 jmptab[1] >= zpad_jmptab[1]);
					goto next2;
				}
				if (jmptab[10] != 0 &&
				    log2k > floor (2.0 * bits_per_word)) {
					ASSERTG (zpad_jmptab == NULL ||
						 jmptab[1] >= zpad_jmptab[1]);
					goto next2;
				}
				break;
			}
		}

/* Move to next jmptable entry */

next2:		for (jmptab += 14; *jmptab; jmptab++);
		jmptab++;
	}

/* If the zero pad FFT length is less than the DWT FFT length, then use */
/* the zero pad FFT length. */

	if (zpad_jmptab != NULL && zpad_jmptab[0] &&
	    (jmptab[0] == 0 || zpad_jmptab[1] < jmptab[1])) {
		INFT[0] = zpad_jmptab;
		return (TRUE);
	}

/* If we found a DWT table entry then return the address in INFT. */

	if (jmptab[0]) {
		INFT[0] = jmptab;
		return (FALSE);
	}

/* Error - neither method could handle this huge number */

	INFT[0] = NULL;
	return (GWERROR_TOO_LARGE);
}


/* Allocate memory and initialize assembly code for arithmetic */
/* modulo k*b^n+c */

int gwsetup (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c,		/* C in K*B^N+C. */
	unsigned long fftlen)	/* Specific FFT size to use (or zero) */
{
	int	gcd, error_code;

/* Sanity check the k value */

	if (k < 1.0) return (GWERROR_K_TOO_SMALL);
	if (k > 18014398509481983.0) return (GWERROR_K_TOO_LARGE);

/* Our code fast code fails if k and c are not relatively prime.  This */
/* is because we cannot calculate 1/k.  Although the user shouldn't call */
/* us with this case, we handle it anyway by reverting to the slow general */
/* purpose multiply routines. */

	if (k == 1.0 || abs (c) == 1)
		gcd = 1;
	else {
		giant	kg, cg;
		kg = newgiant (4);
		cg = newgiant (4);
		dbltog (k, kg);
		itog (abs (c), cg);
		gcdg (kg, cg);
		gcd = cg->n[0];
		free (kg);
		free (cg);
	}

/* Call the internal setup routine when we can.  B must be 2, the */
/* gcd (k, c) must be 1, n must big enough so that their aren't too few */
/* bits per FFT word, also k * mulbyconst and c * mulbyconst cannot be too */
/* large.  Turn off flag indicating general-purpose modulos are being */
/* performed. */

	if (b == 2 && gcd == 1 && n >= 350 &&
	    k * GWMAXMULBYCONST <= MAX_ZEROPAD_K &&
	    abs (c) * GWMAXMULBYCONST <= MAX_ZEROPAD_C) {
		error_code = internal_gwsetup (k, b, n, c, fftlen);
		if (error_code) return (error_code);
		GENERAL_MOD = 0;
	}

/* Emulate b != 2, k not relatively prime to c, small n values, and */
/* large k or c values with a call to the general purpose modulo setup code. */

	else {
		double	bits;
		giant	g;

		bits = log ((double) b) / log (2.0) * (double) n;
		g = newgiant (((unsigned long) bits >> 4) + 8);
		ultog (b, g);
		power (g, n);
		dblmulg (k, g);
		iaddg (c, g);
		error_code = gwsetup_general_mod_giant (g, fftlen);
		free (g);
		if (error_code) return (error_code);
	}

/* For future messages, format the input number as a string */

	gw_as_string (GWSTRING_REP, k, b, n, c);

/* Return success */

	return (0);
}

/* This setup routine is for operations modulo an arbitrary binary number. */
/* This is three times slower than the special forms above. */
/* Only choose a specific FFT size if you know what you are doing!! */

int gwsetup_general_mod (
	unsigned long *array,	/* The modulus as an array of longs */
	unsigned long arraylen,	/* The number of longs in the array */
	unsigned long fftlen)	/* Zero or specific FFT size to use. */
{
	giantstruct tmp;
	tmp.sign = arraylen;
	tmp.n = array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	return (gwsetup_general_mod_giant (&tmp, fftlen));
}

int gwsetup_general_mod_giant (
	giant	n,		/* The modulus */
	unsigned long fftlen)	/* Zero or specific FFT size to use. */
{
#define EB	10		/* Extra bits of precision to compute quot. */
	unsigned long len;	/* Bit length of modulus */
	giant	tmp;

/* Setup the FFT code, use an integral number of bits per word if possible. */
/* We reserve some extra bits for extra precision and to make sure we can */
/* zero an integral number of words during copy. */

	len = bitlen (n);
	bit_length = len;
	gwsetup_without_mod (len + len + 2*EB + 64, fftlen);

/* Copy the modulus */

	GW_MODULUS = newgiant ((len >> 4) + 1);
	gtog (n, GW_MODULUS);

/* Remember the modulus.  FFT it for faster use. */

	GW_MODULUS_FFT = gwalloc ();
	gianttogw (n, GW_MODULUS_FFT);
	gwfft (GW_MODULUS_FFT, GW_MODULUS_FFT);

/* Precompute the reciprocal */

	tmp = newgiant ((PARG >> 4) + 1);
	itog (1, tmp);
	gshiftleft (len + len + EB, tmp);
	divg (n, tmp);		/* computes len+EB+1 bits of reciprocal */
	gshiftleft (PARG - len - len - EB, tmp);
				/* shift so gwmul routines wrap */
				/* quotient to lower end of fft */
	GW_RECIP = gwalloc ();
	gianttogw (tmp, GW_RECIP);
	gwfft (GW_RECIP, GW_RECIP);
	free (tmp);

/* Calculate number of words to zero during copy */

	if (len < EB) GW_ZEROWORDSLOW = 0;
	else GW_ZEROWORDSLOW = (unsigned long)
		((double) (len - EB) / fft_bits_per_word);

/* Set flag indicating general-purpose modulo operations are in force */

	GENERAL_MOD = 1;

/* Create dummy string representation. Calling gtoc to get the first */
/* several digits would be better, but it is too slow. */

	sprintf (GWSTRING_REP, "A %ld-bit number", len);

/* Return success */

	return (0);
}

/* This setup routine is for operations without a modulo. In essence, */
/* you are using gwnums as a general-purpose FFT multiply library. */
/* Only choose a specific FFT size if you know what you are doing!! */

int gwsetup_without_mod (
	unsigned long n,	/* Maximum number of bits in OUTPUT numbers. */
	unsigned long fftlen)	/* Zero or specific FFT size to use. */
{
	unsigned long *info, max_exponent, desired_n;
	int	error_code;

/* Call gwinfo and have it figure out the FFT length we will use. */
/* Since the user must zero the upper half of FFT input data, the FFT */
/* outputs will be smaller.  This lets us get about another 0.3 bits */
/* per input word. */

	GWSAFETY_MARGIN -= 0.3;
	error_code = gwinfo (1.0, 2, n, -1, fftlen);
	GWSAFETY_MARGIN += 0.3;

	ASSERTG (error_code != 1);
	if (error_code) return (error_code);

	info = INFT[0];
	max_exponent = info[0];
	fftlen = info[1];

/* If possible, increase n to the next multiple of FFT length.  This is */
/* because rational FFTs are faster than irrational FFTs (no FFT weights). */

	desired_n = ((n + fftlen - 1) / fftlen) * fftlen;
	if (desired_n < max_exponent) n = desired_n;

/* Our FFTs don't handle cases where there are few bits per word because */
/* carries must be propagated over too many words.  Arbitrarily insist */
/* that n is at least 12 * fftlen.  */

	if (n < 12 * fftlen) n = 12 * fftlen;

/* Now setup the assembly code */

	error_code = gwsetup (1.0, 2, n, -1, fftlen);
	ASSERTG (error_code == 0);

/* Set flag indicating general-purpose modulo operations are not in force */

	GENERAL_MOD = 0;

/* Create dummy string representation. */

	strcpy (GWSTRING_REP, "No modulus");

/* Return success */

	return (0);
}


/* Common setup routine for the three different user-visible setup routines */
/* Allocate memory and initialize assembly code for arithmetic */
/* modulo k*b^n+c */

int internal_gwsetup (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. Must be two. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c,		/* C in K*B^N+C. Must be rel. prime to K. */
	unsigned long fftlen)	/* Specific FFT size to use (or zero) */
{
	int	error_code;
	unsigned long mem_needed;
	unsigned long *info;
	double	fft_bit_length;		/* Bit length of the FFT */
	double	*tables;		/* Pointer tables we are building */
	unsigned long pass1_size, pass2_size;

	ASSERTG (FFTLEN == 0);

/* Select the proper FFT size for this k,n,c combination */

	error_code = gwinfo (k, b, n, c, fftlen);
	if (error_code > 1) return (error_code);
	ZERO_PADDED_FFT = error_code;

/* Get pointer to fft info and allocate needed memory */

	fpu_init ();
	info = INFT[0];
	mem_needed = info[3];
	gwnum_memory = tables = (double *) aligned_malloc (mem_needed, 4096);

/* Do a seemingly pointless memset! */
/* The memset will walk through the allocated memory sequentially, which */
/* increases the liklihood that contiguous virtual memory will map to */
/* contiguous physical memory. */

	memset (tables, 0, mem_needed);

/* Setup some useful global variables */

	KARG = k;
	BARG = b;
	PARG = n;
	CARG = c;
	FFTLEN = info[1];

/* Initialize the extended precision code that computes the FFT weights */

	gwfft_weight_setup (ZERO_PADDED_FFT, k, n, c, FFTLEN);

/* Calculate the number of bits in k*2^n.  This will be helpful in */
/* determining how much meory to allocate for giants. */

	bit_length = log (k) / log ((double) 2.0) + n;

/* Calculate the number of bits the underlying FFT computes.  That is, */
/* the point at which data wraps around to the low FFT word.  For a zero */
/* pad FFT, this is simply 2*n.  Otherwise, it is log2(k) + n. */

	fft_bit_length = ZERO_PADDED_FFT ? n * 2.0 : bit_length;

/* Calculate the average number of bits in each FFT word. */

	fft_bits_per_word = fft_bit_length / FFTLEN;

/* Calculate the number of bits in each small FFT word. */

	BITS_PER_WORD = (unsigned long) fft_bits_per_word;

/* Set a flag if this is a rational FFT.  That is, an FFT where all the */
/* weighting factors are 1.0.  This happens when c is -1 and the */
/* fft_bit_length is a multiple of FFTLEN.  The assembly code can make some */
/* obvious optimizations when all the FFT weights are one. */

	RATIONAL_FFT = ((double) BITS_PER_WORD == fft_bits_per_word) && (c == -1);

/* Remember the maximum number of bits per word that this FFT length */
/* supports.  We this in gwnear_fft_limit.  Note that zero padded FFTs */
/* can support an extra 0.3 bits per word because of the all the zeroes. */

	fft_max_bits_per_word = (double) info[0] / (double) FFTLEN;
	if (ZERO_PADDED_FFT) fft_max_bits_per_word += 0.3;

/* Compute extra bits - the number of adds we can tolerate without */
/* a normalization operation. Under normal circumstances, max_bits */
/* will be greater than virtual bits, but playing with the safety margin */
/* or forcing use of a specific FFT length could change that. */

	EXTRA_BITS = (unsigned long)
		pow (2.0, (fft_max_bits_per_word - virtual_bits_per_word ()) / 2.0);

/* This debug code will make sure the memory allocated is appropriate */
/* for the tables we are about to init. */

#ifdef GDEBUG
	{double *t1=tables;
#endif

/* See how many cache lines are grouped in pass 1.  This will affect how */
/* we build the normalization tables.  Note that cache line sizes are */
/* different in the x87 (16 bytes) and SSE2 code (64 bytes). */

	PASS1_CACHE_LINES = (info[8] & 0xFFFF);

/* Determine the pass 1 & pass 2 sizes.  This affects how we build */
/* many of the sin/cos tables. */

	PASS2_LEVELS = info[10];	/* Num FFT levels done in pass2 */
	pass2_size = 1 << PASS2_LEVELS;	/* Complex values in pass2 section */
	pass1_size = FFTLEN / pass2_size; /* Real values in a pass1 section */

/* Remember the size of the scratch area */

	SCRATCH_SIZE = info[9];

/* Initialize tables for the SSE2 assembly code. */

	if (CPU_FLAGS & CPU_SSE2) {

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with the scratch area which is also used in the first pass. */

		if (pass2_size > 1) {
			int	i, carry_table_size;
			double	xmm_bigval;
			((double **)GWPROCPTRS)[15] = tables;
			carry_table_size = FFTLEN / (pass2_size / 2);
			xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
			for (i = 0; i < carry_table_size; i++)
				*tables++ = xmm_bigval;
		}

/* Reserve room for the pass 1 scratch area. */

		((double**)GWPROCPTRS)[16] = tables;
		if (SCRATCH_SIZE)
			tables = (double *) ((char *) tables + SCRATCH_SIZE);

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		((double **)GWPROCPTRS)[12] = tables;
		tables = build_norm_table (tables, pass2_size, 0);

/* Build sin/cos tables used in pass 1.  If FFTLEN is a power of two, */
/* many of the sin/cos tables can be shared. */

		((double **)GWPROCPTRS)[2] = tables;
		tables = build_sin_cos_table (tables, pass1_size, c < 0 || ZERO_PADDED_FFT, pass2_size == 1 ? 2 : 1);

		if (pass2_size > 1 && pass1_size == pow_two_above_or_equal (pass1_size))
			GWPROCPTRS[3] = GWPROCPTRS[2];
		else {
			((double **)GWPROCPTRS)[3] = tables;
			tables = build_sin_cos_table (tables, pass1_size/4, c < 0 || ZERO_PADDED_FFT, 1);
		}

		if (pass1_size == pow_two_above_or_equal (pass1_size)) {
			GWPROCPTRS[4] = GWPROCPTRS[3];
			GWPROCPTRS[5] = GWPROCPTRS[3];
			GWPROCPTRS[6] = GWPROCPTRS[3];
		} else {
			((double **)GWPROCPTRS)[4] = tables;
			tables = build_sin_cos_table (tables, pass1_size/16, c < 0 || ZERO_PADDED_FFT, 1);
			((double **)GWPROCPTRS)[5] = tables;
			tables = build_sin_cos_table (tables, pass1_size/64, c < 0 || ZERO_PADDED_FFT, 1);
			((double **)GWPROCPTRS)[6] = tables;
			tables = build_sin_cos_table (tables, pass1_size/256, c < 0 || ZERO_PADDED_FFT, 1);
		}

/* Build sin/cos and premultiplier tables used in pass 2 of two pass FFTs */
/* Remember that pass2_size is the number of complex values in a pass 2 */
/* section, but build_sin_cos_table wants the number of reals in a section. */
/* Also, in SSE2 the first complex section is handled in the real code */
/* by handling a double-sized section (when compared to the x87 code's */
/* pass 2 real section.  For example, the pass2 8 levels real section */
/* handles 512 doubles in x87, but handles 1024 doubles in SSE2. */

		if (pass2_size > 1) {
			((double **)GWPROCPTRS)[0] = tables;
			tables = build_premult_table (tables, pass2_size);
			((double **)GWPROCPTRS)[1] = tables;
			tables = build_sin_cos_table (tables, pass2_size*2, 0, 1);

			if (c < 0 || ZERO_PADDED_FFT) {
				((double **)GWPROCPTRS)[7] = tables;
				tables = build_sin_cos_table (tables, pass2_size * 4, 1, 2);
				((double **)GWPROCPTRS)[8] = tables;
				tables = build_sin_cos_table (tables, pass2_size, 1, 1);
			}

//			if (pass1_size == pow_two_above_or_equal (pass1_size)) {
				GWPROCPTRS[9] = GWPROCPTRS[8];
				GWPROCPTRS[10] = GWPROCPTRS[8];
				GWPROCPTRS[11] = GWPROCPTRS[8];
//			} else {
//				((double **)GWPROCPTRS)[9] = tables;
//				tables = build_sin_cos_table (tables, pass2_size/4, c < 0 || ZERO_PADDED_FFT, 1);
//				((double **)GWPROCPTRS)[10] = tables;
//				tables = build_sin_cos_table (tables, pass2_size/16, c < 0 || ZERO_PADDED_FFT, 1);
//				((double **)GWPROCPTRS)[11] = tables;
//				tables = build_sin_cos_table (tables, pass2_size/64, c < 0 || ZERO_PADDED_FFT, 1);
//			}
		}

/* Build the plus1-pre-multiplier table (complex weights applied when c > 0 */
/* and we are doing a all-complex FFT rather than emulating it with a */
/* zero-padded FFT. */

		if (c > 0 && !ZERO_PADDED_FFT) {
			((double **)GWPROCPTRS)[17] = tables;
			tables = build_plus1_table (tables, pass2_size);
		}

/* Build the column normalization multiplier table. */

		((double **)GWPROCPTRS)[13] = tables;
		tables = build_norm_table (tables, pass2_size, 1);

/* Build the table of big vs. little flags. */

		((double **)GWPROCPTRS)[14] = tables;
		tables = build_biglit_table (tables, pass2_size);
	}

/* Initialze table for the x87 assembly code. */

	if (! (CPU_FLAGS & CPU_SSE2)) {

/* Allocate a table for carries.  Init with zero.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with the scratch area which is also used in the first pass. */

		if (pass2_size > 1) {
			int	i, carry_table_size;
			((double **)GWPROCPTRS)[15] = tables;
			carry_table_size = FFTLEN / pass2_size;
			for (i = 0; i < carry_table_size; i++) *tables++ = 0.0;
		}

/* Reserve room for the pass 1 scratch area. */

		((double**)GWPROCPTRS)[16] = tables;
		if (SCRATCH_SIZE)
			tables = (double *) ((char *) tables + SCRATCH_SIZE);

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		((double **)GWPROCPTRS)[12] = tables;
		tables = build_x87_norm_table (tables, pass2_size, 0);

/* Build sin/cos tables used in pass 1.  If FFTLEN is a power of two, */
/* many of the sin/cos tables can be shared. */

		((double **)GWPROCPTRS)[2] = tables;
		tables = build_x87_sin_cos_table (tables, pass1_size, c < 0 || ZERO_PADDED_FFT);

		if (pass1_size == pow_two_above_or_equal (pass1_size)) {
			GWPROCPTRS[3] = GWPROCPTRS[2];
			GWPROCPTRS[4] = GWPROCPTRS[2];
			GWPROCPTRS[5] = GWPROCPTRS[2];
			GWPROCPTRS[6] = GWPROCPTRS[2];
		} else {
			((double **)GWPROCPTRS)[3] = tables;
			tables = build_x87_sin_cos_table (tables, pass1_size/4, c < 0 || ZERO_PADDED_FFT);
			((double **)GWPROCPTRS)[4] = tables;
			tables = build_x87_sin_cos_table (tables, pass1_size/16, c < 0 || ZERO_PADDED_FFT);
			((double **)GWPROCPTRS)[5] = tables;
			tables = build_x87_sin_cos_table (tables, pass1_size/64, c < 0 || ZERO_PADDED_FFT);
			((double **)GWPROCPTRS)[6] = tables;
			tables = build_x87_sin_cos_table (tables, pass1_size/256, c < 0 || ZERO_PADDED_FFT);
		}

/* Build sin/cos and premultiplier tables used in pass 2 of two pass FFTs */
/* Remember that pass2_size is the number of complex values in a pass 2 */
/* section, but build_x87_sin_cos_table wants the number of reals in */
/* a section. */

		if (pass2_size > 1) {
			((double **)GWPROCPTRS)[0] = tables;
			tables = build_x87_premult_table (tables, pass2_size);
			((double **)GWPROCPTRS)[1] = tables;
			tables = build_x87_sin_cos_table (tables, pass2_size*2, 0);

			if (c < 0 || ZERO_PADDED_FFT) {
				((double **)GWPROCPTRS)[7] = tables;
				tables = build_x87_sin_cos_table (tables, pass2_size*2, 1);
				GWPROCPTRS[8] = GWPROCPTRS[7];
				GWPROCPTRS[9] = GWPROCPTRS[7];
				GWPROCPTRS[10] = GWPROCPTRS[7];
				GWPROCPTRS[11] = GWPROCPTRS[7];
			}
		}

/* Build the plus1-pre-multiplier table (complex weights applied when c > 0 */
/* and we are doing a all-complex FFT rather than emulating it with a */
/* zero-padded FFT. */

		if (c > 0 && !ZERO_PADDED_FFT) {
			((double **)GWPROCPTRS)[17] = tables;
			tables = build_x87_plus1_table (tables, pass2_size);
		}

/* Build the column normalization multiplier table. */

		((double **)GWPROCPTRS)[13] = tables;
		tables = build_x87_norm_table (tables, pass2_size, 1);

/* Build the table of big vs. little flags. */

		((double **)GWPROCPTRS)[14] = tables;
		tables = build_x87_biglit_table (tables, pass2_size);
	}

/* Finish verifying table size */

#ifdef GDEBUG
	{char buf[80];
	long mem = (int) tables - (int) t1;
	if (mem != mem_needed) {
		sprintf (buf, "%d, mem: %d\n", FFTLEN, mem);
		OutputBoth(buf);}}}
#endif

/* Now call assembly routine to finish off the initialization */

	gwsetup2 ();

/* If the carry must be spread over more than 2 words, then set global */
/* so that assembly code knows this.  In theory, we could study what */
/* values of k and c can also use the 2 word carry propagation.  This */
/* isn't a major performance gain. */

	if (ZERO_PADDED_FFT || (k == 1.0 && abs (c) == 1))
		SPREAD_CARRY_OVER_4_WORDS = FALSE;
	else
		SPREAD_CARRY_OVER_4_WORDS = TRUE;

/* Set some global variables that make life easier in the assembly code */
/* that wraps carry out of top FFT word into the bottom FFT word. */
/* This is needed when k > 1 and we are not doing a zero padded FFT. */

	TOP_CARRY_NEEDS_ADJUSTING = (KARG > 1.0 && !ZERO_PADDED_FFT);
	if (TOP_CARRY_NEEDS_ADJUSTING) {
		unsigned long kbits, kbits_lo;
		unsigned long topwordbits, secondwordbits, thirdwordbits;

/* Invert KARG and split KARG for computing top carry adjustment without */
/* precision problems. */

		INVERSE_KARG = 1.0 / k;
		kbits = (unsigned long) ceil (bit_length) - n;
		kbits_lo = kbits / 2;
		KARG_HI = ((unsigned long) k) & ~((1 << kbits_lo) - 1);
		KARG_LO = ((unsigned long) k) &  ((1 << kbits_lo) - 1);

/* Calculate top carry adjusting constants */

		topwordbits = BITS_PER_WORD;
		if (is_big_word (FFTLEN-1)) topwordbits++;
		secondwordbits = BITS_PER_WORD;
		if (is_big_word (FFTLEN-2)) secondwordbits++;
		thirdwordbits = BITS_PER_WORD;
		if (is_big_word (FFTLEN-3)) thirdwordbits++;

		CARRY_ADJUST1 = (double) (1 << kbits);
		CARRY_ADJUST2 = (double) (1 << topwordbits) / (double) (1 << kbits);
		CARRY_ADJUST3 = gwfft_weight (FFTLEN-1);

/* Get the addr of the top three words.  This is funky because in two-pass */
/* FFTs we want the scratch area offset when normalizing after a multiply, */
/* but the FFT data when normalizing after an add/sub.  For one-pass FFTs, */
/* we always want the FFT data offset. */

		HIGH_WORD1_OFFSET = addr_offset (FFTLEN, FFTLEN-1);
		HIGH_WORD2_OFFSET = addr_offset (FFTLEN, FFTLEN-2);
		HIGH_WORD3_OFFSET = addr_offset (FFTLEN, FFTLEN-3);

		raw_gwsetaddin (FFTLEN-1, 0);
		HIGH_SCRATCH1_OFFSET = ADDIN_OFFSET;
		raw_gwsetaddin (FFTLEN-2, 0);
		HIGH_SCRATCH2_OFFSET = ADDIN_OFFSET;
		raw_gwsetaddin (FFTLEN-3, 0);
		HIGH_SCRATCH3_OFFSET = ADDIN_OFFSET;

/* In two-pass FFTs, we only support tweaking the top two words.  Compute */
/* the necessary constants. */

		if (PASS2_LEVELS) {
			ASSERTG (kbits <= topwordbits + secondwordbits);
			CARRY_ADJUST4 = (double) (1 << secondwordbits) *
							gwfft_weight (FFTLEN-2);
		}

/* In one-pass FFTs, we adjust the top three words.  More adjustment */
/* variables are needed. */

		else {
			ASSERTG (kbits <= topwordbits + secondwordbits + thirdwordbits);
			CARRY_ADJUST4 = (double) (1 << secondwordbits);
			CARRY_ADJUST5 = gwfft_weight (FFTLEN-2);
			CARRY_ADJUST6 = (double) (1 << thirdwordbits) *
							gwfft_weight (FFTLEN-3);
		}
	}

/* Set some global variables that make life easier in the assembly code */
/* that handles zero padded FFTs. */

	if (ZERO_PADDED_FFT) {
		unsigned long kbits, bits0, bits1, bits2, bits3, bits4, bits5;
		double	pow2, bigpow2;

		HIGH_WORD1_OFFSET = addr_offset (FFTLEN, FFTLEN/2-1);
		HIGH_WORD2_OFFSET = addr_offset (FFTLEN, FFTLEN/2-2);
		HIGH_WORD3_OFFSET = addr_offset (FFTLEN, FFTLEN/2-3);

		raw_gwsetaddin (FFTLEN/2-1, 0);
		HIGH_SCRATCH1_OFFSET = ADDIN_OFFSET;
		raw_gwsetaddin (FFTLEN/2-2, 0);
		HIGH_SCRATCH2_OFFSET = ADDIN_OFFSET;
		raw_gwsetaddin (FFTLEN/2-3, 0);
		HIGH_SCRATCH3_OFFSET = ADDIN_OFFSET;

		kbits = (unsigned long) ceil (bit_length) - n;
		bits0 = BITS_PER_WORD; if (is_big_word (0)) bits0++;
		bits1 = BITS_PER_WORD; if (is_big_word (1)) bits1++;
		bits2 = BITS_PER_WORD; if (is_big_word (2)) bits2++;
		bits3 = BITS_PER_WORD; if (is_big_word (3)) bits3++;
		bits4 = BITS_PER_WORD; if (is_big_word (4)) bits4++;
		bits5 = BITS_PER_WORD; if (is_big_word (5)) bits5++;

		ZPAD_SHIFT1 = pow ((double) 2.0, (int) bits0);
		ZPAD_SHIFT2 = pow ((double) 2.0, (int) bits1);
		ZPAD_SHIFT3 = pow ((double) 2.0, (int) bits2);
		ZPAD_SHIFT4 = pow ((double) 2.0, (int) bits3);
		ZPAD_SHIFT5 = pow ((double) 2.0, (int) bits4);
		ZPAD_SHIFT6 = pow ((double) 2.0, (int) bits5);

		if (kbits <= BITS_PER_WORD + 3) ZPAD_TYPE = 1;
		else if (kbits <= 2 * BITS_PER_WORD + 3) ZPAD_TYPE = 2;
		else ZPAD_TYPE = 3;

		if (ZPAD_TYPE == 1) {
			ZPAD_K1_LO = k;
			ZPAD_INVERSE_K1 = 1.0 / k;
		}

		if (ZPAD_TYPE == 2) {
			ZPAD_K1_HI = floor (k / ZPAD_SHIFT1);
			ZPAD_K1_LO = k - ZPAD_K1_HI * ZPAD_SHIFT1;
			ZPAD_INVERSE_K1 = ZPAD_SHIFT1 / k;
			ZPAD_K2_HI = floor (k / ZPAD_SHIFT2);
			ZPAD_K2_LO = k - ZPAD_K2_HI * ZPAD_SHIFT2;
			ZPAD_INVERSE_K2 = ZPAD_SHIFT2 / k;
			ZPAD_K3_HI = floor (k / ZPAD_SHIFT3);
			ZPAD_K3_LO = k - ZPAD_K3_HI * ZPAD_SHIFT3;
			ZPAD_INVERSE_K3 = ZPAD_SHIFT3 / k;
			ZPAD_K4_HI = floor (k / ZPAD_SHIFT4);
			ZPAD_K4_LO = k - ZPAD_K4_HI * ZPAD_SHIFT4;
			ZPAD_INVERSE_K4 = ZPAD_SHIFT4 / k;
			ZPAD_K5_HI = floor (k / ZPAD_SHIFT5);
			ZPAD_K5_LO = k - ZPAD_K5_HI * ZPAD_SHIFT5;
			ZPAD_INVERSE_K5 = ZPAD_SHIFT5 / k;
			ZPAD_K6_HI = floor (k / ZPAD_SHIFT6);
			ZPAD_K6_LO = k - ZPAD_K6_HI * ZPAD_SHIFT6;
			ZPAD_INVERSE_K6 = ZPAD_SHIFT6 / k;
		}

		if (ZPAD_TYPE == 3) {
			pow2 = pow ((double) 2.0, (int) bits0);
			bigpow2 = pow ((double) 2.0, (int) (bits0 + bits1));
			ZPAD_K2_HI = floor (k / bigpow2);
			ZPAD_K2_MID = floor ((k - ZPAD_K2_HI*bigpow2) / pow2);
			ZPAD_K2_LO = k - ZPAD_K2_HI*bigpow2 - ZPAD_K2_MID*pow2;
			ZPAD_INVERSE_K2 = pow2 / k;
			pow2 = pow ((double) 2.0, (int) bits1);
			bigpow2 = pow ((double) 2.0, (int) (bits1 + bits2));
			ZPAD_K3_HI = floor (k / bigpow2);
			ZPAD_K3_MID = floor ((k - ZPAD_K3_HI*bigpow2) / pow2);
			ZPAD_K3_LO = k - ZPAD_K3_HI*bigpow2 - ZPAD_K3_MID*pow2;
			ZPAD_INVERSE_K3 = pow2 / k;
			pow2 = pow ((double) 2.0, (int) bits2);
			bigpow2 = pow ((double) 2.0, (int) (bits2 + bits3));
			ZPAD_K4_HI = floor (k / bigpow2);
			ZPAD_K4_MID = floor ((k - ZPAD_K4_HI*bigpow2) / pow2);
			ZPAD_K4_LO = k - ZPAD_K4_HI*bigpow2 - ZPAD_K4_MID*pow2;
			ZPAD_INVERSE_K4 = pow2 / k;
			pow2 = pow ((double) 2.0, (int) bits3);
			bigpow2 = pow ((double) 2.0, (int) (bits3 + bits4));
			ZPAD_K5_HI = floor (k / bigpow2);
			ZPAD_K5_MID = floor ((k - ZPAD_K5_HI*bigpow2) / pow2);
			ZPAD_K5_LO = k - ZPAD_K5_HI*bigpow2 - ZPAD_K5_MID*pow2;
			ZPAD_INVERSE_K5 = pow2 / k;
			pow2 = pow ((double) 2.0, (int) bits4);
			bigpow2 = pow ((double) 2.0, (int) (bits4 + bits5));
			ZPAD_K6_HI = floor (k / bigpow2);
			ZPAD_K6_MID = floor ((k - ZPAD_K6_HI*bigpow2) / pow2);
			ZPAD_K6_LO = k - ZPAD_K6_HI*bigpow2 - ZPAD_K6_MID*pow2;
			ZPAD_INVERSE_K6 = bigpow2 / k;
		}
	}

/* Point to default normalization routines */

	gwsetnormroutine (0, 0, 0);
	POSTFFT = FALSE;
	raw_gwsetaddin (0, 0);

/* Clear globals */

	MAXERR = 0.0;
	GWERROR = 0;
	COPYZERO[0] = 0;
	GW_RANDOM = NULL;

/* Compute maximum allowable difference for error checking */
/* This error check is disabled for mod 2^N+1 arithmetic */

	if (!ZERO_PADDED_FFT && CARG > 0)
		MAXDIFF = 1.0E80;

/* We have observed that the difference seems to vary based on the size */
/* the FFT result word.  This is two times the number of bits per double. */
/* Subtract 1 from bits per double because one bit is the sign bit. */
/* Add in a percentage of the log(FFTLEN) to account for carries. */
/* We use a different threshold for SSE2 which uses 64-bit instead of */
/* 80-bit doubles during the FFT */

	else {
		double bits_per_double, total_bits, loglen;
		bits_per_double = fft_bits_per_word - 1.0;
		if (!ZERO_PADDED_FFT) bits_per_double += log ((double) -c) / log ((double) 2.0);
		loglen = log ((double) FFTLEN) / log ((double) 2.0);
		loglen *= 0.69;
		total_bits = bits_per_double * 2.0 + loglen * 2.0;
		MAXDIFF = pow ((double) 2.0, total_bits -
				((CPU_FLAGS & CPU_SSE2) ? 47.08 : 47.65));
	}

/* Clear counters */

	fft_count = 0;

/* Default size of gwnum_alloc array is 50 */

	gwnum_alloc = NULL;
	gwnum_alloc_count = 0;
	gwnum_alloc_array_size = 50;
	gwnum_free = NULL;
	gwnum_free_count = 0;

/* Compute alignment for allocated data.  Strangely enough assembly */
/* prefetching works best in pass 1 on a P4 if the data is allocated */
/* on an odd cache line.  An optimal 31 of the 32 cache lines on a 4KB */
/* page will be prefetchable.  Page aligned data would only prefetch */
/* 28 of the 32 cache lines. */

	if (CPU_FLAGS & CPU_SSE2) {
		if (PASS2_LEVELS == 0) {	/* One pass */
			GW_ALIGNMENT = 128;	/* P4 cache line alignment */
			GW_ALIGNMENT_MOD = 0;
		} else if (SCRATCH_SIZE == 0) {	/* Small two passes */
			GW_ALIGNMENT = 4096;	/* Page alignment */
			GW_ALIGNMENT_MOD = 0;
		} else {			/* Large two passes */
			GW_ALIGNMENT = 512;	/* Clmblkdst + 1 cache line */
			GW_ALIGNMENT_MOD = 128;
		}
	} else {
		if (PASS2_LEVELS == 0)		/* One pass */
			GW_ALIGNMENT = 128;	/* P4 cache line alignment */
		else				/* Two passes */
			GW_ALIGNMENT = 4096;	/* Page alignment */
		GW_ALIGNMENT_MOD = 0;
	}

/* Return success */

	return (0);
}


/* Cleanup any memory allocated for multi-precision math */

void gwdone (void)
{
	unsigned int i;

	term_giants ();
	aligned_free (gwnum_memory);
	gwnum_memory = NULL;
	free (gwnum_free);
	gwnum_free = NULL;
	if (gwnum_alloc != NULL) {
		for (i = 0; i < gwnum_alloc_count; i++) {
			char	*p;
			long	freeable;
			p = (char *) gwnum_alloc[i];
			freeable = * (long *) (p - 32);
			if (freeable) aligned_free ((char *) p - GW_HEADER_SIZE);
		}
		free (gwnum_alloc);
		gwnum_alloc = NULL;
	}
	free (GW_MODULUS);
	GW_MODULUS = NULL;
	FFTLEN = 0;
}

/* Routine to allocate aligned memory for our big numbers */
/* Memory is allocated on 128-byte boundaries, with an additional */
/* 32 bytes prior to the data for storing useful stuff */

gwnum gwalloc (void)
{
	unsigned long size, aligned_size;
	char	*p, *q;

/* Return cached gwnum if possible */

	if (gwnum_free_count)
		return (gwnum_free[--gwnum_free_count]);

/* Allocate arrays if necessary */

	if (gwnum_alloc == NULL) {
		gwnum_free = (gwnum *)
			malloc (gwnum_alloc_array_size * sizeof (gwnum));
		if (gwnum_free == NULL) return (NULL);
		gwnum_alloc = (gwnum *)
			malloc (gwnum_alloc_array_size * sizeof (gwnum));
		if (gwnum_alloc == NULL) return (NULL);
	} else if (gwnum_alloc_count == gwnum_alloc_array_size) {
		gwnum_alloc_array_size += gwnum_alloc_array_size >> 1;
		gwnum_free = (gwnum *)
			realloc (gwnum_free,
				 gwnum_alloc_array_size * sizeof (gwnum));
		if (gwnum_free == NULL) return (NULL);
		gwnum_alloc = (gwnum *)
			realloc (gwnum_alloc,
				 gwnum_alloc_array_size * sizeof (gwnum));
		if (gwnum_alloc == NULL) return (NULL);
	}

/* Use addr function on the last FFT value to compute the size. */
/* Allocate 96 extra bytes for header information and align the data */
/* appropriately.  When allocating memory out of the big buffer for */
/* the torture test, then only allocate data on 128 byte boundaries */
/* to maximize the number of gwnums allocated. */

	size = gwnum_size (FFTLEN);
	aligned_size = (size + GW_HEADER_SIZE + 127) & ~127;
	if (GW_BIGBUF_SIZE >= size + aligned_size) {
		p = GW_BIGBUF;
		GW_BIGBUF += aligned_size;
		GW_BIGBUF_SIZE -= aligned_size;
	} else {
		p = (char *) aligned_offset_malloc (
				size + GW_HEADER_SIZE, GW_ALIGNMENT,
				GW_HEADER_SIZE + GW_ALIGNMENT_MOD);
		if (p == NULL) return (NULL);
	}

/* Set byte in header that indicates whether the memory should be freed. */

	q = p + GW_HEADER_SIZE;
	if (GW_BIGBUF_SIZE >= size + aligned_size)
		* (long *) (q - 32) = 0;	/* Don't free this memory */
	else
		* (long *) (q - 32) = 1;	/* Do free this memory */

/* Do a seemingly pointless memset!  This actual is very important. */
/* The memset will walk through the allocated memory sequentially, which */
/* increases the likelihood that contiguous virtual memory will map to */
/* contiguous physical memory.  The FFTs, especially the larger ones, */
/* optimizes L2 cache line collisions on the assumption that the FFT data */
/* is in contiguous physical memory.  Failure to do this results in as */
/* much as a 30% performance hit in an SSE2 2M FFT. */

	memset (q, 0, size);

/* Initialize the header */

	* (unsigned long *) (q - 8) = size;	/* Size in bytes */
	* (unsigned long *) (q - 4) = 0;	/* Unnormalized adds count */
	* (unsigned long *) (q - 28) = 0;	/* Has-been-pre-ffted flag */
	* (double *) (q - 16) = 0.0;
	* (double *) (q - 24) = 0.0;

/* Save pointer for easier cleanup */

	gwnum_alloc[gwnum_alloc_count++] = (gwnum) q;

/* Return the gwnum */

	return ((gwnum) q);
}

/* Free one of our special numbers */

void gwfree (
	gwnum	q)
{
	gwnum_free[gwnum_free_count++] = q;
}

/* Specialized routines that let the giants code share the free */
/* memory pool used by gwnums. */

void gwfree_temporarily (
	gwnum	q)
{
	gwfree (q);
}
void gwrealloc_temporarily (
	gwnum	q)
{
	unsigned long i, j;

	for (i = j = 0; i < gwnum_free_count; i++)
		if (gwnum_free[i] != q) gwnum_free[j++] = gwnum_free[i];
	gwnum_free_count = j;
}

/* Free all of our special numbers */

void gwfreeall (void)
{
	unsigned int i;
	if (gwnum_alloc == NULL) return;
	for (i = 0; i < gwnum_alloc_count; i++)
		gwnum_free[i] = gwnum_alloc[i];
	gwnum_free_count = gwnum_alloc_count;
}


void gwcopy (			/* Copy a gwnum */
	gwnum	s,		/* Source */
	gwnum	d)		/* Dest */
{
	unsigned long free_offset;

/* Load the one piece of information that should not be copied over */

	free_offset = ((unsigned long *) d)[-8];

/* Copy the data and 96-byte header */

	memcpy ((char *) d - GW_HEADER_SIZE,
		(char *) s - GW_HEADER_SIZE,
		((unsigned long *) s)[-2] + GW_HEADER_SIZE);

/* Restore the one piece of information that should not be copied over */

	((unsigned long *) d)[-8] = free_offset;
}

/* To optimize use of the L1 cache we scramble the FFT data. */
/* Note:  The Intel L1 data cache is 8KB two-way set associative with */
/* 32 byte cache lines.  Later CPUs have more cache, but we are prepared */
/* for the worst case.  This tiny cache will require us to perform */
/* three "passes" to perform a large FFT. Each pass must minimize */
/* L1 cache line conflicts - that is have no data at the same address */
/* modulo 4096 */

/* 1) We'd like to do as much work as possible in the final pass (called */
/*    pass 2 in a lot of the code).  Since some cache space is required */
/*    for sine/cosine data, we only use half of the L1 cache for FFT data. */
/*    4KB = 512 values = 256 complex values.  Thus, the final pass will */
/*    perform 8 FFT levels.  Also note that it will be advantageous to */
/*    have the real and imaginary values in the same cache line.  Thus, */
/*    the first cache line contains the 0th, 128th, 256th, and 384th FFT */
/*    data values.  Where the 0th and 128th values comprise a single */
/*    complex number as does the 256th and 384th. */

/* 2) To eliminate cache line conflicts in the middle pass (called pass 1 */
/*    in a lot of this code), 32 bytes is wasted after 4KB of FFT data */
/*    If we did not do this every pass 1 value would try to occupy the */
/*    same L1 cache line! */

/* Putting it all together, for FFTLEN=2^16 you get this memory layout:	*/
/*	0	128	256	384		(32 bytes)		*/
/*	1	129	257	385		(32 bytes)		*/
/*		   etc.							*/
/*	127	255	383	511		(32 bytes)		*/
/*		(32 wasted bytes)					*/
/*	512	640	768	896		(32 bytes)		*/
/*	513	641	769	897		(32 bytes)		*/
/*		   etc.							*/
/*	639	767	895	1023		(32 bytes)		*/
/*		(32 wasted bytes)					*/
/*	1024	1152	1280	1408		(32 bytes)		*/
/*		   etc.							*/

/* Well.... I implemented the above only to discover I had dreadful */
/* performance in pass 1.  How can that be?  The problem is that each  */
/* cache line in pass 1 comes from a different 4KB page.  Therefore, */
/* pass 1 accessed 128 different pages.  This is a problem because the */
/* Pentium chip has only 64 TLBs (translation lookaside buffers) to map */
/* logical page addresses into physical addresses.  So we need to shuffle */
/* the data further so that pass 1 data is on fewer pages while */
/* pass 2 data is spread over more pages. */

/* 1st 4KB page		2nd page	...	18th page	*/
/* 0 128 256 384	waste			waste		*/
/* 512 640 768 896	8 136 264 392		waste		*/
/* ...								*/
/* 7680 ...					waste		*/
/* 1 129 257 385				8192 8320 ...	*/
/* 513 ...            	9 137 265 393		8704 ...	*/
/*               ...						*/
/* 7 ...							*/
/* 519 ...							*/
/*	         ...						*/
/* 7687 ...							*/

/* That is, waste 32 bytes after each 512 FFT data values (4KB). */
/* Except after 8192 FFT data values go to the next 4KB page and waste */
/* the first 16*32 bytes.  If you look carefully at the above, you'll see */
/* that in pass 2 the FFT data (values 0 through 511) comes from the first */
/* 16 4KB pages (actually the waste bytes make this 17 4KB pages).  Similarly,
/* the pass 1 data (values 0 up to 65536 stepping by 256) comes from 16 */
/* different 4KB pages  ---  and there are no L1 cache line conflicts!!! */
/* Furthermore, when accessing pages, the pages are an odd number apart */
/* (1 page apart in pass 2, 17 pages apart in pass 1).  This is good in */
/* distributing the pages uniformly among the 4-way set-associative */
/* TLB cache. */

/* How does the above scheme work for the three pass case?  As you might */
/* imagine, more adjustments are necessary.  When doing a 1M FFT we will */
/* work in three passes.  Pass 2 looks at 0 up to 512 step 1, pass 1 */
/* looks at 0 up to 65536 step 256, and pass 0 looks at 0 up to 1048576 */
/* step 32768.  This corresponds to 5 levels in pass 0, 7 in pass 1, 8 in */
/* pass 2.  Notice above that both pass 1 and pass 2 look at the values */
/* 0 and 128 thus they should be on the same 4KB page.  Likewise, values */
/* 0 and 32768 are both used in pass 0 and pass 1 and should be on the */
/* same 4KB page.  After analyzing the various FFT sizes and TLB hit */
/* patterns, I settled on this memory layout: */

/* 1st 4KB page		*/
/* 0 128 256 384	*/
/* 512 640 768 896	*/
/* 16K 16K+128 ...	*/
/* 16K+512...		*/
/* ...			*/
/* 7*16K+512 ...	*/
/* 1 129 257 385	*/
/* 513 ...            	*/
/* ...			*/
/* 7*16K+512+7 ...	*/

/* This is much like the previous layout except that instead of 16 cache */
/* lines that are 512 apart, there are only 2 cache lines that are 512 apart */
/* and 8 that are 16K apart. */

/* To eliminate the cache line conflicts in pass 2, 32 bytes are wasted */
/* every 4KB.  To eliminate cache line conflicts in pass 1 and to keep */
/* the TLB hits uniform, after 16 pages we waste the rest of the 17th */
/* page and the first 8 cache lines of the 18th page.  To eliminate cache */
/* line conflicts in pass 0 and keep the TLB hits uniform, after 16 sets of */
/* 17 pages we move to the next 4KB page and waste the first 16 cache lines */
/* then after wasting 8 sets of 16, we waste another 32 bytes. */

/* It is now getting hard to visualize the FFT, so this program will print */
/* out the cache lines and TLB distributions for FFTS above 64K. */

#ifdef INCLUDED_PROGRAM
#include <stdio.h> 
unsigned long FFTLEN = 0;
/* Copy the addr function here */
void xmain (int incr, int endpt) { 
long	i, x, tlbs[16], lines[128]; 
for (i = 0; i <= 15; i++) tlbs[i] = 0; 
for (i = 0; i <= 127; i++) lines[i] = 0; 
printf ("\n\nTest fftlen: %d, incr: %d, endpt: %d\n", FFTLEN, incr, endpt);
for (i = 0; i < endpt; i += incr) { 
	x = (long) addr((long*)(32*19), i); 
	printf ("i: %d, addr: %d, page: %d, tlb line: %d, cache line: %d\n", 
		i, x, x >> 12, (x >> 12) & 15, (x >> 5) & 127); 
	tlbs[(x >> 12) & 15]++; lines[(x >> 5) & 127]++; 
} 
printf ("\n\nTLBS:"); for (i=0; i<=15; i++) printf (" %d", tlbs[i]); 
printf ("\n\nCache Lines:"); 
for (i = 0; i <= 127; i++) printf (" %d", lines[i]); 
printf ("\n"); 
} 
int main (int argc, char **argv) { 
FFTLEN = 65536 * 2;  xmain (1, 512); xmain (256, 32768); xmain (16384, FFTLEN);
FFTLEN = 65536 * 4;  xmain (1, 512); xmain (256, 65536); xmain (32768, FFTLEN);
FFTLEN = 65536 * 8;  xmain (1, 512); xmain (256, 32768); xmain (16384, FFTLEN);
FFTLEN = 65536 * 16; xmain (1, 512); xmain (256, 65536); xmain (32768, FFTLEN);
FFTLEN = 65536 * 32; xmain (1, 512); xmain (256, 32768); xmain (16384, FFTLEN);
FFTLEN = 65536 * 64; xmain (1, 512); xmain (256, 65536); xmain (32768, FFTLEN);
}
#endif

/* Below is a table of FFT sizes, FFT levels done in each of the three */
/* passes, L1 cache lines used, logical pages touched, and actual pages */
/* touched.  The logical and actual pages touched can be different because */
/* the waste bytes cause "spillage" of data from one 4KB page onto the next. */

/* FFT    FFT levels	L1 cache    Logical pages  Actual pages	*/
/* size   in each pass	lines used  accessed	   accessed	*/
/* ----   ------------	----------  -------------  ------------	*/
/* 4096K  7/7/8		128/128/128 32/16/16	   33/17/18	*/
/* 2048K  7/6/8		128/64/128  16/16/16	   17/17/18	*/
/* 1024K  5/7/8		32/128/128  8/16/16	   9/17/18	*/
/* 512K	  5/6/8		32/64/128   4/16/16	   5/17/18	*/
/* 256K	  3/7/8		8/128/128   2/16/16	   3/17/18	*/
/* 128K	  3/6/8		8/64/128    1/16/16	   2/17/18	*/
/* 64K	  8/8		128/128	    8/16	   9/18		*/
/* 32K	  7/8		64/128	    4/16	   5/18		*/
/* 16K	  6/8		32/128	    2/16	   3/18		*/
/* 8K	  7/6		128/32	    16/1 (flat memory model)	*/
/* 4K	  6/6		64/32	    8/1 (flat memory model)	*/
/* 2K	  5/6		32/32	    4/1 (flat memory model)	*/
/* 1K	  10		256	    2 (flat memory model)	*/
/* 512	  9		128	    1 (flat memory model)	*/
/* 256	  8		64	    1 (flat memory model)	*/

/* NOTE: I once had the brilliant idea of interleaving the sin/cos data */
/* with the FFT data.  That is, the data occupies the even cache */
/* lines and the sin/cos data is in the odd cache lines.  At first */
/* this seems counter productive, as only 4K of FFT data will now fit */
/* in the 8K L1 cache.  However, if you look at how an FFT operates */
/* you'll see loading FFT data, multiply by sin/cos data, store FFT data, */
/* load next block of FFT data, multiply by sin/cos data, store FFT data, */
/* etc.  By storing the sin/cos data in the odd cache lines, loading the */
/* next block of FFT data will toss out the previous block of FFT data */
/* rather than the reusable sin/cos data. */
/* For some reason, however, interleaving resulted in slower performance. */

unsigned long addr_offset (unsigned long fftlen, unsigned long i)
{
	unsigned long addr, i1, i2, i3, i6;

/* P4 uses a different memory layout - more suitable to SSE2 */

	if (CPU_FLAGS & CPU_SSE2) {
		unsigned long sets, pfa, temp;

/* Small FFTs use one pass, not very convoluted.  This the example for	*/
/* a length 2048 FFT:							*/
/*	0	512	1	513	1024	1536	1025	1537	*/
/*	2	...							*/
/*	...								*/
/*	510								*/
/* PFA-style FFTs are a little tricker.  See assembly code for example.	*/

		if (PASS2_LEVELS == 0) {
			sets = fftlen >> 3;
			if (i >= (fftlen >> 1)) {
				i6 = 1;
				i -= (fftlen >> 1);
			} else
				i6 = 0;
			i1 = i & 1; i >>= 1;
			i3 = 0;
			for (pfa = sets; pfa > 8; pfa >>= 1);
			if (pfa == 5) {
				temp = sets / 5;
				if (i < temp * 2) {
					sets = temp;
				} else {
					i3 = temp; i -= temp * 2;
					sets = temp * 4;
				}
			} else if (pfa == 7) {
				temp = sets / 7;
				if (i < temp * 2) {
					sets = temp;
				} else if (i < temp * 6) {
					i3 = temp; i -= temp * 2;
					sets = temp * 2;
				} else {
					i3 = temp * 3; i -= temp * 6;
					sets = temp * 4;
				}
			}
			i3 += i % sets; i /= sets;
			addr = (((((i3 << 1) + i6) << 1) + i1) << 1) + i;
			addr = addr * sizeof (double);
		}

/* Larger FFTs use two passes.  This the example for a length 64K FFT:	*/
/*	0	1K	16K	17K	32K	33K	48K	49K	*/
/*	1	...							*/
/*	...								*/
/*	1023	...							*/
/*	2K	...							*/
/*	...								*/

		else if (PASS2_LEVELS == 8) {
			sets = fftlen >> 10;
			if (i >= (fftlen >> 1)) {
				i6 = 1;
				i -= (fftlen >> 1);
			} else
				i6 = 0;
			i1 = i & 127; i >>= 7;
			i2 = i & 1; i >>= 1;
			i3 = 0;
			for (pfa = sets; pfa > 8; pfa >>= 1);
			if (pfa == 5) {
				temp = sets / 5;
				if (i < temp * 2) {
					sets = temp;
				} else {
					i3 = temp; i -= temp * 2;
					sets = temp * 4;
				}
			} else if (pfa == 7) {
				temp = sets / 7;
				if (i < temp * 2) {
					sets = temp;
				} else if (i < temp * 6) {
					i3 = temp; i -= temp * 2;
					sets = temp * 2;
				} else {
					i3 = temp * 3; i -= temp * 6;
					sets = temp * 4;
				}
			}
			i3 += i % sets; i /= sets;
			addr = (((((((i3 * 130) + i1) << 1) + i6) << 1) + i) << 1) + i2;
			addr = addr * sizeof (double);
		} else {
			sets = fftlen >> 13;
			if (i >= (fftlen >> 1)) {
				i6 = 1;
				i -= (fftlen >> 1);
			} else
				i6 = 0;
			i1 = i & 1023; i >>= 10;
			i2 = i & 1; i >>= 1;
			i3 = 0;
			for (pfa = sets; pfa > 8; pfa >>= 1);
			if (pfa == 5) {
				temp = sets / 5;
				if (i < temp * 2) {
					sets = temp;
				} else {
					i3 = temp; i -= temp * 2;
					sets = temp * 4;
				}
			} else if (pfa == 7) {
				temp = sets / 7;
				if (i < temp * 2) {
					sets = temp;
				} else if (i < temp * 6) {
					i3 = temp; i -= temp * 2;
					sets = temp * 2;
				} else {
					i3 = temp * 3; i -= temp * 6;
					sets = temp * 4;
				}
			}
			i3 += i % sets; i /= sets;
			if (SCRATCH_SIZE == 0)
				addr = i3 * 1090;
			else
				addr = i3 * (1024 + PASS1_CACHE_LINES*2);
			addr = ((((((addr + i1) << 1) + i6) << 1) + i) << 1) + i2;
			addr = addr * sizeof (double);
		}
	}

/* One pass x87 FFTs use a near flat memory model. */

	else if (PASS2_LEVELS == 0) {
		if (i >= (fftlen >> 1)) {
			i2 = 1;
			i -= (fftlen >> 1);
		} else
			i2 = 0;
		addr = i * 16 + i2 * 8;
	}

/* Two pass x87 FFTs use a near flat memory model.  Waste 64 bytes */
/* between 4KB.  Waste 64 bytes between every block (4KB, 16KB, or 64KB). */

	else {
		if (i >= (fftlen >> 1)) {
			i2 = 1;
			i -= (fftlen >> 1);
		} else
			i2 = 0;
		addr = i * 16 + i2 * 8 + (i >> 8) * 64 + (i >> PASS2_LEVELS) * 64;
	}

/* Return the offset */

	return (addr);
}

/* Return the address of ith element in the FFT array */

double *addr (gwnum g, unsigned long i)
{
	return ((double *) ((char *) g + addr_offset (FFTLEN, i)));
}

/* Return the size of a gwnum */

unsigned long gwnum_size (unsigned long fftlen)
{
	return (addr_offset (fftlen, fftlen - 1) + sizeof (double));
}

/* Each FFT word is multiplied by a two-to-phi value.  These */
/* routines set and get the FFT value without the two-to-phi */
/* multiplier. */

int get_fft_value (
	gwnum	g,
	unsigned long i,
	long	*retval)
{
	double	val;

/* Get the FFT data and validate it */

	val = * addr (g, i);
	if (! is_valid_double (val)) return (GWERROR_BAD_FFT_DATA);

/* Handle the rational FFT case quickly */

	if (RATIONAL_FFT) {
		*retval = (long) val;
	}

/* Multiply by two-to-minus-phi to generate an integer. */

	else {
		val = val * gwfft_weight_inverse_sloppy (i);
		if (val < -0.5)
			*retval = (long) (val - 0.5);
		else
			*retval = (long) (val + 0.5);
	}

/* Return success */

	return (0);
}

void set_fft_value (
	gwnum	g,
	unsigned long i,
	long	val)
{

/* Handle the rational FFT case quickly */

	if (RATIONAL_FFT || val == 0) {
		* addr (g, i) = val;
		return;
	}

/* Multiply by two-to-phi to generate the proper double. */

	* addr (g, i) = val * gwfft_weight_sloppy (i);
}

/* Some words in the FFT data contain floor(p/N), some words contain */
/* floor(p/N)+1 bits.  This function returns TRUE in the latter case. */

int is_big_word (
	unsigned long i)
{
	unsigned long base, next_base;

/* Compute the number of bits in this word.  It is a big word if */
/* the number of bits is more than BITS_PER_WORD. */

	base = gwfft_base (i);
	next_base = gwfft_base (i+1);
	return ((next_base - base) > BITS_PER_WORD);
}

/* Routine map a bit number into an FFT word and bit within that word */

void bitaddr (
	unsigned long bit,
	unsigned long *word,
	unsigned long *bit_in_word)
{

/* What word is the bit in? */

	*word = (unsigned long) ((double) bit / fft_bits_per_word);
	if (*word >= FFTLEN) *word = FFTLEN - 1;

/* Compute the bit within the word. */

	*bit_in_word = bit - gwfft_base (*word);
}

/* Return a description of the FFT type chosen */

void gwfft_description (
	char	*buf)		/* Buffer to return string in */
{
	sprintf (buf, "%sFFT length %lu%s",
		 ZERO_PADDED_FFT ? "zero-padded " :
		 GENERAL_MOD ? "generic reduction " : "",
		 FFTLEN >= 4096 ? FFTLEN / 1024 : FFTLEN,
		 FFTLEN >= 4096 ? "K" : "");
}

/* Return a string representation of a k/b/n/c combination */

void gw_as_string (
	char	*buf,		/* Buffer to return string in */
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c)		/* C in K*B^N+C */
{
	if (k != 1.0)
		sprintf (buf, "%.0f*%lu^%lu%c%lu", k, b, n,
			 c < 0 ? '-' : '+', abs (c));
	else if (b == 2 && c == -1)
		sprintf (buf, "M%lu", n);
	else
		sprintf (buf, "%lu^%lu%c%lu", b, n,
			 c < 0 ? '-' : '+', abs (c));
}

/* Return TRUE if we are operating near the limit of this FFT length */
/* Input argument is the percentage to consider as near the limit. */
/* For example, if percent is 1.0 and the FFT can handle 20 bits per word, */
/* then if there are more than 19.98 bits per word this function will */
/* return TRUE. */

int gwnear_fft_limit (
	double	pct)
{

/* Return TRUE if the virtual bits per word is near the maximum bits */
/* per word. */

	return (virtual_bits_per_word () >
			(100.0 - pct) / 100.0 * fft_max_bits_per_word);
}

/* Compute the virtual bits per word.  That is, the mersenne-mod-equivalent */
/* bits that this k,c combination uses.  For a non-zero-padded FFT */
/* log2(k) / 2 and log2(c) extra bits of precision are required.  This */
/* virtual value can tell us how close we are to this FFT length's limit. */

double virtual_bits_per_word ()
{
	double	logk, logc;

	if (ZERO_PADDED_FFT)
		return ((double) (PARG + PARG) / (double) FFTLEN);
	else {
		logk = log (KARG) / log ((double) 2.0);
		logc = log ((double) abs (CARG)) / log ((double) 2.0);
		return ((double) (logk + PARG) / (double) FFTLEN +
			0.5 * logk + 0.75 * logc);
	}
}

/* Given k,b,n,c determine the fft length */

unsigned long gwmap_to_fftlen (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. Must be two. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	unsigned long *info;

/* Handle exponents larger than SSE2 FFTs can handle */

	if (n > MAX_PRIME_SSE2) return (MAX_FFTLEN);

/* Get pointer to fft info and return the FFT length */

	gwinfo (k, b, n, c, 0);
	info = INFT[0];
	return (info[1]);
}

/* Given an fft length, determine the maximum allowable exponent */

unsigned long map_fftlen_to_max_exponent (
	unsigned long fftlen)
{
	unsigned long *info;

/* Get pointer to fft info and return the FFT length */

	gwinfo (1.0, 2, 0, -1, fftlen);
	info = INFT[0];
	return (info[0]);
}

/* Given an fft length, determine how much memory is used for */
/* normalization and sin/cos tables */

unsigned long gwmap_to_memused (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. Must be two. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	unsigned long *info;

/* Get pointer to fft info and return the memory used */

	gwinfo (k, b, n, c, 0);
	info = INFT[0];
	return (info[3]);
}

/* Make a guess as to how long a squaring will take. */

double gwmap_to_timing (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. Must be two. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c,		/* C in K*B^N+C. Must be rel. prime to K. */
	int	cpu_type)
{
	double	timing;
	unsigned long *info;

/* Get pointer to fft info */

	gwinfo (k, b, n, c, 0);
	info = INFT[0];

/* Use my PII-400 or P4-1400 timings as a guide. */

	timing = ((float *) info)[2];

/* Since the program is about 10% memory bound, the program will not */
/* speed up linearly with increase in chip speed.  Note, no attempt is */
/* made to differentiate between 66 MHz memory and 100 MHz memory - we're */
/* just returning an educated guess here. */

	if (CPU_FLAGS & CPU_SSE2) {
		timing = 0.10 * timing + 0.90 * timing * 1400.0 / CPU_SPEED;
	} else {
		timing = 0.10 * timing + 0.90 * timing * 400.0 / CPU_SPEED;
		if (cpu_type <= 4) timing *= REL_486_SPEED;
		if (cpu_type == 5) timing *= REL_PENT_SPEED;
		if (cpu_type == 7) timing *= REL_K6_SPEED;
		if (cpu_type == 11) timing *= REL_K7_SPEED;
		if (CPU_FLAGS & CPU_PREFETCH) timing *= 0.80;
	}
	return (timing);
}


/* Internal routine to help gwcopyzero */

void calc8ptrs (
	unsigned long n,
	unsigned long *ptrs)
{
	unsigned long i, j, k;

/* This is a grossly inefficient way to do this.  However, it should */
/* be called rarely. */

	for (i = 0; i < 8; i++) ptrs[i] = 0;
	for (i = 0; i < n; i++) {
		j = addr_offset (FFTLEN, i);
		k = (j & 63) >> 3;
		if (j >= ptrs[k]) ptrs[k] = j - (k << 3) + 64;
	}
}


/* Routine that sets up and calls assembly code to copy a gwnum from */
/* source to dest while zeroing some lower FFT words */

void gwcopyzero (
	gwnum	s,
	gwnum	d,
	unsigned long n)
{
static	unsigned long saved_n = 0;

	ASSERTG (((long *) s)[-7] == 0);	// Number not partially FFTed?

/* Handle case where no words are zeroed.  Some of the assembly routines */
/* do not like a word count of zero. */

	if (n == 0) {
		gwcopy (s, d);
		return;
	}

/* Call assembly language copy-and-zeroing routine */

	SRCARG = s;
	DESTARG = d;
	NUMARG = (long)(n);
	if ((CPU_FLAGS & CPU_SSE2) && (COPYZERO[0] == 0 || n != saved_n)) {
		saved_n = n;
		calc8ptrs (n, (unsigned long *) COPYZERO);
	}
	gw_copyzero ();

/* Copy the unnormalized add counter and clear the */
/* has been partially FFTed flag. */

	((long *) d)[-1] = ((long *) s)[-1];
	((long *) d)[-7] = 0;
}


/* Add a small constant at the specified bit position after the */
/* next multiplication.  This only works if k=1. */

void gwsetaddinatbit (
	long	value,
	unsigned long bit)
{
	unsigned long word, bit_in_word;

	ASSERTG (KARG == 1.0);

/* Tell assembly code to add the shifted value to the multiplication result. */

	bitaddr (bit, &word, &bit_in_word);
	raw_gwsetaddin (word, value << bit_in_word);
}

/* Routine that tells the assembly code to add a small value to the */
/* results of each multiply */

void gwsetaddin (
	long	value)
{
	unsigned long word, bit_in_word;

	ASSERTG (KARG == 1.0 || (BARG == 2 && abs (CARG) == 1));

/* In a zero-padded FFT, the value is added into ZPAD0 */

	if (ZERO_PADDED_FFT) {
		ADDIN_VALUE = (double) value;
		return;
	}

/* If value is even, shift it right and increment bit number.  This */
/* will ensure that we modify the proper FFT word. */

	for (bit_in_word = 0; value && (value & 1) == 0; value >>= 1)
		bit_in_word++;

/* Convert the input value to 1/k format.  Case 1 (2^n+/-1: Inverse of k */
/* is 1.  Case 2 (k*2^n-1): Inverse of k is 2^n.  Case 3 (k*2^n+1): Inverse */
/* of k is -2^n.  No other cases can be handled. */

	if (KARG == 1.0) {
		bitaddr (bit_in_word, &word, &bit_in_word);
	}
	else if (BARG == 2 && CARG == -1) {
		bitaddr (PARG + bit_in_word, &word, &bit_in_word);
	}
	else if (BARG == 2 && CARG == 1) {
		bitaddr (PARG + bit_in_word, &word, &bit_in_word);
		value = -value;
	}

/* Tell assembly code to add the shifted value to the multiplication result. */

	raw_gwsetaddin (word, value << bit_in_word);
}

/* Routine that tells the assembly code to add a small value to the */
/* results of each multiply */

void raw_gwsetaddin (
	unsigned long word,
	long	val)
{
	unsigned long row;

/* Compute the offset to the FFT data value */

	ADDIN_OFFSET = addr_offset (FFTLEN, word);

/* If this is a two-pass SSE2 FFT, then we need to tell the assembly code */
/* the affected "row", that is which set of pass 1 data the add-in will */
/* take place */

	if (CPU_FLAGS & CPU_SSE2) {
		if (PASS2_LEVELS == 0) {
			row = ADDIN_OFFSET & 31;
			if (row == 8) ADDIN_OFFSET += 8;
			if (row == 16) ADDIN_OFFSET -= 8;
		} else if (PASS2_LEVELS == 8) {
			row = (word & 127) / PASS1_CACHE_LINES;
			ADDIN_ROW = 128 / PASS1_CACHE_LINES - row;
			ADDIN_OFFSET -= row * PASS1_CACHE_LINES * 64;
		}

/* Factor in the blkdst value in xfft3.mac to compute the two pass */
/* SSE2 addin_offset. */

		else {
			row = (word & 1023) / PASS1_CACHE_LINES;
			ADDIN_ROW = (1024 / PASS1_CACHE_LINES - row) * 256;
			ADDIN_OFFSET -= row * PASS1_CACHE_LINES * 64;

/* This case is particularly nasty as we have to convert the FFT data offset */
/* into a scratch area offset.  In assembly language terms, this means */
/* subtracting out multiples of blkdst and adding in multiples of clmblkdst */
/* and clmblkdst8. */

			if (SCRATCH_SIZE) {
				row = ADDIN_OFFSET /
					(65536 + PASS1_CACHE_LINES * 128);
				ADDIN_OFFSET -= row *
					(65536 + PASS1_CACHE_LINES * 128);
				ADDIN_OFFSET +=
					row * (PASS1_CACHE_LINES * 64) +
					(row >> 3) * 128;
			}
		}
	}

/* And now x87 FFTs also can use a scratch area.  Like the SSE2 code */
/* we have to convert the FFT data offsets for two-pass FFTs. */

	if (! (CPU_FLAGS & CPU_SSE2) && PASS2_LEVELS) {
		unsigned long num_cache_lines, cache_line;

		num_cache_lines = (1 << (PASS2_LEVELS - 1));
		cache_line = ((word >> 1) & (num_cache_lines - 1));

		ADDIN_ROW = ((num_cache_lines>>7) - (cache_line>>7)) * 65536 +
			    (128 / PASS1_CACHE_LINES -
			     (cache_line & 127) / PASS1_CACHE_LINES) * 256;
		ADDIN_OFFSET -= (cache_line >> 7) * 64 +
				(cache_line / PASS1_CACHE_LINES) *
				PASS1_CACHE_LINES * 32;

/* This case is particularly nasty as we have to convert the FFT data offset */
/* into a scratch area offset.  In assembly language terms, this means */
/* subtracting out multiples of blkdst and adding in multiples of clmblkdst */
/* and clmblkdst32. */

		if (SCRATCH_SIZE) {
			unsigned long blkdst;

			blkdst = addr_offset (FFTLEN, 1 << PASS2_LEVELS);
			row = ADDIN_OFFSET / blkdst;
			ADDIN_OFFSET -= row * blkdst;
			ADDIN_OFFSET += row * (PASS1_CACHE_LINES * 32);

/* Handle the FFTs where clmblkdst32 is used */

			if ((FFTLEN >> (PASS2_LEVELS+1)) >= 128)
				ADDIN_OFFSET += (row >> 5) * 64;
		}
	}

/* Set the addin value - multiply it by two-to-phi and FFTLEN/2/KARG. */

	ADDIN_VALUE = (double) val * gwfft_weight_sloppy (word) * FFTLEN * 0.5 / KARG;
}


/* Routine to add a small number (-255 to 255) to a gwnum.  Some day, */
/* I might optimize this routine for the cases where just one or two */
/* doubles need to be modified in the gwnum */

void gwaddsmall (
	gwnum	g,		/* Gwnum to add a value into */
	int	addin)		/* Small value to add to g */
{
	gwnum	tmp;

/* A simple brute-force implementation */

	tmp = gwalloc ();
	if (addin >= 0) {
		dbltogw ((double) addin, tmp);
		gwaddquick (tmp, g);
	} else {
		dbltogw ((double) -addin, tmp);
		gwsubquick (tmp, g);
	}
	gwfree (tmp);
}

/********************************************************/
/* Routines to convert between gwnums and other formats */
/********************************************************/

void specialmodg (giant	g);

/* Convert a double to a gwnum */

void dbltogw (double d, gwnum g)
{
	giantstruct tmp;
	unsigned long tmparray[2];

	tmp.n = (unsigned long *) &tmparray;
	setmaxsize (&tmp, 2);
	dbltog (d, &tmp);
	gianttogw (&tmp, g);
}

/* Convert a binary value to a gwnum */

void binarytogw (
	unsigned long *array,	/* Array containing the binary value */
	unsigned long arraylen,	/* Length of the array */
	gwnum	n)		/* Destination gwnum */
{
	giantstruct tmp;
	tmp.sign = arraylen;
	tmp.n = array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	gianttogw (&tmp, n);
}

/* Convert a giant to gwnum FFT format.  Giant must be a positive number. */

void gianttogw (
	giant	a,
	gwnum	g)
{
	giant	newg;
	unsigned long i, mask1, mask2, e1len;
	int	bits1, bits2, bits_in_next_binval;
	unsigned long *e1, binval, carry;

/* To make the mod k*b^n+c step faster, gwnum's are pre-multiplied by 1/k */
/* If k is greater than 1, then we calculate the inverse of k, multiply */
/* the giant by the inverse of k, and do a mod k*b^n+c. */

	if (KARG > 1) {
		newg = popg ((((unsigned long) bit_length >> 5) + 1) * 2);

		/* Easy case 1 (k*2^n-1): Inverse of k is 2^n */

		if (BARG == 2 && CARG == -1) {
			gtog (a, newg);
			gshiftleft (PARG, newg);
		}

		/* Easy case 2 (k*2^n+1): Inverse of k is -2^n */

		else if (BARG == 2 && CARG == 1) {
			gtog (a, newg);
			negg (newg);
			gshiftleft (PARG, newg);
		}

		else {				/* General inverse case */
			giant	n;
			n = popg (((unsigned long) bit_length >> 5) + 1);
			ultog (BARG, n);	/* Compute k*b^n+c */
			power (n, PARG);
			dblmulg (KARG, n);
			iaddg (CARG, n);
			dbltog (KARG, newg);	/* Compute 1/k */
			invg (n, newg);
			ASSERTG (newg->sign > 0);  /* Assert inverse found */
			mulg (a, newg);		/* Multiply input num by 1/k */
			pushg (1);
		}

		specialmodg (newg);
		a = newg;
	}

/* Now convert the giant to FFT format */

	ASSERTG (a->sign >= 0);
	e1len = a->sign;
	e1 = a->n;

	bits1 = BITS_PER_WORD;
	bits2 = bits1 + 1;
	mask1 = (1L << bits1) - 1;
	mask2 = (1L << bits2) - 1;
	if (e1len) {binval = *e1++; e1len--; bits_in_next_binval = 32;}
	else binval = 0;
	carry = 0;
	for (i = 0; i < FFTLEN; i++) {
		int	big_word, bits;
		long	value, mask;
		big_word = is_big_word (i);
		bits = big_word ? bits2 : bits1;
		mask = big_word ? mask2 : mask1;
		if (i == FFTLEN - 1) value = binval;
		else value = binval & mask;
		value = value + carry;
		if (value > (mask >> 1) && bits > 1 && i != FFTLEN - 1) {
			value = value - (mask + 1);
			carry = 1;
		} else {
			carry = 0;
		}
		set_fft_value (g, i, value);

		binval >>= bits;
		if (e1len == 0) continue;
		if (bits_in_next_binval < bits) {
			if (bits_in_next_binval)
				binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
			bits -= bits_in_next_binval;
			e1++; e1len--; bits_in_next_binval = 32;
			if (e1len == 0) continue;
		}
		if (bits) {
			binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
			bits_in_next_binval -= bits;
		}
	}
	((long *) g)[-1] = 1;	/* Set unnormalized add counter */
	((long *) g)[-7] = 0;	/* Clear has been partially FFTed flag */

/* Free allocated memory */

	if (KARG > 1.0) pushg (1);
}

/* Convert a gwnum to a binary value.  Returns the number of longs */
/* written to the array.  The array is NOT zero-padded.  Returns a */
/* negative number if an error occurs during the conversion.  An error */
/* can happen if the FFT data contains a NaN or infinity value. */

long gwtobinary (
	gwnum	n,		/* Source gwnum */
	unsigned long *array,	/* Array to contain the binary value */
	unsigned long arraylen)	/* Maximum size of the array */
{
	giant	tmp;
	int	err_code;

	tmp = popg ((PARG >> 5) + 5);
	err_code = gwtogiant (n, tmp);
	if (err_code < 0) return (err_code);
	ASSERTG ((unsigned long) tmp->sign >= arraylen);
	memcpy (array, tmp->n, tmp->sign * sizeof (unsigned long));
	pushg (1);
	return (tmp->sign);
}

/* Convert a gwnum to a giant.  WARNING: Caller must allocate an array that */
/* is several words larger than the maximum result that can be returned. */
/* This is a gross kludge that lets gwtogiant use the giant for intermediate */
/* calculations. */

int gwtogiant (
	gwnum	gg,
	giant	v)
{
	long	val;
	int	j, bits, bitsout, carry, err_code;
	unsigned long i, limit, *outptr;

	ASSERTG (((long *) gg)[-7] == 0);	// Number not partially FFTed?

/* If this is a general-purpose mod, then only convert the needed words */
/* which will be less than half the FFT length.  If this is a zero padded */
/* FFT, then only convert a little more than half of the FFT data words. */
/* For a DWT, convert all the FFT data. */

	if (GENERAL_MOD) limit = GW_ZEROWORDSLOW + 3;
	else if (ZERO_PADDED_FFT) limit = FFTLEN / 2 + 4;
	else limit = FFTLEN;

/* Collect bits until we have all of them */

	carry = 0;
	bitsout = 0;
	outptr = v->n;
	*outptr = 0;
	for (i = 0; i < limit; i++) {
		err_code = get_fft_value (gg, i, &val);
		if (err_code) return (err_code);
		bits = BITS_PER_WORD;
		if (is_big_word (i)) bits++;
		val += carry;
		for (j = 0; j < bits; j++) {
			*outptr >>= 1;
			if (val & 1) *outptr += 0x80000000;
			val >>= 1;
			bitsout++;
			if (bitsout == 32) {
				outptr++;
				bitsout = 0;
			}
		}
		carry = val;
	}

/* Finish outputting the last word and any carry data */

	while (bitsout || (carry != -1 && carry != 0)) {
		*outptr >>= 1;
		if (carry & 1) *outptr += 0x80000000;
		carry >>= 1;
		bitsout++;
		if (bitsout == 32) {
			outptr++;
			bitsout = 0;
		}
	}

/* Set the length */

	v->sign = (outptr - v->n);
	while (v->sign && v->n[v->sign-1] == 0) v->sign--;

/* If carry is -1, the gwnum is negative.  Ugh.  Flip the bits and sign. */

	if (carry == -1) {
		for (j = 0; j < v->sign; j++) v->n[j] = ~v->n[j];
		while (v->sign && v->n[v->sign-1] == 0) v->sign--;
		iaddg (1, v);
		v->sign = -v->sign;
	}

/* The gwnum is not guaranteed to be smaller than k*b^n+c.  Handle this */
/* possibility.  This also converts negative values to positive. */

	specialmodg (v);

/* Since all gwnums are premultiplied by the inverse of k, we must now */
/* multiply by k to get the true result. */

	if (KARG > 1) {
		giant	newg;
		newg = popg (((unsigned long) bit_length >> 5) + 3);
		dbltog (KARG, newg);
		mulg (v, newg);
		specialmodg (newg);
		gtog (newg, v);
		pushg (1);
	}

/* Return success */

	return (0);
}

/* Special modg.  This is a fast implementation of mod k*2^n+c using just */
/* shifts, adds, and divide and mul by small numbers.  All others moduli */
/* call the slow giants code. */

void specialmodg (
	giant	g)
{
	int	neg, count;
	giant	n;

/* If the modulus is a general-purpose number, then let the giants code */
/* do the work. */

	if (GENERAL_MOD) {
		modg (GW_MODULUS, g);
		return;
	}

/* Calculate the modulo number - k*b^n+c */

	n = popg (((unsigned long) bit_length >> 5) + 1);
	ultog (BARG, n);
	power (n, PARG);
	dblmulg (KARG, n);
	iaddg (CARG, n);

/* If b is not 2 let the giants code do the work. */

	if (BARG != 2) {
		modg (n, g);
		pushg (1);
		return;
	}

/* Do the quick modulus code twice because in the case where */
/* abs(c) > k once won't get us close enough. */

	neg = FALSE;
	for (count = 0; count < 2; count++) {

/* Handle negative input values */

	    neg ^= (g->sign < 0);
	    g->sign = abs (g->sign);

/* If number is bigger than the modulus, do a mod using shifts and adds */
/* This will get us close to the right answer. */

	    if (gcompg (g, n) > 0) {
		giant	tmp1, tmp2, tmp3;

/* Allocate temporaries */

		tmp1 = popg (((unsigned long) bit_length >> 5) + 5);
		tmp2 = popg ((((unsigned long) bit_length >> 5) + 5) * 2);
		tmp3 = popg (((unsigned long) bit_length >> 5) + 5);

/* Calculate the modulo by dividing the upper bits of k, multiplying by */
/* c and subtracting that from the bottom bits. */
		
		gtogshiftright (PARG, g, tmp1);	// Upper bits
		gmaskbits (PARG, g);		// Lower bits

		gtog (tmp1, tmp2);
		dbltog (KARG, tmp3);
		divg (tmp3, tmp1);		// Upper bits over K
		mulg (tmp1, tmp3);
		subg (tmp3, tmp2);		// Upper bits mod K

		gshiftleft (PARG, tmp2);
		addg (tmp2, g);			// Upper bits mod K+lower bits

		imulg (CARG, tmp1);		// Upper bits over K times C
		subg (tmp1, g);

		pushg (3);
	    }
	}

/* Add or subtract n until the g is between 0 and n-1 */

	while (g->sign < 0) addg (n, g);
	while (gcompg (g, n) >= 0) subg (n, g);

/* If input was negative, return k*b^n+c - g */

	if (neg && g->sign) {
		g->sign = -g->sign;
		addg (n, g);
	}

/* Free memory */

	pushg (1);
}

/******************************************************************/
/* Wrapper routines for the multiplication assembly code routines */
/******************************************************************/

/* Internal wrapper routine to call fftmul assembly code */

void raw_gwfftmul (
	gwnum	s,
	gwnum	d)
{
	unsigned long norm_count1, norm_count2;
	double	sumdiff;

/* Get the unnormalized add count for later use */

	norm_count1 = ((unsigned long *) s)[-1];
	norm_count2 = ((unsigned long *) d)[-1];

/* Call the assembly code */

	SRCARG = s;
	DESTARG = d;
	gw_mul ();
	if (! is_valid_double (gwsumout (d))) GWERROR |= 1;
	fftinc (2);

/* Adjust if necessary the SUM(INPUTS) vs. SUM(OUTPUTS).  If norm_count */
/* is more than one, then the sums will be larger than normal.  This */
/* could trigger a spurious MAXDIFF warning.  Shrink the two SUMS to */
/* compensate. */

	if (norm_count1 != 1 || norm_count2 != 1) {
		double	adjustment;
		adjustment = 1.0 / ((double)norm_count1 * (double)norm_count2);
		gwsuminp (d) *= adjustment;
		gwsumout (d) *= adjustment;
	}

/* Test SUM(INPUTS) vs. SUM(OUTPUTS) */

	sumdiff = gwsuminp (d) - gwsumout (d);
	if (fabs (sumdiff) > MAXDIFF) GWERROR |= 2; 

/* Reset the unnormalized add count */

	((unsigned long *) d)[-1] = 1;
}

/* Common code to emulate the modulo with two multiplies in the */
/* general purpose case */

void emulate_mod (
	gwnum	s)		/* Source and destination */
{
	gwnum	tmp;
	double	saved_addin_value;

/* Save and clear the addin value */

	saved_addin_value = ADDIN_VALUE;
	ADDIN_VALUE = 0.0;

/* Copy the number and zero out the low words. */

	tmp = gwalloc ();
	gwcopyzero (s, tmp, GW_ZEROWORDSLOW);

/* Multiply by the reciprocal that has been carefully shifted so that the */
/* integer part of the result wraps to the lower FFT words.  Adjust the */
/* normalization routine so that the FFT code zeroes the high FFT words */
/* and we are left with just the quotient! */

	NORMRTN = GWPROCPTRS[norm_routines + 4 + (NORMNUM & 1)];
	raw_gwfftmul (GW_RECIP, tmp);

/* Muliply quotient and modulus.  Select normalization routine that does */
/* not zero the high FFT words. */

	NORMRTN = GWPROCPTRS[norm_routines + (NORMNUM & 1)];
	raw_gwfftmul (GW_MODULUS_FFT, tmp);

/* Subtract from the original number to get the remainder */

	gwsub (tmp, s);
	gwfree (tmp);

/* Restore the addin value */

	ADDIN_VALUE = saved_addin_value;
}

/* User-visible routines */

void gwfft (			/* Forward FFT */
	gwnum	s,		/* Source number */
	gwnum	d)		/* Destination (can overlap source) */
{

/* Copy the unnormalized add count */

	((unsigned long *) d)[-1] = ((unsigned long *) s)[-1];

/* Call the assembly code */

	SRCARG = s;
	DESTARG = d;
	gw_fft ();
	fftinc (1);
}

void gwsquare (			/* Square a number */
	gwnum	s)		/* Source and destination */
{
	unsigned long norm_count;
	double	sumdiff;

/* Get the unnormalized add count for later use */

	norm_count = ((unsigned long *) s)[-1];

/* Call the assembly code */

	NORMRTN = GWPROCPTRS[norm_routines + NORMNUM];
	DESTARG = s;
	gw_square ();
	if (! is_valid_double (gwsumout (s))) GWERROR |= 1;
	fftinc (2);

/* Adjust if necessary the SUM(INPUTS) vs. SUM(OUTPUTS).  If norm_count */
/* is more than one, then the sums will be larger than normal.  This */
/* could trigger a spurious MAXDIFF warning.  Shrink the two SUMS to */
/* compensate. */

	if (norm_count != 1) {
		double	adjustment;
		adjustment = 1.0 / ((double) norm_count * (double) norm_count);
		gwsuminp (s) *= adjustment;
		gwsumout (s) *= adjustment;
	}

/* Test SUM(INPUTS) vs. SUM(OUTPUTS) */

	sumdiff = gwsuminp (s) - gwsumout (s);
	if (fabs (sumdiff) > MAXDIFF) GWERROR |= 2; 

/* Reset the unnormalized add count */

	((unsigned long *) s)[-1] = 1;

/* Emulate mod with 2 multiplies case */

	if (GENERAL_MOD) emulate_mod (s);
}

void gwfftmul (			/* Multiply already FFTed source with dest */
	gwnum	s,		/* Already FFTed source number */
	gwnum	d)		/* Non-FFTed source. Also destination */
{

/* Call the assembly code */

	NORMRTN = GWPROCPTRS[norm_routines + NORMNUM];
	raw_gwfftmul (s, d);

/* Emulate mod with 2 multiplies case */

	if (GENERAL_MOD) emulate_mod (d);
}

void gwfftfftmul (		/* Multiply two already FFTed sources */
	gwnum	s,		/* Already FFTed source number */
	gwnum	s2,		/* Already FFTed source number */
	gwnum	d)		/* Destination (can overlap sources) */
{
	unsigned long norm_count1, norm_count2;
	double	sumdiff;

/* Get the unnormalized add count for later use */

	norm_count1 = ((unsigned long *) s)[-1];
	norm_count2 = ((unsigned long *) s2)[-1];

/* Call the assembly code */

	NORMRTN = GWPROCPTRS[norm_routines + NORMNUM];
	SRCARG = s;
	SRC2ARG = s2;
	DESTARG = d;
	gw_mulf ();
	if (! is_valid_double (gwsumout (d))) GWERROR |= 1;
	fftinc (1);

/* Adjust if necessary the SUM(INPUTS) vs. SUM(OUTPUTS).  If norm_count */
/* is more than one, then the sums will be larger than normal.  This */
/* could trigger a spurious MAXDIFF warning.  Shrink the two SUMS to */
/* compensate. */

	if (norm_count1 != 1 || norm_count2 != 1) {
		double	adjustment;
		adjustment = 1.0 / ((double)norm_count1 * (double)norm_count2);
		gwsuminp (d) *= adjustment;
		gwsumout (d) *= adjustment;
	}

/* Test SUM(INPUTS) vs. SUM(OUTPUTS) */

	sumdiff = gwsuminp (d) - gwsumout (d);
	if (fabs (sumdiff) > MAXDIFF) GWERROR |= 2; 

/* Reset the unnormalized add count */

	((unsigned long *) d)[-1] = 1;

/* Emulate mod with 2 multiplies case */

	if (GENERAL_MOD) emulate_mod (d);
}

void gwmul (			/* Multiply source with dest */
	gwnum	s,		/* Source number (changed to FFTed source!) */
	gwnum	d)		/* Source and destination */
{
	gwfft (s,s);
	gwfftmul (s,d);
}

void gwsafemul (		/* Multiply source with dest */
	gwnum	s,		/* Source number (not changed) */
	gwnum	d)		/* Source and destination */
{
	gwnum	qqq;

	qqq = gwalloc ();
	gwfft (s, qqq);
	gwfftmul (qqq, d);
	gwfree (qqq);
}

/* Generate random FFT data */

void gw_random_number (
	gwnum	x)
{
	giant	g;
	unsigned long i, len;

/* Generate the random number */

	len = (((unsigned long) bit_length) >> 5) + 1;
	g = popg (len);
	for (i = 0; i < len; i++) {
		g->n[i] = ((unsigned long) rand() << 20) +
			  ((unsigned long) rand() << 10) +
			  (unsigned long) rand();
	}
	g->sign = len;
	specialmodg (g);
	gianttogw (g, x);
	pushg (1);
}

/* Square a number using a slower method that will have reduced */
/* round-off error on non-random input data.  Caller must make sure the */
/* input number has not been partially or fully FFTed. */

void gwsquare_carefully (
	gwnum	s)		/* Source and destination */
{
	gwnum	tmp1, tmp2;
	double	saved_addin_value;

/* Generate a random number, if we have't already done so */

	if (GW_RANDOM == NULL) {
		GW_RANDOM = gwalloc ();
		gw_random_number (GW_RANDOM);
	}

/* Save and clear the addin value */

	saved_addin_value = ADDIN_VALUE;
	ADDIN_VALUE = 0.0;

/* Now do the squaring using three multiplies and adds */

	tmp1 = gwalloc ();
	tmp2 = gwalloc ();
	gwstartnextfft (0);			/* Disable POSTFFT */
	gwadd3 (s, GW_RANDOM, tmp1);		/* Compute s+random */
	gwfft (GW_RANDOM, tmp2);
	gwfftmul (tmp2, s);			/* Compute s*random */
	gwfftfftmul (tmp2, tmp2, tmp2);		/* Compute random^2 */
	ADDIN_VALUE = saved_addin_value;	/* Restore the addin value */
	gwsquare (tmp1);			/* Compute (s+random)^2 */
	gwsubquick (tmp2, tmp1);		/* Calc s^2 from 3 results */
	gwaddquick (s, s);
	gwsub3 (tmp1, s, s);

/* Free memory and return */

	gwfree (tmp1);
	gwfree (tmp2);
}

/*********************************************************/
/* Wrapper routines for the add and sub assembly code    */
/*********************************************************/

void gwadd3quick (		/* Add two numbers without normalizing */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{

	ASSERTG (((unsigned long *) s1)[-7] == ((unsigned long *) s2)[-7]);

/* Update the count of unnormalized adds and subtracts */

	((unsigned long *) d)[-1] =
		((unsigned long *) s1)[-1] + ((unsigned long *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((unsigned long *) d)[-7] = ((unsigned long *) s1)[-7];

/* Now do the add */

	SRCARG = s1;
	SRC2ARG = s2;
	DESTARG = d;
	gw_addq ();
}

void gwsub3quick (		/* Compute s1 - s2 without normalizing */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{

	ASSERTG (((unsigned long *) s1)[-7] == ((unsigned long *) s2)[-7]);

/* Update the count of unnormalized adds and subtracts */

	((unsigned long *) d)[-1] =
		((unsigned long *) s1)[-1] + ((unsigned long *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((unsigned long *) d)[-7] = ((unsigned long *) s1)[-7];

/* Now do the subtract */

	SRCARG = s2;
	SRC2ARG = s1;
	DESTARG = d;
	gw_subq ();
}

void gwaddsub4quick (		/* Add & sub two numbers without normalizing */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2)		/* Destination #2 */
{

	ASSERTG (((unsigned long *) s1)[-7] == ((unsigned long *) s2)[-7]);

/* Update the counts of unnormalized adds and subtracts */

	((unsigned long *) d1)[-1] =
	((unsigned long *) d2)[-1] =
		((unsigned long *) s1)[-1] + ((unsigned long *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((unsigned long *) d1)[-7] =
	((unsigned long *) d2)[-7] = ((unsigned long *) s1)[-7];

/* Now do the add & subtract */

	SRCARG = s1;
	SRC2ARG = s2;
	DESTARG = d1;
	DEST2ARG = d2;
	gw_addsubq ();
}


void gwadd3 (			/* Add two numbers normalizing if needed */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	unsigned long normcnt1, normcnt2;

	ASSERTG (((unsigned long *) s1)[-7] == 0);
	ASSERTG (((unsigned long *) s2)[-7] == 0);

/* Get counts of unnormalized adds and subtracts */

	normcnt1 = ((unsigned long *) s1)[-1];
	normcnt2 = ((unsigned long *) s2)[-1];

/* Set the has-been-partially-FFTed flag */

	((unsigned long *) d)[-7] = 0;

/* Now do the add */

	SRCARG = s1;
	SRC2ARG = s2;
	DESTARG = d;
	if (normcnt1 + normcnt2 <= EXTRA_BITS) {
		gw_addq ();
		((unsigned long *) d)[-1] = normcnt1 + normcnt2;
	} else {
		gw_add ();
		((unsigned long *) d)[-1] = 1;
	}
}

void gwsub3 (			/* Compute s1 - s2 normalizing if needed */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	unsigned long normcnt1, normcnt2;

	ASSERTG (((unsigned long *) s1)[-7] == 0);
	ASSERTG (((unsigned long *) s2)[-7] == 0);

/* Get counts of unnormalized adds and subtracts */

	normcnt1 = ((unsigned long *) s1)[-1];
	normcnt2 = ((unsigned long *) s2)[-1];

/* Set the has-been-partially-FFTed flag */

	((unsigned long *) d)[-7] = 0;

/* Now do the subtract */

	SRCARG = s2;
	SRC2ARG = s1;
	DESTARG = d;
	if (normcnt1 + normcnt2 <= EXTRA_BITS) {
		gw_subq ();
		((unsigned long *) d)[-1] = normcnt1 + normcnt2;
	} else {
		gw_sub ();
		((unsigned long *) d)[-1] = 1;
	}
}

void gwaddsub4 (		/* Add & sub two nums normalizing if needed */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2)		/* Destination #2 */
{
	unsigned long normcnt1, normcnt2;

	ASSERTG (((unsigned long *) s1)[-7] == 0);
	ASSERTG (((unsigned long *) s2)[-7] == 0);

/* Get counts of unnormalized adds and subtracts */

	normcnt1 = ((unsigned long *) s1)[-1];
	normcnt2 = ((unsigned long *) s2)[-1];

/* Set the has-been-partially-FFTed flag */

	((unsigned long *) d1)[-7] = ((unsigned long *) d2)[-7] = 0;

/* Now do the add & subtract */

	SRCARG = s1;
	SRC2ARG = s2;
	DESTARG = d1;
	DEST2ARG = d2;
	if (normcnt1 + normcnt2 <= EXTRA_BITS) {
		gw_addsubq ();
		((unsigned long *) d1)[-1] =
		((unsigned long *) d2)[-1] = normcnt1 + normcnt2;
	} else {
		gw_addsub ();
		((unsigned long *) d1)[-1] =
		((unsigned long *) d2)[-1] = 1;
	}
}


void gwfftadd3 (		/* Add two FFTed numbers */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{

	ASSERTG (((unsigned long *) s1)[-7] == ((unsigned long *) s2)[-7]);

/* Update the count of unnormalized adds and subtracts */

	((unsigned long *) d)[-1] =
		((unsigned long *) s1)[-1] + ((unsigned long *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((unsigned long *) d)[-7] = ((unsigned long *) s1)[-7];

/* If this is a zero-padded FFT, then also add the 7 copied doubles in */
/* the gwnum header */

	if (ZERO_PADDED_FFT) {
		d[-5] = s1[-5] + s2[-5];
		d[-6] = s1[-6] + s2[-6];
		d[-7] = s1[-7] + s2[-7];
		d[-8] = s1[-8] + s2[-8];
		d[-9] = s1[-9] + s2[-9];
		d[-10] = s1[-10] + s2[-10];
		d[-11] = s1[-11] + s2[-11];
	}

/* Now do the add */

	SRCARG = s1;
	SRC2ARG = s2;
	DESTARG = d;
	gw_addf ();
}

void gwfftsub3 (		/* Compute FFTed s1 - FFTed s2 */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{

	ASSERTG (((unsigned long *) s1)[-7] == ((unsigned long *) s2)[-7]);

/* Update the count of unnormalized adds and subtracts */

	((unsigned long *) d)[-1] =
		((unsigned long *) s1)[-1] + ((unsigned long *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((unsigned long *) d)[-7] = ((unsigned long *) s1)[-7];

/* If this is a zero-padded FFT, then also subtract the 7 copied doubles in */
/* the gwnum header */

	if (ZERO_PADDED_FFT) {
		d[-5] = s1[-5] - s2[-5];
		d[-6] = s1[-6] - s2[-6];
		d[-7] = s1[-7] - s2[-7];
		d[-8] = s1[-8] - s2[-8];
		d[-9] = s1[-9] - s2[-9];
		d[-10] = s1[-10] - s2[-10];
		d[-11] = s1[-11] - s2[-11];
	}

/* Now do the subtract */

	SRCARG = s2;
	SRC2ARG = s1;
	DESTARG = d;
	gw_subf ();
}

void gwfftaddsub4 (		/* Add & sub two FFTed numbers */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2)		/* Destination #2 */
{

	ASSERTG (((unsigned long *) s1)[-7] == ((unsigned long *) s2)[-7]);

/* Update the counts of unnormalized adds and subtracts */

	((unsigned long *) d1)[-1] =
	((unsigned long *) d2)[-1] =
		((unsigned long *) s1)[-1] + ((unsigned long *) s2)[-1];

/* Copy the has-been-partially-FFTed flag */

	((unsigned long *) d1)[-7] =
	((unsigned long *) d2)[-7] = ((unsigned long *) s1)[-7];

/* If this is a zero-padded FFT, then also add & sub the 7 copied doubles in */
/* the gwnum header */

	if (ZERO_PADDED_FFT) {
		d1[-5] = s1[-5] + s2[-5];	d2[-5] = s1[-5] - s2[-5];
		d1[-6] = s1[-6] + s2[-6];	d2[-6] = s1[-6] - s2[-6];
		d1[-7] = s1[-7] + s2[-7];	d2[-7] = s1[-7] - s2[-7];
		d1[-8] = s1[-8] + s2[-8];	d2[-8] = s1[-8] - s2[-8];
		d1[-9] = s1[-9] + s2[-9];	d2[-9] = s1[-9] - s2[-9];
		d1[-10] = s1[-10] + s2[-10];	d2[-10] = s1[-10] - s2[-10];
		d1[-11] = s1[-11] + s2[-11];	d2[-11] = s1[-11] - s2[-11];
	}

/* Now do the add & subtract */

	SRCARG = s1;
	SRC2ARG = s2;
	DESTARG = d1;
	DEST2ARG = d2;
	gw_addsubf ();
}
