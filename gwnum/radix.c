/*----------------------------------------------------------------------
| radix.c
|
| This file contains the C routines for radix conversion when required
| by gianttogw or gwtogiant.
| 
|  Copyright 2020 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "gwnum.h"
#include "gwutil.h"
#include "radix.h"

// IDEAS:  Use libdivide library to increase BRUTE_WORDS to 4

/* Forward declarations */

void radix_gwcopy (gwhandle *src_gwdata, gwhandle *dst_gwdata, gwnum src, gwnum dst, uint64_t big_word_flags, int num_big_word_flags);
void brute_convert (gwhandle *gwdata, giant g, gwnum x, int offset, uint64_t big_word_flags);
__inline int intgcd(int k, int c) {
	stackgiant(kg,2); stackgiant(cg,2);itog(k,kg);itog(c,cg);gcdg(kg,cg);
	return (cg->n[0]);
}

/* Internal routine to convert a giant to gwnum FFT format. */

void nonbase2_gianttogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	g,		/* Input giant */
	gwnum	x)		/* Output gwnum */
{
	gwhandle new_gwdata;		/* A smaller, k=1 FFT to speed our radix conversion */
	gwnum	pow2_multiplier, t1, t3;
	int	i, num_chunks, fft_words_per_mult, num_big_word_flags;
	uint64_t big_word_flags;

#define BRUTE_FORCE_WORDS		2
#define giant_extract(s,off,len,d)	{if ((s)->sign <= (off)) (d)->sign = 0; \
					 else { (d)->sign = min ((s)->sign - (off), (len)); \
						memcpy ((d)->n, (s)->n+(off), (d)->sign * sizeof (uint32_t)); \
						while((d)->sign && (d)->n[(d)->sign-1] == 0) (d)->sign--;}}

	ASSERTG (g->sign >= 0);		/* We only handle positive numbers */

// Brute force conversion for small inputs

	if (g->sign <= BRUTE_FORCE_WORDS) {
		stackgiant(tmpg,BRUTE_FORCE_WORDS);
		gtog (g, tmpg);
		dbltogw (gwdata, 0.0, x);
		for (i = 0, big_word_flags = 0; i < 64; i++) if (is_big_word (gwdata, i)) big_word_flags |= ((uint64_t) 1) << i;
		brute_convert (gwdata, tmpg, x, 0, big_word_flags);
		return;
	}

// We do a new gwsetup to do our multiplications here.  In some cases, we could use the existing gwdata, but if k != 1 the
// caller's gwdata normalization code does unwanted multiplications by k.  In addition, handling the smallword / bigword differences
// would be tedious.

// The smallest FFT we can get away must hold all our initial partial results.  This FFT size likely won't work because binary
// converted data will not fit exactly in radix-b FFT words.  Try larger and larger FFT sizes until we find an FFT that suits our needs.
// We don't restrict ourselves to rational FFTs to avoid bigword/smallword complexities, instead we require a less restrictive
// rule that every partial multiplication result uses the same bigword/smallword sequence.

	if (1) {		// Later this will test for a cached new_gwdata
		int	exp, fftlen;			/* Exponent and fft length used to setup new_gwdata */
		int	num_pairs;
		int	b_per_mult;			/* Number of radix b that must appear in each half-pair multiply result */

		// Calc how big each partial result is in radix-b (with no borrow from next result)
		// This is the number of bits in a multiplication result plus one to avoid borrow divided by bits in radix b.
		b_per_mult = (int) ceil ((double) (2 * BRUTE_FORCE_WORDS * 32 + 1) / log2 (gwdata->b));
		// Calc number of chunks and pairs
		num_chunks = (int) ceil (gwdata->bit_length / (BRUTE_FORCE_WORDS * 32));	// Total num chunks
		num_pairs = (num_chunks + 1) / 2;						// Number of pairs
		// Loop until we find a good FFT to use
		for (fftlen = 0; ; fftlen++) {
			double	words_per_mult;		/* Number of FFT words a specific FFT uses to hold a partial multiply result */
			int	res, option;

			gwinit (&new_gwdata);
			gwset_minimum_fftlen (&new_gwdata, fftlen);
			res = gwinfo (&new_gwdata, 1.0, gwdata->b, exp = num_pairs * b_per_mult, -1);
			if (res) {			// On error, try for a small rational FFT
				ASSERTG (fftlen == 0);
				fftlen = exp - 1;
				continue;
			}
			
			fftlen = new_gwdata.FFTLEN;
			words_per_mult = (double) b_per_mult * (double) fftlen / (double) exp;

			// Do we have an FFT we can use?  FFT words per mult must be an integer so that every partial multiplication
			// result has the same bigword/smallword pattern.  And the FFT must be large enough to hold all our pairs.
			if (words_per_mult == ceil (words_per_mult) && fftlen >= num_pairs * (int) words_per_mult) break;

			// To make words_per_mult an integer we have 3 options.  We can increase b_per_mult or exp to
			// reduce the words_per_mult.  Or we can increase word_per_mult as long as we can still store
			// num_pairs partial multiply results.  Or we can try a larger FFT length.

			// Option 1:  try smaller words_per_mult, increase b_per_word and exp to reduce words_per_mult to an integer
			// Option 2:  try larger words_per_mult as long as FFT can still hold num_pairs.
			// If either of these options results in the same FFT length, we have a winner
			for (option = 1; option <= 2; option++) {
				int newwpm = (int) floor (words_per_mult) + (option - 1);
				int newbpm = round_up_to_multiple_of (b_per_mult, newwpm / intgcd (fftlen, newwpm));
				if (newwpm > newbpm) newwpm = newbpm;		// When base is large, use rational FFT one b per FFT word
				if (fftlen >= num_pairs * newwpm) {
					gwinit (&new_gwdata);
					res = gwinfo (&new_gwdata, 1.0, gwdata->b, exp = newbpm * fftlen / newwpm, -1);
					if (res == 0 && fftlen >= (int) new_gwdata.FFTLEN) option = 99;	// Winner!
				}
			}
			if (option >= 99) break;

			// Option 3: Try a larger FFT length
		}
		gwinit (&new_gwdata);
		if (gwdata->num_threads > 1) gwset_num_threads (&new_gwdata, gwdata->num_threads);
		gwset_specific_fftlen (&new_gwdata, fftlen);
		gwsetup (&new_gwdata, 1.0, gwdata->b, exp, -1);
	}

// Allocate temporaries for our calculations

	pow2_multiplier = gwalloc (&new_gwdata);
	t1 = gwalloc (&new_gwdata);
	t3 = gwalloc (&new_gwdata);

// Calculate number of input chunks and number of FFT words in upper * pow2_multiplier + lower
// The number of FFT words must end on an integral number of b.

	num_chunks = divide_rounding_up (g->sign, BRUTE_FORCE_WORDS);
	fft_words_per_mult = (int) ceil (
			(2 * BRUTE_FORCE_WORDS * 32 + 1)	// Number of bits in multiplication result plus one to avoid borrow
			/ (log2 (gwdata->b) * new_gwdata.avg_num_b_per_word));	// Over bits per FFT word
	fft_words_per_mult = round_up_to_multiple_of (fft_words_per_mult, new_gwdata.FFTLEN / intgcd (new_gwdata.n, new_gwdata.FFTLEN));

// Pre-calculate the bigword vs. smallword flags so that brute_convert does not need to do that repeatedly.

	ASSERTG (fft_words_per_mult <= 64);
	big_word_flags = 0;
	num_big_word_flags = fft_words_per_mult;
	for (i = 0; i < num_big_word_flags; i++) {
		if (is_big_word (&new_gwdata, i)) big_word_flags |= ((uint64_t) 1) << i;
	}

// Brute force convert the first pow2_multiplier 2^(32*BRUTE_FORCE_WORDS)

	{
		stackgiant(tmpg,BRUTE_FORCE_WORDS+1);
		setone (tmpg);
		gshiftleft (32*BRUTE_FORCE_WORDS, tmpg);
		ulsubg (1, tmpg);
		dbltogw (&new_gwdata, 0.0, pow2_multiplier);
		brute_convert (&new_gwdata, tmpg, pow2_multiplier, 0, big_word_flags);
		*pow2_multiplier += 1.0;
	}
#ifdef GDEBUG
	gwsetnormroutine (&new_gwdata,0,1,0);
#endif

// In the first round we convert from binary (BRUTE_FORCE_WORDS * 32 bits) to radix b.  We do this in lower / upper pairs that we combine
// by calculating upper * pow2_multiplier + lower.  We store the result in t1.

	int	giant_offset, gwnum_offset, num_left;

	num_chunks = (num_chunks + 1) / 2;			// This now represents the number of pairs to process
	giant_offset = 0;
	dbltogw (&new_gwdata, 0.0, t3);				// We'll convert the lower halves to t3
	dbltogw (&new_gwdata, 0.0, t1);				// We convert the upper halves to t1

	gwnum_offset = 0;
	for (num_left = num_chunks; num_left; num_left--) {
		stackgiant(tmpg,BRUTE_FORCE_WORDS);

		giant_extract (g, giant_offset, BRUTE_FORCE_WORDS, tmpg);
	        brute_convert (&new_gwdata, tmpg, t3, gwnum_offset, big_word_flags);

		giant_extract (g, giant_offset + BRUTE_FORCE_WORDS, BRUTE_FORCE_WORDS, tmpg);
	        brute_convert (&new_gwdata, tmpg, t1, gwnum_offset, big_word_flags);

		gwnum_offset += fft_words_per_mult;
		giant_offset += BRUTE_FORCE_WORDS * 2;
	}

	gwmul (&new_gwdata, pow2_multiplier, t1);		// Multiply the upper halves by the power of two multiplier
	gwaddquick (&new_gwdata, t3, t1);			// Add the lower halves to the multiplied upper halves

// Now we take our gwnum (t1) holding data converted to radix b and combine lower / upper pairs until we
// get down to just one fully radix-converted value.

	while (num_chunks > 1) {
		int	j, num_left, gwnum_offset;

		// Double the exponent of the pow2_multiplier
		gwstartnextfft (&new_gwdata, 1);
		gwsquare (&new_gwdata, pow2_multiplier);
		gwstartnextfft (&new_gwdata, 0);

		// Process pairs of chunks.  Move FFT words from upper chunk of pair to t3.
		num_chunks = (num_chunks + 1) / 2;			// This now represents the number of pairs in t1
		dbltogw (&new_gwdata, 0.0, t3);
		for (gwnum_offset = 0, num_left = num_chunks; num_left; num_left--, gwnum_offset += fft_words_per_mult * 2) {
			for (j = gwnum_offset; j < gwnum_offset + fft_words_per_mult; j++) {
				double	*upper_word;
				if (j + fft_words_per_mult >= (int) new_gwdata.FFTLEN) break;
				upper_word = addr (&new_gwdata, t1, j + fft_words_per_mult);
				* addr (&new_gwdata, t3, j) = *upper_word;
				*upper_word = 0.0;
			}
		}

		gwmul (&new_gwdata, pow2_multiplier, t3);		// Apply the power of two multiplier to the upper halves
		gwaddquick (&new_gwdata, t3, t1);			// Add the multiplied upper and the lower
		fft_words_per_mult = fft_words_per_mult * 2;
	}
	ASSERTG (gw_get_maxerr(&new_gwdata) < 0.43);

/* Cleanup a little */

	gwfree (&new_gwdata, pow2_multiplier);
	gwfree (&new_gwdata, t3);

// Finally, convert the converted result (t1) to the caller's gwdata.

	radix_gwcopy (&new_gwdata, gwdata, t1, x, big_word_flags, num_big_word_flags);

/* Finish cleanup and return */

	gwfree (&new_gwdata, t1);
	gwdone (&new_gwdata);
}


// Brute force conversion from binary (stored in a giant) to part of a gwnum using a different radix

void brute_convert (
	gwhandle *gwdata,
	giant	g,			// Input giant (not preserved)
	gwnum	x,			// Output  gwnum
	int	offset,			// FFT word to start storing the converted giant
	uint64_t big_word_flags)	// Big word vs. small word flags
{
	int	mask1, mask2;

	mask1 = intpow ((int) gwdata->b, (int) gwdata->NUM_B_PER_SMALL_WORD);
	mask2 = gwdata->b * mask1;

/* Convert the giant to FFT format using uint64_t */

#if BRUTE_FORCE_WORDS <= 2

	if (g->sign == 0) return;
	ASSERTG (g->sign <= 2);
	{
		uint64_t g64, q;

		if (g->sign == 2) g64 = (((uint64_t) g->n[1]) << 32) + g->n[0];
		else g64 = g->n[0];
		while (g64) {
			int	mask;
			long	value;

//			mask = is_big_word (gwdata, offset) ? mask2 : mask1;
			mask = (big_word_flags & 1) ? mask2 : mask1;
			big_word_flags >>= 1;

			q = g64 / mask;				// quotient
			value = (long) (g64 - q * mask);	// value = remainder
			set_fft_value (gwdata, x, offset++, value);
			g64 = q;				// g64 = quotient
		}
	}

/* Convert the giant to FFT format */

#else
	{
		stackgiant(tmpg,BRUTE_FORCE_WORDS+1);
		stackgiant(q,BRUTE_FORCE_WORDS+1);
		while (g->sign) {
			int	mask;
			long	value;

//			mask = is_big_word (gwdata, offset) ? mask2 : mask1;
			mask = (big_word_flags & 1) ? mask2 : mask1;
			big_word_flags >>= 1;

			itog (mask, tmpg);			// tmpg = divisor
			gtog (g, q);
			divg (tmpg, q);				// q = quotient
			mulgi (&gwdata->gdata, q, tmpg);
			subg (tmpg, g);				// g = remainder
			value = (g->sign) ? g->n[0] : 0;	// value = remainder
			set_fft_value (gwdata, x, offset++, value);
			gtog (q, g);				// g = quotient
		}
	}
#endif

}

// Copy gwnum's from one gwhandle to another gwhandle using the same base

void radix_gwcopy (
	gwhandle *src_gwdata,
	gwhandle *dst_gwdata,
	gwnum	src,
	gwnum	dst,
	uint64_t big_word_flags,
	int	num_big_word_flags)
{
	int	srci, dsti, i, last_dsti;
	long	srcval, dstval, carry;
	int	srcb, dstb;
	int	radix_powers[32];

	radix_powers[0] = 1;
	radix_powers[1] = dst_gwdata->b;
	for (i = 2; i <= (int) dst_gwdata->NUM_B_PER_SMALL_WORD + 1; i++) radix_powers[i] = radix_powers[1] * radix_powers[i-1];

	srci = 0;
	srcb = 0;
	carry = 0;
	last_dsti = dst_gwdata->FFTLEN-1;
	if (dst_gwdata->NUM_B_PER_SMALL_WORD == 0) while (!is_big_word (dst_gwdata, last_dsti)) last_dsti--;
	for (dsti = 0; dsti < (int) dst_gwdata->FFTLEN; dsti++) {
		int	b_we_have;

		dstb = dst_gwdata->NUM_B_PER_SMALL_WORD + !!is_big_word (dst_gwdata, dsti);
		dstval = 0;
		b_we_have = 0;

		while (b_we_have < dstb) {
			int	b_to_get, remainder;

			// Get more src data if we do not have any
			if (srcb == 0) {
				// If no more src words, we're done
				if (srci == src_gwdata->FFTLEN) break;
				// Read in next source value
				get_fft_value (src_gwdata, src, srci, &srcval);
				srcb = src_gwdata->NUM_B_PER_SMALL_WORD + ((big_word_flags >> (srci % num_big_word_flags)) & 1);
				srci++;
			}

			// Get part or all of source word
			b_to_get = dstb - b_we_have;
			if (b_to_get < srcb) {
				remainder = srcval % radix_powers[b_to_get];
				srcval = srcval / radix_powers[b_to_get];
				srcb -= b_to_get;
			} else {
				remainder = srcval;
				b_to_get = srcb;
				srcb = 0;
			}

			// Add the just fetched src data to the destination
			dstval += remainder * radix_powers[b_we_have];
			b_we_have += b_to_get;
		}

		// Store in balanced representation except for the last destination FFT word.
		dstval += carry;
		if (dsti == last_dsti) {			// NOTE: This occurs before the last FFT word if avg_num_b_per_word < 1
			if (dstval < 0) dstval += radix_powers[dstb];
		}
		else if (dstval < -radix_powers[dstb] / 2) {
			dstval += radix_powers[dstb];
			carry = -1;
		}
		else if (dstval > radix_powers[dstb] / 2) {
			dstval -= radix_powers[dstb];
			carry = 1;
		} else
			carry = 0;
		set_fft_value (dst_gwdata, dst, dsti, dstval);
	}
}
