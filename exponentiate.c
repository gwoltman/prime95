/*----------------------------------------------------------------------
| Copyright 2020-2021 Mersenne Research, Inc.  All rights reserved
|
| Auxiliary routines to exponentiate gwnums
+---------------------------------------------------------------------*/

/* Includes */

#include "common.h"		// Included in all GIMPS sources
#include "exponentiate.h"

// This version is a simple left-to-right powering

// Raise a gwnum to a power
void simple_exponentiate (gwhandle *gwdata, gwnum x, uint64_t power)
{
	uint64_t current_bit;
	int	first_squaring;
	gwnum	scratch;		// Scratch gwnum for code to use

	// Handle simple cases
	if (power == 0) { dbltogw (gwdata, 1.0, x); return; }
	if (power == 1) return;

	scratch = gwalloc (gwdata);
//BUG	if (scratch == NULL) goto oom;

	// Find the topmost bit in power
	for (current_bit = 0x8000000000000000ULL; power < current_bit; current_bit >>= 1);

	// Exponentiate
	first_squaring = 1;
	for (current_bit >>= 1; current_bit; current_bit >>= 1) {
		if (first_squaring) {
			gwfft (gwdata, x, scratch);
			gwmul3 (gwdata, scratch, scratch, x, (current_bit != 1 || (power & 1) == 0) ? GWMUL_STARTNEXTFFT : 0);
			first_squaring = 0;
		} else
			gwsquare2 (gwdata, x, x, (current_bit != 1 || (power & 1) == 0) ? GWMUL_STARTNEXTFFT : 0);
		if (power & current_bit)
			gwmul3 (gwdata, scratch, x, x, current_bit != 1 ? GWMUL_STARTNEXTFFT : 0);
	}

	// Free allocated memory
	gwfree (gwdata, scratch);
}


// Utility routine to help with bit testing

#define is_bit_set(array,bitnum)	((array)[(bitnum)>>6] & (1ULL << ((bitnum) & 63)))

// Raise a gwnum to a power.  This version uses a 3-bit windowing scheme to reduce multiplies a little

void exponentiate_window3 (gwhandle *gwdata, gwnum x, uint64_t *power, int arraylen)
{
	int	current_bit, top3, first_squaring_done;
	gwnum	x1, x3, x5, x7;		// Scratch gwnums for code to use
	int	x7_is_really_x2;

	// Allocate space for temporaries
	x1 = gwalloc (gwdata);
	x3 = gwalloc (gwdata);
	x5 = gwalloc (gwdata);
	x7 = gwalloc (gwdata);
//BUG	if (x1 == NULL || x3 == NULL || x5 == NULL || x7 == NULL) goto oom;

	// Pre-compute x^1, x^3, x^5.  We delay computing x^7 until we are sure we will need it.
	gwfft (gwdata, x, x1);
	gwmul3 (gwdata, x1, x1, x7, GWMUL_STARTNEXTFFT); x7_is_really_x2 = 1;
	gwmul3 (gwdata, x7, x1, x3, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
	gwmul3 (gwdata, x7, x3, x5, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);

	// Find the topmost bit in power
	current_bit = arraylen * 64 - 1;
	while (!is_bit_set (power, current_bit)) current_bit--;

	// Special case the first 3 bits of the exponent
	// 100:  Normal would be 1,2,4,8.  We do 1,2,3,5,8		Normal is better only if x3,x5 never used
	// 101:  Normal would be 1,2,4,5.  We do 1,2,3,5
	// 110:  Normal would be 1,2,3,6,12.  We do 1,2,3,5,7,12	Normal is better only if x5,x7 never used
	// 111:  Normal would be 1,2,3,6,7.  We do 1,2,3,5,7
	top3 = 4;
	if (is_bit_set (power, current_bit - 1)) top3 += 2;
	if (is_bit_set (power, current_bit - 2)) top3 += 1;
	if (top3 == 4) gwmul3 (gwdata, x3, x5, x, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	if (top3 == 5) gwmul3 (gwdata, x5, x5, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	if (top3 == 6) {
		gwmul3 (gwdata, x7, x5, x7, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT); x7_is_really_x2 = 0;
		gwmul3 (gwdata, x5, x7, x, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	}
	if (top3 == 7) {
		gwmul3 (gwdata, x7, x5, x7, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT); x7_is_really_x2 = 0;
		gwmul3 (gwdata, x7, x7, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	}
	first_squaring_done = 1;	// Set flag indicating we've already done the first squaring needed in exponentiate loop

	// Sliding window exponentiate
	for (current_bit -= 3; current_bit >= 0; current_bit--) {
		if (first_squaring_done) first_squaring_done = 0;
		else gwsquare2 (gwdata, x, x, current_bit || is_bit_set (power, 0) ? GWMUL_STARTNEXTFFT : 0);
		if (!is_bit_set (power, current_bit)) continue;
		// We have a window starting with a set bit, see if we want to apply x^1, x^3, x^5, or x^7
		if (current_bit >= 2 && is_bit_set (power, current_bit - 2)) {
			gwsquare2 (gwdata, x, x, GWMUL_STARTNEXTFFT);
			gwsquare2 (gwdata, x, x, GWMUL_STARTNEXTFFT);
			if (is_bit_set (power, current_bit - 1)) {	// Apply x^7
				if (x7_is_really_x2) { gwmul3 (gwdata, x7, x5, x7, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT); x7_is_really_x2 = 0; }
				gwmul3 (gwdata, x7, x, x, GWMUL_FFT_S1 | (current_bit != 2 ? GWMUL_STARTNEXTFFT : 0));
			} else {					// Apply x^5
				gwmul3 (gwdata, x5, x, x, GWMUL_FFT_S1 | (current_bit != 2 ? GWMUL_STARTNEXTFFT : 0));
			}
			current_bit -= 2;
		} else {
			if (current_bit >= 1 && is_bit_set (power, current_bit - 1)) {	// Apply x^3
				gwsquare2 (gwdata, x, x, GWMUL_STARTNEXTFFT);
				gwmul3 (gwdata, x3, x, x, current_bit != 1 ? GWMUL_STARTNEXTFFT : 0);
				current_bit -= 1;
			} else {				// Apply x^1
				gwmul3 (gwdata, x1, x, x, current_bit != 1 ? GWMUL_STARTNEXTFFT : 0);
			}
		}
	}

	// Free allocated memory
	gwfree (gwdata, x1);
	gwfree (gwdata, x3);
	gwfree (gwdata, x5);
	gwfree (gwdata, x7);
}

// Raise a gwnum to a power.  This version uses a 4-bit windowing scheme to reduce multiplies even more.

void exponentiate_window4 (gwhandle *gwdata, gwnum x, uint64_t *power, int arraylen)
{
	gwnum	x1, x3, x5, x7, x9, x11, x13, x15;		// Scratch gwnums for code to use
	int	current_bit, top4;
	bool	first_squaring_done;

	// Allocate space for temporaries
	x1 = gwalloc (gwdata);
	x3 = gwalloc (gwdata);
	x5 = gwalloc (gwdata);
	x7 = gwalloc (gwdata);
	x9 = gwalloc (gwdata);
	x11 = gwalloc (gwdata);
	x13 = gwalloc (gwdata);
	x15 = gwalloc (gwdata);
//BUG	if (x1 == NULL || x3 == NULL || x5 == NULL || x7 == NULL || x9 == NULL || x11 == NULL || x13 == NULL || x15 == NULL) goto oom;

	// Pre-compute x^1, x^3, x^5, ... x^15
	gwfft (gwdata, x, x1);
	gwmul3 (gwdata, x1, x1, x15, GWMUL_STARTNEXTFFT);
	gwmul3 (gwdata, x15, x1, x3, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
	gwmul3 (gwdata, x15, x3, x5, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	gwmul3 (gwdata, x15, x5, x7, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	gwmul3 (gwdata, x15, x7, x9, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	gwmul3 (gwdata, x15, x9, x11, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	gwmul3 (gwdata, x15, x11, x13, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	gwmul3 (gwdata, x15, x13, x15, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);

	// Find the topmost bit in power
	current_bit = arraylen * 64 - 1;
	while (!is_bit_set (power, current_bit)) current_bit--;

	// Special case the first 4 bits of the exponent.
	top4 = 8;
	if (is_bit_set (power, current_bit - 1)) top4 += 4;
	if (is_bit_set (power, current_bit - 2)) top4 += 2;
	if (is_bit_set (power, current_bit - 3)) top4 += 1;
	if (top4 == 8) gwmul3 (gwdata, x7, x9, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// x^16
	if (top4 == 9) gwmul3 (gwdata, x9, x9, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// x^18
	if (top4 == 10) gwmul3 (gwdata, x9, x11, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// x^20
	if (top4 == 11) gwmul3 (gwdata, x11, x11, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// x^22
	if (top4 == 12) gwmul3 (gwdata, x11, x13, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// x^24
	if (top4 == 13) gwmul3 (gwdata, x13, x13, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// x^26
	if (top4 == 14) gwmul3 (gwdata, x13, x15, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// x^28
	if (top4 == 15) gwmul3 (gwdata, x15, x15, x, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// x^30
	first_squaring_done = 1;	// Set flag indicating we've already done the first squaring needed in exponentiate loop

	// Sliding window exponentiate
	for (current_bit -= 4; current_bit >= 0; current_bit--) {
		if (first_squaring_done) first_squaring_done = 0;
		else gwsquare2 (gwdata, x, x, current_bit || is_bit_set (power, 0) ? GWMUL_STARTNEXTFFT : 0);
		if (!is_bit_set (power, current_bit)) continue;
		// We have a window starting with a set bit, see if we want to apply x^1, x^3, x^5, x^7, x^9, x^11, x^13, or x^15
		if (current_bit >= 3 && is_bit_set (power, current_bit - 3)) {
			gwsquare2 (gwdata, x, x, GWMUL_STARTNEXTFFT);
			gwsquare2 (gwdata, x, x, GWMUL_STARTNEXTFFT);
			gwsquare2 (gwdata, x, x, GWMUL_STARTNEXTFFT);
			if (is_bit_set (power, current_bit - 2)) {	// Apply x^11 or x^15
				if (is_bit_set (power, current_bit - 1)) {	// Apply x^15
					gwmul3 (gwdata, x15, x, x, GWMUL_FFT_S1 | (current_bit != 3 ? GWMUL_STARTNEXTFFT : 0));
				} else {					// Apply x^11
					gwmul3 (gwdata, x11, x, x, GWMUL_FFT_S1 | (current_bit != 3 ? GWMUL_STARTNEXTFFT : 0));
				}
			}
			else if (is_bit_set (power, current_bit - 1)) {	// Apply x^13
				gwmul3 (gwdata, x13, x, x, GWMUL_FFT_S1 | (current_bit != 3 ? GWMUL_STARTNEXTFFT : 0));
			} else {					// Apply x^9
				gwmul3 (gwdata, x9, x, x, GWMUL_FFT_S1 | (current_bit != 3 ? GWMUL_STARTNEXTFFT : 0));
			}
			current_bit -= 3;
		}
		else if (current_bit >= 2 && is_bit_set (power, current_bit - 2)) {
			gwsquare2 (gwdata, x, x, GWMUL_STARTNEXTFFT);
			gwsquare2 (gwdata, x, x, GWMUL_STARTNEXTFFT);
			if (is_bit_set (power, current_bit - 1)) {	// Apply x^7
				gwmul3 (gwdata, x7, x, x, GWMUL_FFT_S1 | (current_bit != 2 ? GWMUL_STARTNEXTFFT : 0));
			} else {					// Apply x^5
				gwmul3 (gwdata, x5, x, x, GWMUL_FFT_S1 | (current_bit != 2 ? GWMUL_STARTNEXTFFT : 0));
			}
			current_bit -= 2;
		} else {
			if (current_bit >= 1 && is_bit_set (power, current_bit - 1)) {	// Apply x^3
				gwsquare2 (gwdata, x, x, GWMUL_STARTNEXTFFT);
				gwmul3 (gwdata, x3, x, x, current_bit != 1 ? GWMUL_STARTNEXTFFT : 0);
				current_bit -= 1;
			} else {				// Apply x^1
				gwmul3 (gwdata, x1, x, x, current_bit != 1 ? GWMUL_STARTNEXTFFT : 0);
			}
		}
	}

	// Free allocated memory
	gwfree (gwdata, x1);
	gwfree (gwdata, x3);
	gwfree (gwdata, x5);
	gwfree (gwdata, x7);
	gwfree (gwdata, x9);
	gwfree (gwdata, x11);
	gwfree (gwdata, x13);
	gwfree (gwdata, x15);
}

// Internal routine to call an exponentiate routine using the appropriate window size

void exponentiate_array (gwhandle *gwdata, gwnum x, uint64_t *power, int arraylen)
{
	// Handle simple cases.  In theory, window size 3 breaks even with simple exponentiation at around 9 bits.
	while (arraylen >= 1 && power[arraylen-1] == 0) arraylen--;
	if (arraylen <= 0) { simple_exponentiate (gwdata, x, 0); return; }
	if (arraylen == 1 && power[0] <= 256) { simple_exponentiate (gwdata, x, *power); return; }
	// In theory, window size 4 breaks even with window size 3 at about 64 bits
	if (arraylen == 1) exponentiate_window3 (gwdata, x, power, arraylen);
	else exponentiate_window4 (gwdata, x, power, arraylen);
}

// Raise a gwnum to a power
void exponentiate (gwhandle *gwdata, gwnum x, uint64_t power)
{
	exponentiate_array (gwdata, x, &power, 1);
}

// Raise a gwnum to a mpz power
void exponentiate_mpz (gwhandle *gwdata, gwnum x, mpz_t power)
{
	uint64_t *array;
	size_t	len;

	array = (uint64_t *) malloc (divide_rounding_up (mpz_sizeinbase (power, 2), 64) * sizeof (uint64_t));
	mpz_export (array, &len, -1, sizeof (uint64_t), 0, 0, power);
	exponentiate_array (gwdata, x, array, (int) len);
	free (array);
}

