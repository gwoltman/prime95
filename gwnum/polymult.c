/*----------------------------------------------------------------------
| polymult.c
|
| This file contains the C routines for fast polynomial multiplication
| 
|  Copyright 2021 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "cpuid.h"
#include "gwnum.h"
#include "gwutil.h"
#include "polymult.h"
#if defined (AVX) || defined (FMA) || defined (AVX512)
#include <immintrin.h>
#endif

// TO DO:
// 1) AVX-512, FMA, SSE2? versions
// 2) Zero padded gwnums.  Test all FFT sizes and FFT types -- especially funky PFA ffts.
// 7) Variable GW_HEADER_SIZE.
// 9) Is a fast convert-gwnum-to-larger-FFT feasible? necessary?
// 10) Tune for Brute / Karatsuba / FFT crossovers.
// 12) Tune routine that returns needed safety margin needed for given poly sizes
// 13) multi-threading
// 14) Compressed gwnums.  Does it pay off?
// 15) Add polymult to gwtest or as a standalone program
// 16) Handle out-of-memory gracefully
// 17) Sin/cos accuracy, PFA FFTs, 2-pass FFTs
// 18) DJB FFTs
// 19) For large ffts, use Pavel's idea of copying the complex values to a gwnum and use an unweighted FFT with a special no-carry-propagation norm routine.
//     Perhaps only do this for most common architectures, such as AVX-512 and FMA 64-bit.
// 20) Radix-4 non-vector FFT.
// 21) Should polymult be its own library or a part of the gwnum library?
// 22) Can we add asserts that polymult isn't violating EXTRA_BITS?
// 25) Optimize.  Read/write several cache lines at a time (fill up L2 cache?).  Prefetch future cache lines.
// 28) Optimize brute_force to do less compute when outvec entries are NULL.
// 28) Optimize brute_force to do less compute when one or both invec entries are RLPs.
// 29) monic polynomials are tricky if FFT(1) is not all ones.  Use GW_FFT1 and treat as non-monic?
// 30) option to negate outputs??
// 32) If monic input and polymult is doing an FFT and there is zero padding going on, then pad with one (or GW_FFT1) rather than use monic_add.  Be careful
//	that this does not cause any unwanted wrapping of data!
// 33) Are there cases with brute/Karatsuba would benefit from using a bigger poly rather than using monic_add?
// 34) Support GWMUL_BY_CONST, GWMUL_ADD_CONST option in gwunfft2.
// 35) A monic 1024 times monic 1024 using length 2048 FFT could apply the two ones in the FFT vectors and patch outvec[0] (possibly with gwsetaddin(-1)).
//	A monic 1023 times monic 1023 using length 2048 FFT could apply the two ones in the FFT vectors with no patching rather than monic  post-processing.
//	That should be a little faster.  Other monic combinations have similar optimizations.
// 36) For cached FFTs, is there a way to avoid the transpose(write_line/read_line) step?
// 37) If transpose is expensive, can we offer a way to "pre-transpose" invec1 for use multiple times?  With compression of exponent field in double??
//	Compression still needs semi-random/sequential access to invec1 for multi-threading.
// 38)  For poly #2, is there a quick block compression we can do where the extra read/write cost will be more than compensated by the few more extra
//	gwnums we can add to poly #2?
// 39)  Think more on unbalanced poly mults.  Is a 100 x 200 mult better done as two 100 x 100 multiplies?  Answer: Never for brute force, maybe for FFT.
// 40) Is invec1_negate code OK when GW_FFT1 is not all ones?  I would think so.
// 41) Add callback routine to check for ^C/ESC
// 42) Transparent / or easy way to store gwnums on disk as part of input and/or output vectors.  Picture using 12GB memory and 50GB disk for really huge polys.
// 43) Add callback to set affinity or use gwnum's callback
// 44) Cache created threads?

// Forward declarations required for the AVX, FMA, AVX512 compiles
int next_line (pmhandle *pmdata);
unsigned long cache_line_offset (gwhandle *gwdata, unsigned long i);
void generate_sincos (pmhandle *pmdata, const unsigned int size);
void polymult_line_avx (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, gwnum *gw_invec2, int invec2_size, gwnum *gw_outvec, int outvec_size, int options);

#define INVEC1_OPTIONS	(POLYMULT_INVEC1_MONIC | POLYMULT_INVEC1_RLP | POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC1_FFTED)
#define INVEC2_OPTIONS	(POLYMULT_INVEC2_MONIC | POLYMULT_INVEC2_RLP | POLYMULT_INVEC2_NEGATE | POLYMULT_INVEC2_FFTED)

#if !defined (AVX) && !defined (FMA) && !defined (AVX512)

// Get the needed safety_margin required for an invec1_size by invec2_size polymult
double polymult_safety_margin (int invec1_size, int invec2_size) {
	int	n = (invec1_size < invec2_size) ? invec1_size : invec2_size;
	// The number of EXTRA_BITS (output FFT bits) consumed by a polymult ~= log2(n)/2.  However, gwinit's safety margin parameter refers to
	// input FFT bits.  Thus, divide by another 2.
	return (log2 ((double) n) / 4.0);
}

// Get the FFT size that will be used for an n = invec1_size + invec2_size polymult
int polymult_fft_size (int n)
{
	int	fft_sizes[] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576,
				2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912};
	for (int i = 0; i < sizeof (fft_sizes) / sizeof (int); i++) if (fft_sizes[i] >= n) return (fft_sizes[i]);
	return (1073741824);
}

// Initialize a polymult handle
void polymult_init (
	pmhandle *pmdata,		// Handle for polymult library that needs initializing
	gwhandle *gwdata)		// Handle for gwnum FFT library
{
	memset (pmdata, 0, sizeof (pmhandle));
	pmdata->gwdata = gwdata;

	// Initialize multithreading mutex and default number of threads
	gwmutex_init (&pmdata->poly_mutex);
	pmdata->num_threads = gwget_num_threads (gwdata);

	// Initialize default breakeven points
	if (pmdata->gwdata->cpu_flags & CPU_AVX) {
		pmdata->KARAT_BREAK = 32;			//GW: Fix me
		pmdata->FFT_BREAK = 64;				//GW: Fix me
	}
	else {
		pmdata->KARAT_BREAK = 32;			//GW: Fix me
		pmdata->FFT_BREAK = 64;				//GW: Fix me
	}
}

// Terminate use of a polymult handle.  Free up memory.
void polymult_done (
	pmhandle *pmdata)		// Handle for polymult library
{
	gwmutex_destroy (&pmdata->poly_mutex);
	free (pmdata->twiddle_real);
	free (pmdata->twiddle_imag);
	free (pmdata->twiddle_real4);
	free (pmdata->twiddle_imag4);
}


// Return the next line that needs processing
int next_line (
	pmhandle *pmdata)		// Handle for polymult library
{
	int	return_value;

	// Use mutex to guard next_line counter, return next available line
	gwmutex_lock (&pmdata->poly_mutex);
	return_value = pmdata->next_line++;
	gwmutex_unlock (&pmdata->poly_mutex);
	return (return_value);
}

//GW -- move to gwnum.c?

/* Get offset to nth "cache line" from an FFTed gwnum.  Where a cache line is 128 bytes for an AVX-512 FFT, 64 bytes otherwise. */

unsigned long cache_line_offset (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long i)
{
	unsigned long addr, grps;

/* Memory layouts for AVX-512 machines */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		if (gwdata->PASS1_SIZE == 0) {
			addr = i * 128;
			/* Now optionally add 64 pad bytes every 4KB */
			if (gwdata->FOURKBGAPSIZE) addr = addr + (addr / (gwdata->FOURKBGAPSIZE << 6)) * 64;
		}
		else {
			grps = i / (gwdata->PASS2_SIZE / 8);
			addr = i * 128;
			/* Now optionally add bytes for every 4KB page and every pass 2 block */
			addr = addr + (addr >> 12) * gwdata->FOURKBGAPSIZE + grps * gwdata->PASS2GAPSIZE;
		}
	}

/* Memory layouts for AVX machines */

	else if (gwdata->cpu_flags & CPU_AVX) {
		if (gwdata->PASS2_SIZE == 0) {
			addr = i * 64;
			/* Now optionally add 64 pad bytes every 1KB, 2KB or 4KB */
			if (gwdata->FOURKBGAPSIZE) addr = addr + (addr / (gwdata->FOURKBGAPSIZE << 6)) * 64;
		}
		else {
			grps = i / (gwdata->PASS2_SIZE / 4);
			addr = i * 64;
			/* Now add 64 bytes every 4KB and one pass2gapsize for every pass 2 block. */
			addr = addr + (addr >> 12) * gwdata->FOURKBGAPSIZE + grps * gwdata->PASS2GAPSIZE;
		}
	}

/* Memory layouts for SSE2 machines */

	else if (gwdata->cpu_flags & CPU_SSE2) {

/* PFA-style FFTs are a little tricker.  See assembly code for example.	*/
		if (gwdata->PASS2_SIZE == 0) {
			addr = i * 64;
		}
//GW: BUG!?
/* and PFA layouts are even funkier.					*/
		else if (gwdata->FFT_TYPE == FFT_TYPE_HOME_GROWN) {
			grps = i / (gwdata->PASS2_SIZE / 2);
			addr = i * 64;
			/* Now add 128 bytes every 8KB and one pass2gapsize for every pass 2 block. */
			addr = addr + (addr >> 13) * 128 + grps * gwdata->PASS2GAPSIZE;
		}
/* Newer traditional radix-4 large FFTs use don't have a special layout for PFA. */
		else {
			grps = i / (gwdata->PASS2_SIZE / 2);
			addr = i * 64;
			/* Now add 128 bytes every 8KB and one pass2gapsize for every pass 2 block. */
			addr = addr + (addr >> 13) * 128 + grps * gwdata->PASS2GAPSIZE;
		}
	}

/* One pass x87 FFTs use a near flat memory model. */

//BUG --- complex vals are 8 away from reals (at least in integer space) ---- not an issue if polymult is only for x64
	else if (gwdata->PASS2_SIZE == 0) {
		addr = i * 64;
	}
	else {
		grps = (i * 8) / gwdata->PASS2_SIZE;
		addr = i * 64 + grps * 64;
	}

/* Return the offset */

	return (addr);
}


// Helper routines for FFT code

// Bit reverse a number
unsigned int bitrev (unsigned int i, unsigned int size)
{
	unsigned int res;
	for (res = 0; size > 1; size >>= 1) {
		res = (res << 1) + (i & 1);
		i >>= 1;
	}
	return (res);
}

// Generate radix-2 and radix-4 twiddle factors
void generate_sincos (
	pmhandle *pmdata,		// Handle for polymult library
	const unsigned int size)
{
	const double twopi = 2.0 * 3.1415926535897932384626433;

	// Return quickly if twiddles have already been initialized
	if (pmdata->twiddles_initialized && pmdata->twiddles_initialized == size) return;

	// Make sure only one thread initializes the twiddles
	gwmutex_lock (&pmdata->poly_mutex);

	// If some other thread beat us to the initialization, then return
	if (pmdata->twiddles_initialized && pmdata->twiddles_initialized == size) {
		gwmutex_unlock (&pmdata->poly_mutex);
		return;
	}

	// Free the previous twiddles
	if (pmdata->twiddles_initialized) {
//GW: accept size < pmdata->twiddles_initialized (at least for power of 2)
		free (pmdata->twiddle_real), pmdata->twiddle_real = NULL;
		free (pmdata->twiddle_imag), pmdata->twiddle_imag = NULL;
		free (pmdata->twiddle_real4), pmdata->twiddle_real4 = NULL;
		free (pmdata->twiddle_imag4), pmdata->twiddle_imag4 = NULL;
	}

//GW: Consider long doubles for better accuracy
	pmdata->twiddle_real = (double *) malloc (size/2 * sizeof (double));
	pmdata->twiddle_imag = (double *) malloc (size/2 * sizeof (double));
	for (unsigned int i = 0; i < size/2; i++) {
		double angle = twopi * (double) bitrev (i, size/2) / (double) size;
		pmdata->twiddle_real[i] = cos (angle);
		pmdata->twiddle_imag[i] = sin (angle);
	}

//GW: Consider DJB twiddles to reduce table size by 33%
	pmdata->twiddle_real4 = (double *) malloc (size/4 * 3 * sizeof (double));
	pmdata->twiddle_imag4 = (double *) malloc (size/4 * 3 * sizeof (double));
	for (unsigned int i = 0; i < size/4; i++) {
		double angle = twopi * (double) bitrev (i, size/4) / (double) size;
		pmdata->twiddle_real4[i*3] = cos (angle);
		pmdata->twiddle_imag4[i*3] = sin (angle);
		pmdata->twiddle_real4[i*3+1] = cos (2.0*angle);
		pmdata->twiddle_imag4[i*3+1] = sin (2.0*angle);
		pmdata->twiddle_real4[i*3+2] = cos (3.0*angle);
		pmdata->twiddle_imag4[i*3+2] = sin (3.0*angle);
	}

	pmdata->twiddles_initialized = size;
	gwmutex_unlock (&pmdata->poly_mutex);
}

#endif


//--------------------------------------------------------------------
//		Not-vector-optimized implementation
//--------------------------------------------------------------------

#if !defined (AVX) && !defined (FMA) && !defined (AVX512)

bool read_line (gwhandle *gwdata, gwnum *data, int size, int index, double *vec, int options)
{
	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += 2*(size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// If we've read all the gwnum data, return FALSE
	if (index * 2 >= (int) gwdata->FFTLEN) return (FALSE);

	// No special instructions implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(double)); continue; }
			vec[2*j] = data[j][0];
//bug? return 0.0 instead of -11 ---- fix when we support FFTing invec
			vec[2*j+1] = 0.0; //data[j][-11];
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(double)); continue; }
			if (gwdata->PASS2_SIZE == 0) vec[2*j] = data[j][1];
			else vec[2*j] = data[j][2];
//bug? return 0.0 instead of -10.  How are we going to implement FFTed poly read?  Its own routine? It's own FFT state? Compressed?
			vec[2*j+1] = 0.0; //data[j][-10];
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 4) / sizeof (double);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(double)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
// bug???  Index assuming what?  skip every other set of 4?  Examine ancient FFT sources!
			if (gwdata->PASS2_SIZE == 0 && index < 8 && !gwdata->ALL_COMPLEX_FFT) {
				vec[2*j] = data[j][offset+(index&3)*2];
				vec[2*j+1] = data[j][offset+(index&3)*2+1];
			} else {
				vec[2*j] = data[j][offset+(index&2)*2+(index&1)];
				vec[2*j+1] = data[j][offset+(index&2)*2+(index&1)+2];
			}
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		for (int j = 0; j < size; j++) {
			vec[2*j] = -vec[2*j];
			vec[2*j+1] = -vec[2*j+1];
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) {
			vec[-2*j] = vec[2*j];
			vec[-2*j+1] = vec[2*j+1];
		}
	}

	// Return success
	return (TRUE);
}

void write_line (gwhandle *gwdata, gwnum *data, int size, int index, double *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments_avx
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) || ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP)))
		vec -= 2;

	// If both inputs are RLPs, then the output is an RLP too.  Skip writing half of the data.
	if (options & POLYMULT_INVEC1_RLP && options & POLYMULT_INVEC2_RLP) vec += 2*(size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// No special instructions implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			data[j][0] = vec[2*j];
//			data[j][-11] = vec[2*j+1];  //GW:bug?
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {		// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (gwdata->PASS2_SIZE == 0) data[j][1] = vec[2*j];
			else data[j][2] = vec[2*j];
//			data[j][-10] = vec[2*j+1];	//GW:bug?
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 4) / sizeof (double);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (gwdata->PASS2_SIZE == 0 && index < 8 && !gwdata->ALL_COMPLEX_FFT) {
				data[j][offset+(index&3)*2] = vec[2*j];
				data[j][offset+(index&3)*2+1] = vec[2*j+1];
			} else {
				data[j][offset+(index&2)*2+(index&1)] = vec[2*j];
				data[j][offset+(index&2)*2+(index&1)+2] = vec[2*j+1];
			}
		}
	}
}

// A brute force O(N^2) implementation of polynomial multiplication
void brute_line (
	double	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	double	*invec2,		// Second input poly
	int	invec2_size,		// Size of the first input polynomial
	double	*outvec)		// Output poly
{
	// No special instructions implementation
	for (int i = 0; i < invec1_size; i++) {
		double real_x = invec1[2*i];
		double imag_x = invec1[2*i+1];
		for (int j = 0; j < invec2_size; j++) {
			double real_y = invec2[2*j];
			double imag_y = invec2[2*j+1];
			double real_tmp = real_x * real_y - imag_x * imag_y;
			double imag_tmp = real_x * imag_y + imag_x * real_y;
			if (i == 0 || j == invec2_size-1) {
				outvec[2*(i+j)] = real_tmp;
				outvec[2*(i+j)+1] = imag_tmp;
			} else {
				outvec[2*(i+j)] += real_tmp;
				outvec[2*(i+j)+1] += imag_tmp;
			}
		}
	}
}

// A recursive Karatsuba O(N^1.565) implementation of polynomial multiplication
//	a	b	(input vector 1)
//	c	d	(input vector 2)
//-----------------   
//ac   ad+bc   bd
//
//(a+b)*(c+d) = ac + ad + bc + bd
// ad+bc = (a+b)*(c+d) - ac - bd
void karatsuba_line (
	pmhandle *pmdata,		// Handle for polymult library
	double	*invec1,		// Coefficients of first input poly
	int	invec1_size,		// Size of the first input polynomial
	double	*invec2,		// Coefficients of second input poly
	int	invec2_size,		// Size of the first input polynomial
	double	*tmp,			// Temporary space
	double	*outvec)		// Output poly coefficients
{
	int	i, a_size, b_size, c_size, d_size, a_plus_b_size, c_plus_d_size, bd_size, ac_size, a_plus_b_times_c_plus_d_size;
	double	*a, *b, *c, *d, *ac, *bd, *a_plus_b, *c_plus_d, *a_plus_b_times_c_plus_d;

	// Handle end of recursion
	if (invec1_size == 1 || invec2_size == 1 || invec1_size + invec2_size < pmdata->KARAT_BREAK) {
		brute_line (invec1, invec1_size, invec2, invec2_size, outvec);
		return;
	}

	// To simplify code, always make invec2 the larger poly.
	if (invec1_size > invec2_size) {
		int tmpsize = invec1_size; invec1_size = invec2_size; invec2_size = tmpsize;
		double *tmpvec = invec1; invec1 = invec2; invec2 = tmpvec;
	}

	// Split larger vector roughly in half.  For readability, create pointers to two halves.
	d_size = (invec2_size + 1) / 2;
	c_size = invec2_size - d_size;
	c = invec2 + 2*d_size, d = invec2;

	// Handle case where vectors are of wildly unequal lengths
	if (d_size >= invec1_size) {
		// Compute bd
		bd = outvec;
		karatsuba_line (pmdata, invec1, invec1_size, d, d_size, tmp, bd);
		bd_size = invec1_size + d_size - 1;
		// Compute bc and store it in tmp
		double *bc = tmp;
		int bc_size = invec1_size + c_size - 1;
		tmp += 2*bc_size;
		karatsuba_line (pmdata, invec1, invec1_size, c, c_size, tmp, bc);
		// add bc to output
		for (i = 0; i < 2*bc_size; i++) outvec[i + 2*d_size] += bc[i];
		return;
	}

	// Split smaller vector, the size of the two lower halves must be the same.  For readability, create pointers to two halves.
	b_size = d_size;
	a_size = invec1_size - b_size;
	a = invec1 + 2*b_size, b = invec1;

	// Compute (a+b)
	a_plus_b = outvec;
	for (i = 0; i < 2*a_size && i < 2*b_size; i++) a_plus_b[i] = a[i] + b[i];
//	for ( ; i < 2*a_size; i++) a_plus_b[i] = a[i];		// In case a_size > b_size
	for ( ; i < 2*b_size; i++) a_plus_b[i] = b[i];		// In case b_size > a_size
	a_plus_b_size = i / 2;
	// Compute (c+d)
	c_plus_d = &outvec[i];
	for (i = 0; i < 2*c_size && i < 2*d_size; i++) c_plus_d[i] = c[i] + d[i];
//	for ( ; i < 2*c_size; i++) c_plus_d[i] = c[i];		// In case c_size > d_size
	for ( ; i < 2*d_size; i++) c_plus_d[i] = d[i];		// In case d_size > c_size
	c_plus_d_size = i / 2;
	// Compute (a+b) * (c+d) and store it in tmp
	a_plus_b_times_c_plus_d = tmp;
	a_plus_b_times_c_plus_d_size = a_plus_b_size + c_plus_d_size - 1;
	tmp += 2*a_plus_b_times_c_plus_d_size;
	karatsuba_line (pmdata, a_plus_b, a_plus_b_size, c_plus_d, c_plus_d_size, tmp, a_plus_b_times_c_plus_d);

	// Compute bd
	bd = outvec;
	karatsuba_line (pmdata, b, b_size, d, d_size, tmp, bd);
	bd_size = b_size + d_size - 1;
	// Output the complex value 0+0i
	outvec[2*bd_size] = 0;
	outvec[2*bd_size+1] = 0;
	// Compute ac
	ac = &outvec[2*(bd_size+1)];
	karatsuba_line (pmdata, a, a_size, c, c_size, tmp, ac);
	ac_size = a_size + c_size - 1;

	// subtract ac from middle
	for (i = 0; i < 2*ac_size; i++) a_plus_b_times_c_plus_d[i] -= ac[i];
	// subtract bd from middle
	for (i = 0; i < 2*bd_size; i++) a_plus_b_times_c_plus_d[i] -= bd[i];
	// add middle to output
	for (i = 0; i < 2*a_plus_b_times_c_plus_d_size; i++) outvec[i + 2*b_size] += a_plus_b_times_c_plus_d[i];
}

void fft_line (
	pmhandle *pmdata,		// Handle for polymult library
	double	*vec,			// Array of complex values
	const unsigned int size)	// The fft size
{
	generate_sincos (pmdata, size);

	// Do first two levels.  We know the high half of the vector is zero.
	ASSERTG (size > 2);
	unsigned int stride = size/4;
	for (unsigned int group_member = 0; group_member < stride; group_member++) {
		double real_a = vec[2*(group_member)];
		double imag_a = vec[2*(group_member) + 1];
		double real_b = vec[2*(group_member + stride)];
		double imag_b = vec[2*(group_member + stride) + 1];
		vec[2*(group_member)] = real_a + real_b;
		vec[2*(group_member) + 1] = imag_a + imag_b;
		vec[2*(group_member + stride)] = real_a - real_b;
		vec[2*(group_member + stride) + 1] = imag_a - imag_b;
		vec[2*(group_member + 2*stride)] = real_a - imag_b;
		vec[2*(group_member + 2*stride) + 1] = imag_a + real_b;
		vec[2*(group_member + 3*stride)] = real_a + imag_b;
		vec[2*(group_member + 3*stride) + 1] = imag_a - real_b;
	}
	stride >>= 1;

	// Radix 2 FFT
	for ( ; stride >= 1; stride >>= 1) {
		double *twiddle_real_ptr = pmdata->twiddle_real;
		double *twiddle_imag_ptr = pmdata->twiddle_imag;
		for (unsigned int group = 0; group < size; group += stride * 2, twiddle_real_ptr++, twiddle_imag_ptr++) {
			for (unsigned int group_member = 0; group_member < stride; group_member++) {
				double *real_x = vec + 2*(group + group_member);
				double *imag_x = vec + 2*(group + group_member) + 1;
				double *real_y = vec + 2*(group + group_member + stride);
				double *imag_y = vec + 2*(group + group_member + stride) + 1;
				double real_twiddled_y = *real_y * *twiddle_real_ptr - *imag_y * *twiddle_imag_ptr;
				double imag_twiddled_y = *imag_y * *twiddle_real_ptr + *real_y * *twiddle_imag_ptr;
				*real_y = *real_x - real_twiddled_y;
				*imag_y = *imag_x - imag_twiddled_y;
				*real_x = *real_x + real_twiddled_y;
				*imag_x = *imag_x + imag_twiddled_y;
			}
		}
	}
}

void invfft_line (
	pmhandle *pmdata,		// Handle for polymult library
	double	*vec,			// Array of complex values
	const unsigned int size)	// The fft size
{
	generate_sincos (pmdata, size);

	unsigned int stride = 1;		// First stride
	for ( ; stride < size; stride <<= 1) {
		double *twiddle_real_ptr = pmdata->twiddle_real;
		double *twiddle_imag_ptr = pmdata->twiddle_imag;
		for (unsigned int group = 0; group < size; group += stride * 2, twiddle_real_ptr++, twiddle_imag_ptr++) {
			for (unsigned int group_member = 0; group_member < stride; group_member++) {
				double *real_x = vec + 2*(group + group_member);
				double *imag_x = vec + 2*(group + group_member) + 1;
				double *real_y = vec + 2*(group + group_member + stride);
				double *imag_y = vec + 2*(group + group_member + stride) + 1;
				double xr = *real_x + *real_y;
				double xi = *imag_x + *imag_y;
				double yr = *real_x - *real_y;
				double yi = *imag_x - *imag_y;
				*real_x = xr;
				*imag_x = xi;
				*real_y = yr * *twiddle_real_ptr - yi * -*twiddle_imag_ptr;
				*imag_y = yr * -*twiddle_imag_ptr + yi * *twiddle_real_ptr;
			}
		}
	}
}

// out = in1 * in2 pointwise also divide by FFT length
void pointmul_line (double *invec1, double *invec2, double *outvec, int size)
{
	double inv_fft_size = 1.0 / (double) size;
	for (int j = 0; j < size; j++) {
		double real_x = invec1[2*j];
		double imag_x = invec1[2*j+1];
		double real_y = invec2[2*j];
		double imag_y = invec2[2*j+1];
		double tmp = real_x * real_y - imag_x * imag_y;
		outvec[2*j+1] = (real_x * imag_y + imag_x * real_y) * inv_fft_size;
		outvec[2*j] = tmp * inv_fft_size;
	}
}

// outvec += input vector shifted the proper amount to create a monic polymult result
void monic_line_add_hi (double *invec, int invec_size, bool invec_monic_rlp, double *outvec, int outvec_size)
{
	// Add words starting from the highest word
	for (invec += 2*(invec_size-1), outvec += 2*(outvec_size-1); ; invec -= 2, outvec -= 2, invec_size--) {
		if (invec_size == 0) {
			if (invec_monic_rlp) outvec[0] += 1.0;
			break;
		}
		outvec[0] += invec[0];
		outvec[1] += invec[1];
	}
}

// outvec += input vector shifted to the low part of outvec
void monic_line_add_lo (double *invec, int invec_size, double *outvec, int outvec_size)
{
	for ( ; invec_size; invec += 2, outvec += 2, invec_size--) {
		outvec[0] += invec[0];
		outvec[1] += invec[1];
	}
}

// Make adjustments to outvec for any monic input vectors
void monic_line_adjustments (double *invec1, int invec1_size, double *invec2, int invec2_size, double *outvec, int outvec_size, int options)
{
	// If neither input is monic, no adjustments are necessary
	if (! (options & POLYMULT_INVEC1_MONIC) && ! (options & POLYMULT_INVEC2_MONIC)) return;

	// If either input is a monic RLP, create a new entry before the output vector calculated thusfar
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) ||
	    ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))) {
		outvec -= 2;
		outvec_size++;
		memset (&outvec[0], 0, 2*sizeof(double));
	}

	// Create a new entry after the output vector calculated thusfar so that we can add to it
	outvec_size++;
	memset (&outvec[2*(outvec_size-1)], 0, 2*sizeof(double));

	// Adjust for an implied high coefficient of one -- add in the other poly in the high words
	if (options & POLYMULT_INVEC1_MONIC)
		monic_line_add_hi (invec2, invec2_size, (options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP), outvec, outvec_size);
	if (options & POLYMULT_INVEC2_MONIC)
		monic_line_add_hi (invec1, invec1_size, (options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP), outvec, outvec_size);

	// Adjust for an implied low coefficient of one -- add in the other poly to the low words
	if ((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP))
		monic_line_add_lo (invec2, invec2_size, outvec, outvec_size);
	if ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))
		monic_line_add_lo (invec1, invec1_size, outvec, outvec_size);
}

// Select from among polymult versions
void polymult_line (
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*gw_invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	gwnum	*gw_invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial
	gwnum	*gw_outvec,		// Output poly of size invec1_size + invec2_size - 1
	int	outvec_size,		// Size of the output polynomial.  Should be invec1_size + invec2_size (less one if not monic).
	int	options)
{
	int	must_fft, fft_size, adjusted_invec1_size, adjusted_invec2_size, adjusted_outvec_size;
	double	*invec1, *invec2, *outvec, *tmpvec;

	// Adjust sizes due to RLPs.  Calculate the output size assuming no monic inputs.
	adjusted_invec1_size = (options & POLYMULT_INVEC1_RLP) ? 2 * invec1_size - 1 : invec1_size;
	adjusted_invec2_size = (options & POLYMULT_INVEC2_RLP) ? 2 * invec2_size - 1 : invec2_size;
	adjusted_outvec_size = adjusted_invec1_size + adjusted_invec2_size - 1;

	// Polymult using no vector instructions.  Allocate memory.  For the output vector allocate two extra in case there are any monic RLP inputs.
	must_fft = (options & (POLYMULT_FORWARD_FFT | POLYMULT_INVEC1_FFTED | POLYMULT_INVEC2_FFTED));
	if (!must_fft && adjusted_outvec_size < pmdata->KARAT_BREAK) {
		invec1 = (double *) malloc (2 * adjusted_invec1_size * sizeof (double));
		invec2 = (double *) malloc (2 * adjusted_invec2_size * sizeof (double));
		outvec = (double *) malloc (2 * (adjusted_outvec_size + 2) * sizeof (double)) + 2;
		tmpvec = NULL;
	}
	else if (!must_fft && adjusted_outvec_size < pmdata->FFT_BREAK) {
		invec1 = (double *) malloc (2 * adjusted_invec1_size * sizeof (double));
		invec2 = (double *) malloc (2 * adjusted_invec2_size * sizeof (double));
		outvec = (double *) malloc (2 * (adjusted_outvec_size + 2) * sizeof (double)) + 2;
		tmpvec = (double *) malloc (2 * (adjusted_outvec_size + 32) * sizeof (double));		// Could need one extra for each recursion
	}
	else {
		fft_size = (options & POLYMULT_FORWARD_FFT) ? outvec_size :
			   (options & POLYMULT_INVEC1_FFTED) ? invec1_size :
			   (options & POLYMULT_INVEC2_FFTED) ? invec2_size : polymult_fft_size (adjusted_outvec_size);
		tmpvec = (double *) malloc (2 * fft_size * sizeof (double));		// FFT of invec1 goes here
		outvec = (double *) malloc (2 * (fft_size + 2) * sizeof (double)) + 2;	// FFT of invec2 goes here as does outvec
		if (options & POLYMULT_INVEC2_MONIC) invec1 = (double *) malloc (2 * fft_size * sizeof (double));
		else invec1 = tmpvec;
		if (options & POLYMULT_INVEC1_MONIC) invec2 = (double *) malloc (2 * fft_size * sizeof (double));
		else invec2 = outvec;
	}

	// Grab one complex value from each gwnum coefficient.  Work that data set.  Then move onto next complex val in gwnums.
	for ( ; ; ) {
		int line = next_line (pmdata);
		if (!read_line (pmdata->gwdata, gw_invec1, invec1_size, line, invec1, options & INVEC1_OPTIONS)) break;
		read_line (pmdata->gwdata, gw_invec2, invec2_size, line, invec2, options & INVEC2_OPTIONS);

		if (!must_fft && adjusted_outvec_size < pmdata->KARAT_BREAK) {
			brute_line (invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, outvec);
		}
		else if (!must_fft && adjusted_outvec_size < pmdata->FFT_BREAK) {
			karatsuba_line (pmdata, invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, tmpvec, outvec);
		}
		else {
			if (!(options & POLYMULT_INVEC1_FFTED)) {
				if (invec1 != tmpvec) memcpy (tmpvec, invec1, 2*adjusted_invec1_size*sizeof(double));
				memset (tmpvec + 2*adjusted_invec1_size, 0, 2*(fft_size-adjusted_invec1_size)*sizeof(double));
				fft_line (pmdata, tmpvec, fft_size);
			}
			if (options & POLYMULT_FORWARD_FFT) {
				write_line (pmdata->gwdata, gw_outvec, fft_size, line, tmpvec, options);
				continue;
			}
			if (!(options & POLYMULT_INVEC2_FFTED)) {
				if (invec2 != outvec) memcpy (outvec, invec2, 2*adjusted_invec2_size*sizeof(double));
				memset (outvec + 2*adjusted_invec2_size, 0, 2*(fft_size-adjusted_invec2_size)*sizeof(double));
				fft_line (pmdata, outvec, fft_size);
			}
			pointmul_line (tmpvec, outvec, outvec, fft_size);
			invfft_line (pmdata, outvec, fft_size);
		}
		monic_line_adjustments (invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, outvec, adjusted_outvec_size, options);
		write_line (pmdata->gwdata, gw_outvec, outvec_size, line, outvec, options);
	}
//for (int j=0; j < gwfftlen(pmdata->gwdata); j++){
//int off=addr_offset(pmdata->gwdata,j)/8;
//if (off % 8 >= 4) continue;	
//if (off == 0 || off == 2) ASSERTG(gw_outvec[0][off] == gw_invec1[0][off]* gw_invec2[0][off]);
//else if (off % 4 == 0) {
//	double r = gw_invec1[0][off] * gw_invec2[0][off] - gw_invec1[0][off+2] * gw_invec2[0][off+2];
//	double i = gw_invec1[0][off]* gw_invec2[0][off+2] + gw_invec1[0][off+2]* gw_invec2[0][off];
//ASSERTG(gw_outvec[0][off] == r);
//ASSERTG(gw_outvec[0][off+2] == i);
//}
//}
	// Free memory
	if (invec1 != tmpvec) free (invec1);
	if (invec2 != outvec) free (invec2);
	free (outvec - 2);
	free (tmpvec);
}

#endif

//--------------------------------------------------------------------
//		AVX optimized implementation
//--------------------------------------------------------------------

#ifdef AVX

bool read_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, __m256d *vec, int options)
{
	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += 2*(size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// If we've read all the gwnum data, return FALSE
	if (index * 8 >= (int) gwdata->FFTLEN) return (FALSE);

	// AVX implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(__m256d)); continue; }
			vec[2*j] = _mm256_set_pd (0.0, 0.0, 0.0, data[j][0]);
//bug? return 0.0 instead of -11
//			vec[2*j+1] = _mm256_set_pd (0.0, 0.0, 0.0, data[j][-11]);
			vec[2*j+1] = _mm256_set_pd (0.0, 0.0, 0.0, 0.0);
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(__m256d)); continue; }
			vec[2*j] = _mm256_set_pd (data[j][3], data[j][2], data[j][1], data[j][4]);
//bug? return 0.0 instead of -10
//			vec[2*j+1] = _mm256_set_pd (data[j][7], data[j][6], data[j][5], data[j][-10]);
			vec[2*j+1] = _mm256_set_pd (data[j][7], data[j][6], data[j][5], 0.0);
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
//BUG - handle various FFT paddings - addr_offset may not be right as the asm code is allowed to store data in the gaps.
// It might be easier/faster to just FFT the gap data.
// Or create new gwnum cache_line_offset that returns cache lines in order skipping the gaps
		unsigned long offset = cache_line_offset (gwdata, index) / sizeof (__m256d);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(__m256d)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			vec[2*j] = ((__m256d *) data[j])[offset];			// Real data
			vec[2*j+1] = ((__m256d *) data[j])[offset+1];			// Imaginary data
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		__m256d zero = _mm256_set_pd (0.0, 0.0, 0.0, 0.0);
		for (int j = 0; j < size; j++) {
			vec[2*j] = _mm256_sub_pd (zero, vec[2*j]);
			vec[2*j+1] = _mm256_sub_pd (zero, vec[2*j+1]);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) {
			vec[-2*j] = vec[2*j];
			vec[-2*j+1] = vec[2*j+1];
		}
	}

	// Return success
	return (TRUE);
}

void write_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, __m256d *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments_avx
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) || ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP)))
		vec -= 2;

	// If both inputs are RLPs, then the output is an RLP too.  Skip writing half of the data.
	if ((options & POLYMULT_INVEC1_RLP) && (options & POLYMULT_INVEC2_RLP)) vec += 2*(size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;
//BUG-Zero-padded not handled yet.  Need bigger gwnum header (if we allow storing ffted polys)
//BUG-current use of header location -10 conflicts with zero-padded use of location -10.

	// AVX implementation
	union {
		__m256d	a;
		double	b[4];
	} x;
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			x.a = vec[2*j], data[j][0] = x.b[0];
//bug? return 0.0 instead of -11
//			x.a = vec[2*j+1], data[j][-11] = x.b[0];
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			x.a = vec[2*j], data[j][4] = x.b[0], data[j][1] = x.b[1], data[j][2] = x.b[2], data[j][3] = x.b[3];
//bug? return 0.0 instead of -10
			x.a = vec[2*j+1], /*data[j][-10] = x.b[0],*/ data[j][5] = x.b[1], data[j][6] = x.b[2], data[j][7] = x.b[3];
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
//BUG - handle various FFT paddings - addr_offset may not be right as the asm code is allowed to store data in the gaps.
// It might be easier/faster to just FFT the gap data.
		unsigned long offset = cache_line_offset (gwdata, index) / sizeof (__m256d);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			((__m256d *) data[j])[offset] = vec[2*j];			// Real data
			((__m256d *) data[j])[offset+1] = vec[2*j+1];			// Imaginary data
		}
	}
}

// A brute force O(N^2) implementation of polynomial multiplication
void brute_line_avx (
	__m256d	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	__m256d	*invec2,		// Second input poly
	int	invec2_size,		// Size of the first input polynomial
	__m256d	*outvec)		// Output poly
{
	// AVX implementation
	for (int i = 0; i < invec1_size; i++) {
		__m256d real_x = invec1[2*i];
		__m256d imag_x = invec1[2*i+1];
		for (int j = 0; j < invec2_size; j++) {
			__m256d real_y = invec2[2*j];
			__m256d imag_y = invec2[2*j+1];
			__m256d real_tmp = _mm256_sub_pd (_mm256_mul_pd (real_x, real_y), _mm256_mul_pd (imag_x, imag_y));
			__m256d imag_tmp = _mm256_add_pd (_mm256_mul_pd (real_x, imag_y), _mm256_mul_pd (imag_x, real_y));
			if (i == 0 || j == invec2_size-1) {
				outvec[2*(i+j)] = real_tmp;
				outvec[2*(i+j)+1] = imag_tmp;
			} else {
				outvec[2*(i+j)] = _mm256_add_pd (outvec[2*(i+j)], real_tmp);
				outvec[2*(i+j)+1] = _mm256_add_pd (outvec[2*(i+j)+1], imag_tmp);
			}
		}
	}
}

void karatsuba_line_avx (
	pmhandle *pmdata,		// Handle for polymult library
	__m256d	*invec1,		// Coefficients of first input poly
	int	invec1_size,		// Size of the first input polynomial
	__m256d	*invec2,		// Coefficients of second input poly
	int	invec2_size,		// Size of the first input polynomial
	__m256d	*tmp,			// Temporary space
	__m256d	*outvec)		// Output poly coefficients
{
	int	i, a_size, b_size, c_size, d_size, a_plus_b_size, c_plus_d_size, bd_size, ac_size, a_plus_b_times_c_plus_d_size;
	__m256d	*a, *b, *c, *d, *ac, *bd, *a_plus_b, *c_plus_d, *a_plus_b_times_c_plus_d;

	// Handle end of recursion
	if (invec1_size == 1 || invec2_size == 1 || invec1_size + invec2_size < pmdata->KARAT_BREAK) {
		brute_line_avx (invec1, invec1_size, invec2, invec2_size, outvec);
		return;
	}

	// To simplify code, always make invec2 the larger poly.
	if (invec1_size > invec2_size) {
		int tmpsize = invec1_size; invec1_size = invec2_size; invec2_size = tmpsize;
		__m256d *tmpvec = invec1; invec1 = invec2; invec2 = tmpvec;
	}


	// Split larger vector roughly in half.  For readability, create pointers to two halves.
	d_size = (invec2_size + 1) / 2;
	c_size = invec2_size - d_size;
	c = invec2 + 2*d_size, d = invec2;

	// Handle case where vectors are of wildly unequal lengths
	if (d_size >= invec1_size) {
		// Compute bd
		bd = outvec;
		karatsuba_line_avx (pmdata, invec1, invec1_size, d, d_size, tmp, bd);
		bd_size = invec1_size + d_size - 1;
		// Compute bc and store it in tmp
		__m256d *bc = tmp;
		int bc_size = invec1_size + c_size - 1;
		tmp += 2*bc_size;
		karatsuba_line_avx (pmdata, invec1, invec1_size, c, c_size, tmp, bc);
		// add bc to output
		for (i = 0; i < 2*bc_size; i++) outvec[i + 2*d_size] = _mm256_add_pd (outvec[i + 2*d_size], bc[i]);
		return;
	}

	// Split smaller vector, the size of the two lower halves must be the same.  For readability, create pointers to two halves.
	b_size = d_size;
	a_size = invec1_size - b_size;
	a = invec1 + 2*b_size, b = invec1;

	// Compute (a+b)
	a_plus_b = outvec;
	for (i = 0; i < 2*a_size && i < 2*b_size; i++) a_plus_b[i] = _mm256_add_pd (a[i], b[i]);
//	for ( ; i < 2*a_size; i++) a_plus_b[i] = a[i];		// In case a_size > b_size
	for ( ; i < 2*b_size; i++) a_plus_b[i] = b[i];		// In case b_size > a_size
	a_plus_b_size = i / 2;
	// Compute (c+d)
	c_plus_d = &outvec[i];
	for (i = 0; i < 2*c_size && i < 2*d_size; i++) c_plus_d[i] = _mm256_add_pd (c[i], d[i]);
//	for ( ; i < 2*c_size; i++) c_plus_d[i] = c[i];		// In case c_size > d_size
	for ( ; i < 2*d_size; i++) c_plus_d[i] = d[i];		// In case d_size > c_size
	c_plus_d_size = i / 2;
	// Compute (a+b) * (c+d) and store it in tmp
	a_plus_b_times_c_plus_d = tmp;
	a_plus_b_times_c_plus_d_size = a_plus_b_size + c_plus_d_size - 1;
	tmp += 2*a_plus_b_times_c_plus_d_size;
	karatsuba_line_avx (pmdata, a_plus_b, a_plus_b_size, c_plus_d, c_plus_d_size, tmp, a_plus_b_times_c_plus_d);

	// Compute bd
	bd = outvec;
	karatsuba_line_avx (pmdata, b, b_size, d, d_size, tmp, bd);
	bd_size = b_size + d_size - 1;
	// Output the complex value 0+0i
	memset (&outvec[2*bd_size], 0, 2*sizeof(__m256d));
	// Compute ac
	ac = &outvec[2*(bd_size+1)];
	karatsuba_line_avx (pmdata, a, a_size, c, c_size, tmp, ac);
	ac_size = a_size + c_size - 1;

	// subtract ac from middle
	for (i = 0; i < 2*ac_size; i++) a_plus_b_times_c_plus_d[i] = _mm256_sub_pd (a_plus_b_times_c_plus_d[i], ac[i]);
	// subtract bd from middle
	for (i = 0; i < 2*bd_size; i++) a_plus_b_times_c_plus_d[i] = _mm256_sub_pd (a_plus_b_times_c_plus_d[i], bd[i]);
	// add middle to output
	for (i = 0; i < 2*a_plus_b_times_c_plus_d_size; i++) outvec[i + 2*b_size] = _mm256_add_pd (outvec[i + 2*b_size], a_plus_b_times_c_plus_d[i]);
}

void fft_line_avx (
	pmhandle *pmdata,		// Handle for polymult library
	__m256d	*vec,			// Array of complex values
	const unsigned int size)	// The fft size
{
	generate_sincos (pmdata, size);

	// Radix-4 FFT, one fewer complex multiply than radix-2, does more work while data is in registers
	unsigned int stride = size;
	while (stride >= 4) {
		stride >>= 2;				// Stride for this radix-4 step
		// First group can be done without twiddles
		for (unsigned int group_member = 0; group_member < stride; group_member++) {
			__m256d real_a = vec[2*(group_member)];
			__m256d imag_a = vec[2*(group_member) + 1];
			__m256d real_b = vec[2*(group_member + stride)];
			__m256d imag_b = vec[2*(group_member + stride) + 1];
			__m256d real_c = vec[2*(group_member + 2*stride)];
			__m256d imag_c = vec[2*(group_member + 2*stride) + 1];
			__m256d real_d = vec[2*(group_member + 3*stride)];
			__m256d imag_d = vec[2*(group_member + 3*stride) + 1];
			__m256d real_new_a = _mm256_add_pd (real_a, real_c);				// a + c
			__m256d imag_new_a = _mm256_add_pd (imag_a, imag_c);
			__m256d real_new_c = _mm256_sub_pd (real_a, real_c);				// a - c
			__m256d imag_new_c = _mm256_sub_pd (imag_a, imag_c);
			__m256d real_new_b = _mm256_add_pd (real_b, real_d);				// b + d
			__m256d imag_new_b = _mm256_add_pd (imag_b, imag_d);
			__m256d real_new_d = _mm256_sub_pd (real_b, real_d);				// b - d
			__m256d imag_new_d = _mm256_sub_pd (imag_b, imag_d);
			vec[2*(group_member)] = _mm256_add_pd (real_new_a, real_new_b);			// a + b
			vec[2*(group_member) + 1] = _mm256_add_pd (imag_new_a, imag_new_b);
			vec[2*(group_member + stride)] = _mm256_sub_pd (real_new_a, real_new_b);	// a - b
			vec[2*(group_member + stride) + 1] = _mm256_sub_pd (imag_new_a, imag_new_b);
			vec[2*(group_member + 2*stride)] = _mm256_sub_pd (real_new_c, imag_new_d);	// c + id
			vec[2*(group_member + 2*stride) + 1] = _mm256_add_pd (imag_new_c, real_new_d);
			vec[2*(group_member + 3*stride)] = _mm256_add_pd (real_new_c, imag_new_d);	// c - id
			vec[2*(group_member + 3*stride) + 1] = _mm256_sub_pd (imag_new_c, real_new_d);
		}
		for (unsigned int group = stride * 4, twiddle_step = 3; group < size; group += stride * 4, twiddle_step += 3) {
			__m256d real_twiddle1 = _mm256_broadcast_sd (pmdata->twiddle_real4 + twiddle_step);
			__m256d imag_twiddle1 = _mm256_broadcast_sd (pmdata->twiddle_imag4 + twiddle_step);
			__m256d real_twiddle2 = _mm256_broadcast_sd (pmdata->twiddle_real4 + twiddle_step + 1);
			__m256d imag_twiddle2 = _mm256_broadcast_sd (pmdata->twiddle_imag4 + twiddle_step + 1);
			__m256d real_twiddle3 = _mm256_broadcast_sd (pmdata->twiddle_real4 + twiddle_step + 2);
			__m256d imag_twiddle3 = _mm256_broadcast_sd (pmdata->twiddle_imag4 + twiddle_step + 2);
			for (unsigned int group_member = 0; group_member < stride; group_member++) {
				__m256d real_a = vec[2*(group + group_member)];
				__m256d imag_a = vec[2*(group + group_member) + 1];
				__m256d real_b = vec[2*(group + group_member + stride)];
				__m256d imag_b = vec[2*(group + group_member + stride) + 1];
				__m256d real_c = vec[2*(group + group_member + 2*stride)];
				__m256d imag_c = vec[2*(group + group_member + 2*stride) + 1];
				__m256d real_d = vec[2*(group + group_member + 3*stride)];
				__m256d imag_d = vec[2*(group + group_member + 3*stride) + 1];
				__m256d real_twiddled_b = _mm256_sub_pd (_mm256_mul_pd (real_b, real_twiddle1), _mm256_mul_pd (imag_b, imag_twiddle1));
				__m256d imag_twiddled_b = _mm256_add_pd (_mm256_mul_pd (imag_b, real_twiddle1), _mm256_mul_pd (real_b, imag_twiddle1));
				__m256d real_twiddled_c = _mm256_sub_pd (_mm256_mul_pd (real_c, real_twiddle2), _mm256_mul_pd (imag_c, imag_twiddle2));
				__m256d imag_twiddled_c = _mm256_add_pd (_mm256_mul_pd (imag_c, real_twiddle2), _mm256_mul_pd (real_c, imag_twiddle2));
				__m256d real_twiddled_d = _mm256_sub_pd (_mm256_mul_pd (real_d, real_twiddle3), _mm256_mul_pd (imag_d, imag_twiddle3));
				__m256d imag_twiddled_d = _mm256_add_pd (_mm256_mul_pd (imag_d, real_twiddle3), _mm256_mul_pd (real_d, imag_twiddle3));
				real_c = _mm256_sub_pd (real_a, real_twiddled_c);
				imag_c = _mm256_sub_pd (imag_a, imag_twiddled_c);
				real_a = _mm256_add_pd (real_a, real_twiddled_c);
				imag_a = _mm256_add_pd (imag_a, imag_twiddled_c);
				real_b = _mm256_add_pd (real_twiddled_b, real_twiddled_d);
				imag_b = _mm256_add_pd (imag_twiddled_b, imag_twiddled_d);
				real_d = _mm256_sub_pd (real_twiddled_b, real_twiddled_d);
				imag_d = _mm256_sub_pd (imag_twiddled_b, imag_twiddled_d);
				vec[2*(group + group_member)] = _mm256_add_pd (real_a, real_b);			// a + b
				vec[2*(group + group_member) + 1] = _mm256_add_pd (imag_a, imag_b);
				vec[2*(group + group_member + stride)] = _mm256_sub_pd (real_a, real_b);	// a - b
				vec[2*(group + group_member + stride) + 1] = _mm256_sub_pd (imag_a, imag_b);
				vec[2*(group + group_member + 2*stride)] = _mm256_sub_pd (real_c, imag_d);	// c + id
				vec[2*(group + group_member + 2*stride) + 1] = _mm256_add_pd (imag_c, real_d);
				vec[2*(group + group_member + 3*stride)] = _mm256_add_pd (real_c, imag_d);	// c - id
				vec[2*(group + group_member + 3*stride) + 1] = _mm256_sub_pd (imag_c, real_d);
			}
		}
	}

	// Radix-2 FFT
	while (stride >= 2) {
		stride >>= 1;				// Stride for this radix-2 step
		double	*twiddle_real_ptr = pmdata->twiddle_real;
		double	*twiddle_imag_ptr = pmdata->twiddle_imag;
		for (unsigned int group = 0; group < size; group += stride * 2, twiddle_real_ptr++, twiddle_imag_ptr++) {
			__m256d real_twiddle = _mm256_broadcast_sd (twiddle_real_ptr);
			__m256d imag_twiddle = _mm256_broadcast_sd (twiddle_imag_ptr);
			for (unsigned int group_member = 0; group_member < stride; group_member++) {
				__m256d *real_x = vec + 2*(group + group_member);
				__m256d *imag_x = vec + 2*(group + group_member) + 1;
				__m256d *real_y = vec + 2*(group + group_member + stride);
				__m256d *imag_y = vec + 2*(group + group_member + stride) + 1;
				__m256d real_twiddled_y = _mm256_sub_pd (_mm256_mul_pd (*real_y, real_twiddle), _mm256_mul_pd (*imag_y, imag_twiddle));
				__m256d imag_twiddled_y = _mm256_add_pd (_mm256_mul_pd (*imag_y, real_twiddle), _mm256_mul_pd (*real_y, imag_twiddle));
				*real_y = _mm256_sub_pd (*real_x, real_twiddled_y);
				*imag_y = _mm256_sub_pd (*imag_x, imag_twiddled_y);
				*real_x = _mm256_add_pd (*real_x, real_twiddled_y);
				*imag_x = _mm256_add_pd (*imag_x, imag_twiddled_y);
			}
		}
	}
}

void invfft_line_avx (
	pmhandle *pmdata,		// Handle for polymult library
	__m256d	*vec,			// Array of complex values
	const unsigned int size)	// The fft size
{
	generate_sincos (pmdata, size);

	// Radix-2 FFT
	unsigned int stride = 1;		// First stride
	// while (stride < size / 4) {		// Code necessary for a radix-2-only inverse FFT
	if (size & 0xAAAAAAAA) {		// If not a power of 4  (radix-2/4 inverse FFT)
		double	*twiddle_real_ptr = pmdata->twiddle_real;
		double	*twiddle_imag_ptr = pmdata->twiddle_imag;
		for (unsigned int group = 0; group < size; group += stride * 2, twiddle_real_ptr++, twiddle_imag_ptr++) {
			double	imag_twiddle_negated = -*twiddle_imag_ptr;
			__m256d real_twiddle = _mm256_broadcast_sd (twiddle_real_ptr);
			__m256d imag_twiddle = _mm256_broadcast_sd (&imag_twiddle_negated);
			for (unsigned int group_member = 0; group_member < stride; group_member++) {
				__m256d *real_x = vec + 2*(group + group_member);
				__m256d *imag_x = vec + 2*(group + group_member) + 1;
				__m256d *real_y = vec + 2*(group + group_member + stride);
				__m256d *imag_y = vec + 2*(group + group_member + stride) + 1;
				__m256d xr = _mm256_add_pd (*real_x, *real_y);
				__m256d xi = _mm256_add_pd (*imag_x, *imag_y);
				__m256d yr = _mm256_sub_pd (*real_x, *real_y);
				__m256d yi = _mm256_sub_pd (*imag_x, *imag_y);
				*real_x = xr;
				*imag_x = xi;
				*real_y = _mm256_sub_pd (_mm256_mul_pd (yr, real_twiddle), _mm256_mul_pd (yi, imag_twiddle));
				*imag_y = _mm256_add_pd (_mm256_mul_pd (yr, imag_twiddle), _mm256_mul_pd (yi, real_twiddle));
			}
		}
		stride <<= 1;			// Next stride amount
	}

	// Radix-4 FFT, one fewer complex multiply than radix-2, does more work while data is in registers
	for ( ; stride < size; stride <<= 2) {
		for (unsigned int group = 0, twiddle_step = 0; group < size; group += stride * 4, twiddle_step += 3) {
			__m256d real_twiddle1 = _mm256_broadcast_sd (pmdata->twiddle_real4 + twiddle_step);
			__m256d imag_twiddle1 = _mm256_broadcast_sd (pmdata->twiddle_imag4 + twiddle_step);
			__m256d real_twiddle2 = _mm256_broadcast_sd (pmdata->twiddle_real4 + twiddle_step + 1);
			__m256d imag_twiddle2 = _mm256_broadcast_sd (pmdata->twiddle_imag4 + twiddle_step + 1);
			__m256d real_twiddle3 = _mm256_broadcast_sd (pmdata->twiddle_real4 + twiddle_step + 2);
			__m256d imag_twiddle3 = _mm256_broadcast_sd (pmdata->twiddle_imag4 + twiddle_step + 2);
			for (unsigned int group_member = 0; group_member < stride; group_member++) {
				__m256d real_a = vec[2*(group + group_member)];
				__m256d imag_a = vec[2*(group + group_member) + 1];
				__m256d real_b = vec[2*(group + group_member + stride)];
				__m256d imag_b = vec[2*(group + group_member + stride) + 1];
				__m256d real_c = vec[2*(group + group_member + 2*stride)];
				__m256d imag_c = vec[2*(group + group_member + 2*stride) + 1];
				__m256d real_d = vec[2*(group + group_member + 3*stride)];
				__m256d imag_d = vec[2*(group + group_member + 3*stride) + 1];

				__m256d next_real_a = _mm256_add_pd (real_a, real_b);			// a + b
				__m256d next_imag_a = _mm256_add_pd (imag_a, imag_b);
				__m256d next_real_b = _mm256_sub_pd (real_a, real_b);			// a - b
				__m256d next_imag_b = _mm256_sub_pd (imag_a, imag_b);
				__m256d next_real_c = _mm256_add_pd (real_c, real_d);			// c + d
				__m256d next_imag_c = _mm256_add_pd (imag_c, imag_d);
				__m256d next_real_d = _mm256_sub_pd (real_c, real_d);			// c - d
				__m256d next_imag_d = _mm256_sub_pd (imag_c, imag_d);

				real_a = _mm256_add_pd (next_real_a, next_real_c);			// a + c
				imag_a = _mm256_add_pd (next_imag_a, next_imag_c);
				real_c = _mm256_sub_pd (next_real_a, next_real_c);			// a - c
				imag_c = _mm256_sub_pd (next_imag_a, next_imag_c);
				real_b = _mm256_add_pd (next_real_b, next_imag_d);			// b - id
				imag_b = _mm256_sub_pd (next_imag_b, next_real_d);
				real_d = _mm256_sub_pd (next_real_b, next_imag_d);			// b + id
				imag_d = _mm256_add_pd (next_imag_b, next_real_d);

				// Special case group zero, no twiddles necessary
				if (group == 0) {
					vec[2*(group_member)] = real_a;
					vec[2*(group_member) + 1] = imag_a;
					vec[2*(group_member + stride)] = real_b;
					vec[2*(group_member + stride) + 1] = imag_b;
					vec[2*(group_member + 2*stride)] = real_c;
					vec[2*(group_member + 2*stride) + 1] = imag_c;
					vec[2*(group_member + 3*stride)] = real_d;
					vec[2*(group_member + 3*stride) + 1] = imag_d;
					continue;
				}

				// Apply the CONJUGATE of the twiddle
				vec[2*(group + group_member)] = real_a;
				vec[2*(group + group_member) + 1] = imag_a;
				vec[2*(group + group_member + stride)] =
					_mm256_add_pd (_mm256_mul_pd (real_b, real_twiddle1), _mm256_mul_pd (imag_b, imag_twiddle1));
				vec[2*(group + group_member + stride) + 1] =
					_mm256_sub_pd (_mm256_mul_pd (imag_b, real_twiddle1), _mm256_mul_pd (real_b, imag_twiddle1));
				vec[2*(group + group_member + 2*stride)] =
					_mm256_add_pd (_mm256_mul_pd (real_c, real_twiddle2), _mm256_mul_pd (imag_c, imag_twiddle2));
				vec[2*(group + group_member + 2*stride) + 1] =
					_mm256_sub_pd (_mm256_mul_pd (imag_c, real_twiddle2), _mm256_mul_pd (real_c, imag_twiddle2));
				vec[2*(group + group_member + 3*stride)] =
					_mm256_add_pd (_mm256_mul_pd (real_d, real_twiddle3), _mm256_mul_pd (imag_d, imag_twiddle3));
				vec[2*(group + group_member + 3*stride) + 1] =
					_mm256_sub_pd (_mm256_mul_pd (imag_d, real_twiddle3), _mm256_mul_pd (real_d, imag_twiddle3));
			}
		}
	}
}

// out = in1 * in2 pointwise also divide by FFT length
void pointmul_line_avx (__m256d *invec1, __m256d *invec2, __m256d *outvec, int size)
{
	double	inv_fft_size = 1.0 / (double) size;
	__m256d inv_fft_size256 = _mm256_broadcast_sd (&inv_fft_size);
	for (int j = 0; j < size; j++) {
		__m256d real_x = invec1[2*j];
		__m256d imag_x = invec1[2*j+1];
		__m256d real_y = invec2[2*j];
		__m256d imag_y = invec2[2*j+1];
		outvec[2*j+1] = _mm256_mul_pd (_mm256_add_pd (_mm256_mul_pd (real_x, imag_y), _mm256_mul_pd (imag_x, real_y)), inv_fft_size256);
		outvec[2*j] = _mm256_mul_pd (_mm256_sub_pd (_mm256_mul_pd (real_x, real_y), _mm256_mul_pd (imag_x, imag_y)), inv_fft_size256);
	}
}

// outvec += input vector shifted the proper amount to create a monic polymult result
void monic_line_add_hi_avx (__m256d *invec, int invec_size, bool invec_monic_rlp, __m256d *outvec, int outvec_size)
{
	// Add words starting from the highest word
	for (invec += 2*(invec_size-1), outvec += 2*(outvec_size-1); ; invec -= 2, outvec -= 2, invec_size--) {
		if (invec_size == 0) {
			if (invec_monic_rlp) outvec[0] = _mm256_add_pd (outvec[0], _mm256_set_pd (1.0, 1.0, 1.0, 1.0));
			break;
		}
		outvec[0] = _mm256_add_pd (outvec[0], invec[0]);
		outvec[1] = _mm256_add_pd (outvec[1], invec[1]);
	}
}

// outvec += input vector shifted to the low part of outvec
void monic_line_add_lo_avx (__m256d *invec, int invec_size, __m256d *outvec, int outvec_size)
{
	for ( ; invec_size; invec += 2, outvec += 2, invec_size--) {
		outvec[0] = _mm256_add_pd (outvec[0], invec[0]);
		outvec[1] = _mm256_add_pd (outvec[1], invec[1]);
	}
}

// Make adjustments to outvec for any monic input vectors
void monic_line_adjustments_avx (__m256d *invec1, int invec1_size, __m256d *invec2, int invec2_size, __m256d *outvec, int outvec_size, int options)
{
	// If neither input is monic, no adjustments are necessary
	if (! (options & POLYMULT_INVEC1_MONIC) && ! (options & POLYMULT_INVEC2_MONIC)) return;

	// If either input is a monic RLP, create a new entry before the output vector calculated thusfar
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) ||
	    ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))) {
		outvec -= 2;
		outvec_size++;
		memset (&outvec[0], 0, 2*sizeof(__m256d));
	}

	// Create a new entry after the output vector calculated thusfar so that we can add to it
	outvec_size++;
	memset (&outvec[2*(outvec_size-1)], 0, 2*sizeof(__m256d));

	// Adjust for an implied high coefficient of one -- add in the other poly in the high words
	if (options & POLYMULT_INVEC1_MONIC)
		monic_line_add_hi_avx (invec2, invec2_size, (options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP), outvec, outvec_size);
	if (options & POLYMULT_INVEC2_MONIC)
		monic_line_add_hi_avx (invec1, invec1_size, (options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP), outvec, outvec_size);

	// Adjust for an implied low coefficient of one -- add in the other poly to the low words
	if ((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP))
		monic_line_add_lo_avx (invec2, invec2_size, outvec, outvec_size);
	if ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))
		monic_line_add_lo_avx (invec1, invec1_size, outvec, outvec_size);
}

// Select from various versions of polymult algorithms
void polymult_line_avx (
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*gw_invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	gwnum	*gw_invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial
	gwnum	*gw_outvec,		// Output poly of size invec1_size + invec2_size - 1
	int	outvec_size,		// Size of the output polynomial.  Should be invec1_size + invec2_size (less one if not monic).
	int	options)
{
	int	must_fft, fft_size, adjusted_invec1_size, adjusted_invec2_size, adjusted_outvec_size;
	__m256d	*invec1, *invec2, *outvec, *tmpvec;

	// Adjust sizes due to RLPs.  Calculate the output size assuming no monic inputs.
	adjusted_invec1_size = (options & POLYMULT_INVEC1_RLP) ? 2 * invec1_size - 1 : invec1_size;
	adjusted_invec2_size = (options & POLYMULT_INVEC2_RLP) ? 2 * invec2_size - 1 : invec2_size;
	adjusted_outvec_size = adjusted_invec1_size + adjusted_invec2_size - 1;

	// Polymult using AVX instructions.  Allocate memory.  For the output vector allocate two extra in case there are any monic RLP inputs.
	must_fft = (options & (POLYMULT_FORWARD_FFT | POLYMULT_INVEC1_FFTED | POLYMULT_INVEC2_FFTED));
	if (!must_fft && adjusted_outvec_size < pmdata->KARAT_BREAK) {
		invec1 = (__m256d *) aligned_malloc (2 * adjusted_invec1_size * sizeof (__m256d), 64);
		invec2 = (__m256d *) aligned_malloc (2 * adjusted_invec2_size * sizeof (__m256d), 64);
		outvec = (__m256d *) aligned_malloc (2 * (adjusted_outvec_size + 2) * sizeof (__m256d), 64) + 2;
		tmpvec = NULL;
	}
	else if (!must_fft && adjusted_outvec_size < pmdata->FFT_BREAK) {
		invec1 = (__m256d *) aligned_malloc (2 * adjusted_invec1_size * sizeof (__m256d), 64);
		invec2 = (__m256d *) aligned_malloc (2 * adjusted_invec2_size * sizeof (__m256d), 64);
		outvec = (__m256d *) aligned_malloc (2 * (adjusted_outvec_size + 2) * sizeof (__m256d), 64) + 2;
		tmpvec = (__m256d *) aligned_malloc (2 * (adjusted_outvec_size + 32) * sizeof (__m256d), 64);	// Could need one extra for each recursion
	}
	else {
		fft_size = (options & POLYMULT_FORWARD_FFT) ? outvec_size :
			   (options & POLYMULT_INVEC1_FFTED) ? invec1_size :
			   (options & POLYMULT_INVEC2_FFTED) ? invec2_size : polymult_fft_size (adjusted_outvec_size);
		tmpvec = (__m256d *) aligned_malloc (2 * fft_size * sizeof (__m256d), 64);		// FFT of invec1 goes here
		outvec = (__m256d *) aligned_malloc (2 * (fft_size + 2) * sizeof (__m256d), 64) + 2;	// FFT of invec2 goes here as does outvec
		if (options & POLYMULT_INVEC2_MONIC) invec1 = (__m256d *) aligned_malloc (2 * fft_size * sizeof (__m256d), 64);
		else invec1 = tmpvec;
		if (options & POLYMULT_INVEC1_MONIC) invec2 = (__m256d *) aligned_malloc (2 * fft_size * sizeof (__m256d), 64);
		else invec2 = outvec;
	}

	// Grab four complex values from each gwnum coefficient.  Work that data set.  Then move onto remaining complex vals in gwnums.
	for ( ; ; ) {
		int line = next_line (pmdata);
		if (!read_line_avx (pmdata->gwdata, gw_invec1, invec1_size, line, invec1, options & INVEC1_OPTIONS)) break;
		read_line_avx (pmdata->gwdata, gw_invec2, invec2_size, line, invec2, options & INVEC2_OPTIONS);

		if (!must_fft && adjusted_outvec_size < pmdata->KARAT_BREAK) {
			brute_line_avx (invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, outvec);
		}
		else if (!must_fft && adjusted_outvec_size < pmdata->FFT_BREAK) {
			karatsuba_line_avx (pmdata, invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, tmpvec, outvec);
		}
		else {
			if (!(options & POLYMULT_INVEC1_FFTED)) {
				if (invec1 != tmpvec) memcpy (tmpvec, invec1, 2*adjusted_invec1_size*sizeof(__m256d));
				memset (tmpvec + 2*adjusted_invec1_size, 0, 2*(fft_size-adjusted_invec1_size)*sizeof(__m256d));
				fft_line_avx (pmdata, tmpvec, fft_size);
			}
			if (options & POLYMULT_FORWARD_FFT) {
				write_line_avx (pmdata->gwdata, gw_outvec, fft_size, line, tmpvec, options);
				continue;
			}
			if (!(options & POLYMULT_INVEC2_FFTED)) {
				if (invec2 != outvec) memcpy (outvec, invec2, 2*adjusted_invec2_size*sizeof(__m256d));
				memset (outvec + 2*adjusted_invec2_size, 0, 2*(fft_size-adjusted_invec2_size)*sizeof(__m256d));
				fft_line_avx (pmdata, outvec, fft_size);
			}
			pointmul_line_avx (tmpvec, outvec, outvec, fft_size);
			invfft_line_avx (pmdata, outvec, fft_size);
		}
		monic_line_adjustments_avx (invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, outvec, adjusted_outvec_size, options);
		write_line_avx (pmdata->gwdata, gw_outvec, outvec_size, line, outvec, options);
	}

	// Free memory
	if (invec1 != tmpvec) aligned_free (invec1);
	if (invec2 != outvec) aligned_free (invec2);
	aligned_free (outvec - 2);
	aligned_free (tmpvec);
}

#endif



#if !defined (AVX) && !defined (FMA) && !defined (AVX512)

struct polymult_thread_info {
	pmhandle *pmdata;		// Handle for polymult library
	gwnum	*invec1;		// First input poly
	int	invec1_size;		// Size of the first input polynomial (or fft size if POLYMULT_INVEC1_FFTED)
	gwnum	*invec2;		// Second input poly
	int	invec2_size;		// Size of the second input polynomial (or fft size if POLYMULT_INVEC2_FFTED)
	gwnum	*outvec;		// Output poly
	int	outvec_size;		// Size of the output polynomial (or fft_size if POLYMULT_FORWARD_FFT)
	int	options;
};

/* Entry point for one thread to do part of the polymult */
void polymult_thread (void *arg)
{
	struct polymult_thread_info *info = (struct polymult_thread_info *) arg;

	// Optimized routines for AVX instructions
	if (info->pmdata->gwdata->cpu_flags & CPU_AVX)
		polymult_line_avx (info->pmdata, info->invec1, info->invec1_size, info->invec2, info->invec2_size, info->outvec, info->outvec_size, info->options);
	// Default routines using no vector instructions
	else
		polymult_line (info->pmdata, info->invec1, info->invec1_size, info->invec2, info->invec2_size, info->outvec, info->outvec_size, info->options);
}

/* Multiply two polynomials producing FFTed gwnums as output.  The output gwnums need to be normalized before use in further gwnum operations. */
/* It is safe to use input gwnums in the output vector.  For a normal polynomial multiply outvec_size must be invec1_size + invec2_size - 1. */
/* For monic polynomial multiply, the leading input coefficients of 1 are omitted as is the leading 1 output coefficient -- thus */
/* outvec_size must be invec1_size + invec2_size. */
void polymult (
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial (or fft size if POLYMULT_INVEC1_FFTED)
	gwnum	*invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial (or fft size if POLYMULT_INVEC2_FFTED)
	gwnum	*outvec,		// Output poly
	int	outvec_size,		// Size of the output polynomial (or fft_size if POLYMULT_FORWARD_FFT)
	int	options)
{
	// Forward FFT not allowed on already FFTed polynomial
	ASSERTG (! ((options & POLYMULT_FORWARD_FFT) && (options & POLYMULT_INVEC1_FFTED)));
	// An already FFTed input cannot be multiplied with a monic polynomial
	ASSERTG (! ((options & POLYMULT_INVEC1_FFTED) && !(options & POLYMULT_INVEC2_MONIC)));
	ASSERTG (! ((options & POLYMULT_INVEC2_FFTED) && !(options & POLYMULT_INVEC1_MONIC)));
	// An already FFTed input must be at least as big as the output vector
	ASSERTG (! ((options & POLYMULT_INVEC1_FFTED) && invec1_size >= outvec_size));
	ASSERTG (! ((options & POLYMULT_INVEC2_FFTED) && invec2_size >= outvec_size));
	// Verify output vector size when we know the input vector sizes (no pre-FFTed inputs)
#define insz1	(invec1_size + !!(options & POLYMULT_INVEC1_MONIC))		// Adjust size for monic poly
#define insz2	(invec2_size + !!(options & POLYMULT_INVEC2_MONIC))
#define insz1a	((options & POLYMULT_INVEC1_RLP) ? 2*insz1-1 : insz1)		// Adjust size for RLP poly
#define insz2a	((options & POLYMULT_INVEC2_RLP) ? 2*insz2-1 : insz2)
#define outsz	(insz1a + insz2a - 1)						// Output poly size
#define outsz1	((options & POLYMULT_INVEC1_RLP && options & POLYMULT_INVEC2_RLP) ? (outsz + 1) / 2 : outsz)	// Output size adjusted for RLP
#define outsz2	(outsz1 - !!(options & POLYMULT_INVEC1_MONIC && options & POLYMULT_INVEC2_MONIC))		// Output size adjusted for monic
	ASSERTG ((options & (POLYMULT_FORWARD_FFT | POLYMULT_INVEC1_FFTED | POLYMULT_INVEC2_FFTED) || outvec_size == outsz2));
#undef insz1
#undef insz2
#undef insz1a
#undef insz2a
#undef outsz
#undef outsz1
#undef outsz2

	// FFT inputs if needed
	if (!(options & POLYMULT_INVEC1_FFTED)) for (int j = 0; j < invec1_size; j++) if (invec1[j] != NULL) gwfft (pmdata->gwdata, invec1[j], invec1[j]);
	if (!(options & POLYMULT_INVEC2_FFTED)) for (int j = 0; j < invec2_size; j++) if (invec2[j] != NULL) gwfft (pmdata->gwdata, invec2[j], invec2[j]);

	// Prepare for polymult in parallel
	pmdata->next_line = 0;
	struct polymult_thread_info info;
	info.pmdata = pmdata;
	info.invec1 = invec1;
	info.invec1_size = invec1_size;
	info.invec2 = invec2;
	info.invec2_size = invec2_size;
	info.outvec = outvec;
	info.outvec_size = outvec_size;
	info.options = options;

	// Fire up each additional thread, then use the current thread too
//GW: bug - fixed array limits num_threads
	gwthread thread_id[512];
	for (int i = 1; i < pmdata->num_threads; i++) gwthread_create_waitable (&thread_id[i], &polymult_thread, (void *) &info);
	polymult_thread (&info);

	// Wait for all threads to complete
	for (int i = 1; i < pmdata->num_threads; i++) gwthread_wait_for_exit (&thread_id[i]);

	// Unfft all output coefficients
	if (!(options & POLYMULT_FORWARD_FFT)) {
		for (int j = 0; j < outvec_size; j++) {
			if (outvec[j] == NULL) continue;
			FFT_state (outvec[j]) = FULLY_FFTed;
			gwunfft2 (pmdata->gwdata, outvec[j], outvec[j], (options & POLYMULT_STARTNEXTFFT) ? GWMUL_STARTNEXTFFT : 0);
		}
	} else {
		// Mark all results poly FFTed
//		for (int j = 0; j < outvec_size; j++)
//GW - set unnorm count too? Don't need to do that if always unffting.
//GW - have a special FFT state for POLYMULT_FORWARD_FFT?
//		ASSERTG (outvec[j] != NULL);
//		FFT_state (outvec[j]) = FULLY_FFTed;
	}
}

#endif

