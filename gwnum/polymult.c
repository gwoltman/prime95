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
#if defined (SSE2) || defined (AVX) || defined (FMA) || defined (AVX512)
#include <immintrin.h>
#endif

// TO DO:
// 1) Test AVX-512, FMA, SSE2 versions
// 2) Zero padded gwnums.  Test all FFT sizes and FFT types -- especially funky PFA ffts.
// 7) Variable GW_HEADER_SIZE.
// 9) Is a fast convert-gwnum-to-larger-FFT feasible? necessary?
// 10) Tune for Brute / Karatsuba / FFT crossovers.
// 12) Tune routine that returns needed safety margin needed for given poly sizes
// 14) Compressed gwnums.  Does it pay off?
// 15) Add polymult to gwtest or as a standalone program
// 16) Handle out-of-memory gracefully
// 17) Sin/cos accuracy, right collection of PFA FFT sizes?, 2-pass FFTs
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
// 43) Stockham auto-sort, read_line/write_line that reads/writes data with stride=1, two-pass cache-friendly ffts
// 44) Cache several sincos tables for different FFT sizes.
// 45) If twiddle_stride gets large, build a smaller twiddle table for smaller strides.
// 46) Look for more FMA opportunities.
// 47) Combine last fft step with pointmul and first invfft step?
// 48) Test if hyperthreading is beneficial
// 49) Problems if FFT(1)!=1!!! The pointwise mul code ought to divide by FFT(1) OR have a special asm entry point for unfft that deals with data that
//	has already had a pointwise mul.  One test case is AVX-512 small FFTs.
//     Furthermore, the monic_adjustments code can't simply add data - it must add FFT(1)*data.


// Forward declarations required for the SSE2, AVX, FMA, AVX512 compiles
int next_line (pmhandle *pmdata);
unsigned long cache_line_offset (gwhandle *gwdata, unsigned long i);
void generate_sincos (pmhandle *pmdata, const unsigned int size);
void polymult_line_sse2 (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, gwnum *gw_invec2, int invec2_size, gwnum *gw_outvec, int outvec_size, int options);
void polymult_line_avx (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, gwnum *gw_invec2, int invec2_size, gwnum *gw_outvec, int outvec_size, int options);
void polymult_line_fma (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, gwnum *gw_invec2, int invec2_size, gwnum *gw_outvec, int outvec_size, int options);
void polymult_line_avx512 (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, gwnum *gw_invec2, int invec2_size, gwnum *gw_outvec, int outvec_size, int options);

#define INVEC1_OPTIONS	(POLYMULT_INVEC1_MONIC | POLYMULT_INVEC1_RLP | POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC1_FFTED)
#define INVEC2_OPTIONS	(POLYMULT_INVEC2_MONIC | POLYMULT_INVEC2_RLP | POLYMULT_INVEC2_NEGATE | POLYMULT_INVEC2_FFTED)

#if defined (AVX512)
#define VLEN		512			// Vector is 256 bits wide
#define VDT		__m512d			// Vector data type
#define VFMA		1			// FMA is supported
#define	read_line	read_line_avx512	// Naming extension for this option
#define	write_line	write_line_avx512
#define	brute_line	brute_line_avx512
#define	karatsuba_line	karatsuba_line_avx512
#define	fft_line	fft_line_avx512
#define	invfft_line	invfft_line_avx512
#define	pointmul_line	pointmul_line_avx512
#define	monic_line_add_hi monic_line_add_hi_avx512
#define	monic_line_add_lo monic_line_add_lo_avx512
#define	monic_line_adjustments monic_line_adjustments_avx512
#define	polymult_line	polymult_line_avx512
#define broadcastsd	_mm512_set1_pd
#define addpd		_mm512_add_pd
#define subpd		_mm512_sub_pd
#define mulpd		_mm512_mul_pd
#define fmaddpd		_mm512_fmadd_pd
#define fmsubpd		_mm512_fmsub_pd
#define fnmaddpd	_mm512_fnmadd_pd
#elif defined (FMA)
#define VLEN		256			// Vector is 256 bits wide
#define VDT		__m256d			// Vector data type
#define VFMA		1			// FMA is supported
#define	read_line	read_line_fma		// Naming extension for this option
#define	write_line	write_line_fma
#define	brute_line	brute_line_fma
#define	karatsuba_line	karatsuba_line_fma
#define	fft_line	fft_line_fma
#define	invfft_line	invfft_line_fma
#define	pointmul_line	pointmul_line_fma
#define	monic_line_add_hi monic_line_add_hi_fma
#define	monic_line_add_lo monic_line_add_lo_fma
#define	monic_line_adjustments monic_line_adjustments_fma
#define	polymult_line	polymult_line_fma
#define broadcastsd	_mm256_set1_pd
#define addpd		_mm256_add_pd
#define subpd		_mm256_sub_pd
#define mulpd		_mm256_mul_pd
#define fmaddpd		_mm256_fmadd_pd
#define fmsubpd		_mm256_fmsub_pd
#define fnmaddpd	_mm256_fnmadd_pd
#elif defined (AVX)
#define VLEN		256			// Vector is 256 bits wide
#define VDT		__m256d			// Vector data type
#define VFMA		0			// FMA not supported
#define	read_line	read_line_avx		// Naming extension for this option
#define	write_line	write_line_avx
#define	brute_line	brute_line_avx
#define	karatsuba_line	karatsuba_line_avx
#define	fft_line	fft_line_avx
#define	invfft_line	invfft_line_avx
#define	pointmul_line	pointmul_line_avx
#define	monic_line_add_hi monic_line_add_hi_avx
#define	monic_line_add_lo monic_line_add_lo_avx
#define	monic_line_adjustments monic_line_adjustments_avx
#define	polymult_line	polymult_line_avx
#define broadcastsd	_mm256_set1_pd
#define addpd		_mm256_add_pd
#define subpd		_mm256_sub_pd
#define mulpd		_mm256_mul_pd
#elif defined (SSE2)
#define VLEN		128			// Vector is 256 bits wide
#define VDT		__m128d			// Vector data type
#define VFMA		0			// FMA not supported
#define	read_line	read_line_sse2		// Naming extension for this option
#define	write_line	write_line_sse2
#define	brute_line	brute_line_sse2
#define	karatsuba_line	karatsuba_line_sse2
#define	fft_line	fft_line_sse2
#define	invfft_line	invfft_line_sse2
#define	pointmul_line	pointmul_line_sse2
#define	monic_line_add_hi monic_line_add_hi_sse2
#define	monic_line_add_lo monic_line_add_lo_sse2
#define	monic_line_adjustments monic_line_adjustments_sse2
#define	polymult_line	polymult_line_sse2
#define broadcastsd	_mm_set1_pd
#define addpd		_mm_add_pd
#define subpd		_mm_sub_pd
#define mulpd		_mm_mul_pd
#else
#define VLEN		64			// Vector is 256 bits wide
#define VDT		double			// Vector data type
#define VFMA		0			// FMA not supported
#define	read_line	read_line_dbl		// Naming extension for this option
#define	write_line	write_line_dbl
#define	brute_line	brute_line_dbl
#define	karatsuba_line	karatsuba_line_dbl
#define	fft_line	fft_line_dbl
#define	invfft_line	invfft_line_dbl
#define	pointmul_line	pointmul_line_dbl
#define	monic_line_add_hi monic_line_add_hi_dbl
#define	monic_line_add_lo monic_line_add_lo_dbl
#define	monic_line_adjustments monic_line_adjustments_dbl
#define	polymult_line	polymult_line_dbl
#define broadcastsd	_mm64_set1_pd
#define addpd		_mm64_add_pd
#define subpd		_mm64_sub_pd
#define mulpd		_mm64_mul_pd
// Emulate vector length 1 using Intel-style instrinsic naming conventions
#define __m64d			double
#define _mm64_set1_pd(a)	(a)
#define _mm64_add_pd(a,b)	((a)+(b))
#define _mm64_sub_pd(a,b)	((a)-(b))
#define _mm64_mul_pd(a,b)	((a)*(b))
#endif

// Emulate FMA on CPUs that do not have an FMA instruction
#if !VFMA
#define fmaddpd(a,b,c)	addpd(mulpd(a,b),c)
#define fmsubpd(a,b,c)	subpd(mulpd(a,b),c)
#define fnmaddpd(a,b,c)	subpd(c,mulpd(a,b))
#endif

#if !defined (SSE2) && !defined (AVX) && !defined (FMA) && !defined (AVX512)

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
//GW: These FFT sizes haven't been timed.  Some may be slower than a larger FFT or more radix-3/5 options might be profitable
	int	pfas[] = {40, 45, 48, 50, 54, 60, 64, 72, 75};		// Provide a variety of radix-3 and radix-5 FFT sizes

	// Power of two when n is small -- brute force is probably faster anyway
	if (n <= 32) for (int i = 2; ; i *= 2) if (i >= n) return (i);

	// Increase n by a power of two until it
	for (int two_multiplier = 1; ; two_multiplier *= 2) {
		if (n > two_multiplier * 75) continue;
		for (int i = 0; i < sizeof (pfas) / sizeof (int); i++) if (two_multiplier * pfas[i] >= n) return (two_multiplier * pfas[i]);
	}
}

// Initialize a polymult handle
void polymult_init (
	pmhandle *pmdata,		// Handle for polymult library that needs initializing
	gwhandle *gwdata)		// Handle for gwnum FFT library
{
	memset (pmdata, 0, sizeof (pmhandle));
	pmdata->gwdata = gwdata;

	// Initialize default number of threads
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
	if (pmdata->thread_ids != NULL) {
		pmdata->termination_in_progress = TRUE;
		gwevent_signal (&pmdata->work_to_do);
		for (int i = 1; i < pmdata->num_threads; i++) gwthread_wait_for_exit (&pmdata->thread_ids[i]);
		free (pmdata->thread_ids);
		pmdata->thread_ids = NULL;
		gwmutex_destroy (&pmdata->poly_mutex);
		gwevent_destroy (&pmdata->work_to_do);
		gwevent_destroy (&pmdata->helpers_done);
	}
	aligned_free (pmdata->twiddles1), pmdata->twiddles1 = NULL;
	aligned_free (pmdata->twiddles2), pmdata->twiddles2 = NULL;
}


// Return the next line that needs processing
int next_line (
	pmhandle *pmdata)		// Handle for polymult library
{
	int	return_value;

	// Use mutex to guard next_line counter, return next available line
	if (pmdata->num_threads > 1) gwmutex_lock (&pmdata->poly_mutex);
	return_value = pmdata->next_line++;
	if (pmdata->num_threads > 1) gwmutex_unlock (&pmdata->poly_mutex);
	return (return_value);
}

//GW -- move to gwnum.c?

/* Get offset to nth "cache line" from an FFTed gwnum.  Where a cache line is 64 bytes. */

unsigned long cache_line_offset (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long i)
{
	unsigned long addr, grps;

/* Memory layouts for AVX-512 machines */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		if (gwdata->PASS1_SIZE == 0) {
			addr = i * 64;
			/* Now optionally add 64 pad bytes every 4KB */
			if (gwdata->FOURKBGAPSIZE) addr = addr + (addr / (gwdata->FOURKBGAPSIZE << 6)) * 64;
		}
		else {
			grps = i / (gwdata->PASS2_SIZE / 4);
			addr = i * 64;
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

// Generate radix-3/4/5 twiddle factors
void generate_sincos (
	pmhandle *pmdata,		// Handle for polymult library
	const unsigned int fft_size)
{
	const double twopi = 2.0 * 3.1415926535897932384626433;
	unsigned int size = fft_size;

	// Return quickly if twiddles have already been initialized
	if (pmdata->twiddles_initialized && pmdata->twiddles_initialized == size) return;

	// Make sure only one thread initializes the twiddles
	if (pmdata->num_threads > 1) gwmutex_lock (&pmdata->poly_mutex);

	// If some other thread beat us to the initialization, then return
	if (pmdata->twiddles_initialized && pmdata->twiddles_initialized == size) {
		if (pmdata->num_threads > 1) gwmutex_unlock (&pmdata->poly_mutex);
		return;
	}

	// Free the previous twiddles
	if (pmdata->twiddles_initialized) {
//GW: accept size < pmdata->twiddles_initialized (at least for power of 2)
//GW: cache several different table sizes?
		aligned_free (pmdata->twiddles1), pmdata->twiddles1 = NULL;
		aligned_free (pmdata->twiddles2), pmdata->twiddles2 = NULL;
	}

//GW: Consider long doubles for better accuracy
	if (size % 3 == 0) {
		pmdata->twiddles1 = (double *) aligned_malloc (size/3 * 2 * sizeof (double), 64);
		for (unsigned int i = 0; i < size/3; i++) {
			double angle = twopi * (double) i / (double) size;
			pmdata->twiddles1[2*i] = cos (angle);
			pmdata->twiddles1[2*i+1] = sin (angle);
		}
		while (size % 3 == 0) size /= 3;
	}
	pmdata->twiddles2 = (double *) aligned_malloc (size/4 * 2 * 2 * sizeof (double), 64);
	for (unsigned int i = 0; i < size/4; i++) {
		double angle = twopi * (double) i / (double) size;
		pmdata->twiddles2[2*(2*i)] = cos (angle);
		pmdata->twiddles2[2*(2*i)+1] = sin (angle);
		pmdata->twiddles2[2*(2*i+1)] = cos (2.0*angle);
		pmdata->twiddles2[2*(2*i+1)+1] = sin (2.0*angle);
	}

	pmdata->twiddles_initialized = fft_size;
	if (pmdata->num_threads > 1) gwmutex_unlock (&pmdata->poly_mutex);
}

#endif

//--------------------------------------------------------------------
//		Read/write cache lines from/to a gwnum
//    Requires intimate knowledge of layout of gwnums in memory
//--------------------------------------------------------------------

#if defined (AVX512)

bool read_line_avx512 (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options, bool post_process_monics)
{
	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
//GW: check if GW_FFT1 should be used instead
		vec[0] = broadcastsd (1.0);
		vec[1] = broadcastsd (0.0);
		vec += 2;
	}

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
	if (index * 16 >= (int) gwdata->FFTLEN) return (FALSE);

	// AVX512 implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			vec[2*j] = _mm512_set_pd (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, data[j][0]);
//bug? return 0.0 instead of -11
//			vec[2*j+1] = _mm512_set_pd (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, data[j][-11]);
			vec[2*j+1] = _mm512_set_pd (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			vec[2*j] = _mm512_set_pd (data[j][7], data[j][6], data[j][5], data[j][4], data[j][3], data[j][2], data[j][1], data[j][8]);
//bug? return 0.0 instead of -10
//			vec[2*j+1] = _mm512_set_pd (data[j][7], data[j][6], data[j][5], data[j][-10]);
			vec[2*j+1] = _mm512_set_pd (data[j][15], data[j][14], data[j][13], data[j][12], data[j][11], data[j][10], data[j][9], 0.0);
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index * 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			vec[2*j] = ((VDT *) data[j])[offset];				// Real data
			vec[2*j+1] = ((VDT *) data[j])[offset+1];			// Imaginary data
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		VDT zero = broadcastsd (0.0);
		for (int j = 0; j < size; j++) {
			vec[2*j] = subpd (zero, vec[2*j]);
			vec[2*j+1] = subpd (zero, vec[2*j+1]);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) {
			vec[-2*j] = vec[2*j];
			vec[-2*j+1] = vec[2*j+1];
		}
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
//GW: check if GW_FFT1 should be used instead
		vec[2*size] = broadcastsd (1.0);
		vec[2*size+1] = broadcastsd (0.0);
	}

	// Return success
	return (TRUE);
}

void write_line_avx512 (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments
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

	// AVX512 implementation
	union {
		VDT	a;
		double	b[8];
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
			x.a = vec[2*j], data[j][7] = x.b[7], data[j][6] = x.b[6], data[j][5] = x.b[5], data[j][4] = x.b[4],
					data[j][3] = x.b[3], data[j][2] = x.b[2], data[j][1] = x.b[1], data[j][8] = x.b[0];
//bug? return 0.0 instead of -10
			x.a = vec[2*j+1], data[j][15] = x.b[7], data[j][14] = x.b[6], data[j][13] = x.b[5], data[j][12] = x.b[4],
					  data[j][11] = x.b[3], data[j][10] = x.b[2], data[j][9] = x.b[1];
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index * 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			((VDT *) data[j])[offset] = vec[2*j];			// Real data
			((VDT *) data[j])[offset+1] = vec[2*j+1];		// Imaginary data
		}
	}
}

#elif defined (FMA)				// FMA and AVX FFTs use the same memory layout

#define read_line_fma	read_line_avx
#define write_line_fma	write_line_avx
bool read_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options, bool post_process_monics);
void write_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options);

#elif defined (AVX)

bool read_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options, bool post_process_monics)
{
	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
//GW: check if GW_FFT1 should be used instead
		vec[0] = broadcastsd (1.0);
		vec[1] = broadcastsd (0.0);
		vec += 2;
	}

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
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			vec[2*j] = _mm256_set_pd (0.0, 0.0, 0.0, data[j][0]);
//bug? return 0.0 instead of -11
//			vec[2*j+1] = _mm256_set_pd (0.0, 0.0, 0.0, data[j][-11]);
			vec[2*j+1] = _mm256_set_pd (0.0, 0.0, 0.0, 0.0);
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			vec[2*j] = _mm256_set_pd (data[j][3], data[j][2], data[j][1], data[j][4]);
//bug? return 0.0 instead of -10
//			vec[2*j+1] = _mm256_set_pd (data[j][7], data[j][6], data[j][5], data[j][-10]);
			vec[2*j+1] = _mm256_set_pd (data[j][7], data[j][6], data[j][5], 0.0);
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			vec[2*j] = ((VDT *) data[j])[offset];				// Real data
			vec[2*j+1] = ((VDT *) data[j])[offset+1];			// Imaginary data
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		VDT zero = broadcastsd (0.0);
		for (int j = 0; j < size; j++) {
			vec[2*j] = subpd (zero, vec[2*j]);
			vec[2*j+1] = subpd (zero, vec[2*j+1]);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) {
			vec[-2*j] = vec[2*j];
			vec[-2*j+1] = vec[2*j+1];
		}
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
//GW: check if GW_FFT1 should be used instead
		vec[2*size] = broadcastsd (1.0);
		vec[2*size+1] = broadcastsd (0.0);
	}

	// Return success
	return (TRUE);
}

void write_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments
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
		VDT	a;
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
			x.a = vec[2*j], data[j][3] = x.b[3], data[j][2] = x.b[2], data[j][1] = x.b[1], data[j][4] = x.b[0];
//bug? return 0.0 instead of -10
			x.a = vec[2*j+1], data[j][7] = x.b[3], data[j][6] = x.b[2], data[j][5] = x.b[1] /*, data[j][-10] = x.b[0]*/ ;
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			((VDT *) data[j])[offset] = vec[2*j];			// Real data
			((VDT *) data[j])[offset+1] = vec[2*j+1];		// Imaginary data
		}
	}
}

#elif defined (SSE2)

bool read_line_sse2 (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options, bool post_process_monics)
{
	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
//GW: check if GW_FFT1 should be used instead
		vec[0] = broadcastsd (1.0);
		vec[1] = broadcastsd (0.0);
		vec += 2;
	}

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
	if (index * 4 >= (int) gwdata->FFTLEN) return (FALSE);

	// Various SSE2 FFT implementations store the imaginary part of a complex number either 8, 16, or 32 bytes above the real part
	int	imag_dist;
	if (gwdata->FFT_TYPE > 0 &&					// r4, r4delay, r4dwpn with 8-complex final levels
	    (gwdata->PASS2_SIZE == 1536 || gwdata->PASS2_SIZE == 2048 || gwdata->PASS2_SIZE == 2560 ||
	     gwdata->PASS2_SIZE == 4608 || gwdata->PASS2_SIZE == 6144 || gwdata->PASS2_SIZE == 7680 ||
	     gwdata->PASS2_SIZE == 8192 || gwdata->PASS2_SIZE == 10240 || gwdata->PASS2_SIZE == 12800)) imag_dist = 32;
	else if (gwdata->ALL_COMPLEX_FFT) imag_dist = 16;
	else if (gwdata->PASS2_SIZE == 0 && index < 4) imag_dist = 8;
	else if (gwdata->PASS2_SIZE != 0 && gwdata->FFT_TYPE == 0 && (index <= 1 || index == 4 || index == 5)) imag_dist = 8;	// hg ffts
	else imag_dist = 16;

	// SSE2 instructions implementation
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			vec[2*j] = _mm_set_pd (0.0, data[j][0]);
//bug? return 0.0 instead of -11 ---- fix when we support FFTing invec
			vec[2*j+1] = _mm_set_pd (0.0, 0.0);
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			if (imag_dist == 8) {
				vec[2*j] = _mm_set_pd (data[j][2], data[j][1]);
				vec[2*j+1] = _mm_set_pd (data[j][3], 0.0); //data[j][-10];
			} else if (imag_dist == 16) {
				vec[2*j] = _mm_set_pd (data[j][1], data[j][2]);
//bug? return 0.0 instead of -10.  How are we going to implement FFTed poly read?  Its own routine? It's own FFT state? Compressed?
				vec[2*j+1] = _mm_set_pd (data[j][3], 0.0); //data[j][-10];
			} else {
				vec[2*j] = _mm_set_pd (data[j][1], data[j][4]);
				vec[2*j+1] = _mm_set_pd (data[j][5], 0.0); //data[j][-10];
			}
		}
	}
	else {				// Process complex values where imaginary part is 8, 16, or 32 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			if (imag_dist == 8) {
				vec[2*j] = _mm_set_pd (data[j][index*4+2], data[j][index*4]);
				vec[2*j+1] = _mm_set_pd (data[j][index*4+3], data[j][index*4+1]);
			} else if (imag_dist == 16) {
				vec[2*j] = ((VDT *) data[j])[offset+(index&1)*2];			// Real data
				vec[2*j+1] = ((VDT *) data[j])[offset+(index&1)*2+1];			// Imaginary data
			} else {
				vec[2*j] = ((VDT *) data[j])[offset+(index&1)];				// Real data
				vec[2*j+1] = ((VDT *) data[j])[offset+(index&1)+2];			// Imaginary data
			}
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		VDT zero = broadcastsd (0.0);
		for (int j = 0; j < size; j++) {
			vec[2*j] = subpd (zero, vec[2*j]);
			vec[2*j+1] = subpd (zero, vec[2*j+1]);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) {
			vec[-2*j] = vec[2*j];
			vec[-2*j+1] = vec[2*j+1];
		}
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
//GW: check if GW_FFT1 should be used instead
		vec[2*size] = broadcastsd (1.0);
		vec[2*size+1] = broadcastsd (0.0);
	}

	// Return success
	return (TRUE);
}

void write_line_sse2 (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments
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

	// Various SSE2 FFT implementations store the imaginary part of a complex number either 8, 16, or 32 bytes above the real part
	int	imag_dist;
	if (gwdata->FFT_TYPE > 0 &&					// r4, r4delay, r4dwpn with 8-complex final levels
	    (gwdata->PASS2_SIZE == 1536 || gwdata->PASS2_SIZE == 2048 || gwdata->PASS2_SIZE == 2560 ||
	     gwdata->PASS2_SIZE == 4608 || gwdata->PASS2_SIZE == 6144 || gwdata->PASS2_SIZE == 7680 ||
	     gwdata->PASS2_SIZE == 8192 || gwdata->PASS2_SIZE == 10240 || gwdata->PASS2_SIZE == 12800)) imag_dist = 32;
	else if (gwdata->ALL_COMPLEX_FFT) imag_dist = 16;
	else if (gwdata->PASS2_SIZE == 0 && index < 4) imag_dist = 8;
	else if (gwdata->PASS2_SIZE != 0 && gwdata->FFT_TYPE == 0 && (index <= 1 || index == 4 || index == 5)) imag_dist = 8;	// hg ffts
	else imag_dist = 16;

	// SSE2 instructions implementation
	union {
		VDT	a;
		double	b[2];
	} x;
	if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			x.a = vec[2*j], data[j][0] = x.b[0];
//			data[j][-11] = vec[2*j+1];  //GW:bug?
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {		// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			;
			if (imag_dist == 8) {
				x.a = vec[2*j], data[j][2] = x.b[1], data[j][1] = x.b[0];
				x.a = vec[2*j+1], data[j][3] = x.b[1];
			} else if (imag_dist == 16) {
				x.a = vec[2*j], data[j][1] = x.b[1], data[j][2] = x.b[0];
				x.a = vec[2*j+1], data[j][3] = x.b[1];
			} else {
				x.a = vec[2*j], data[j][1] = x.b[1], data[j][4] = x.b[0];
				x.a = vec[2*j+1], data[j][5] = x.b[1];
			}
//bug? return 0.0 instead of -10
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (imag_dist == 8) {
				x.a = vec[2*j], data[j][index*4+2] = x.b[1], data[j][index*4] = x.b[0];
				x.a = vec[2*j+1], data[j][index*4+3] = x.b[1], data[j][index*4+1] = x.b[0];
			} else if (imag_dist == 16) {
				((VDT *) data[j])[offset+(index&1)*2] = vec[2*j];			// Real data
				((VDT *) data[j])[offset+(index&1)*2+1] = vec[2*j+1];			// Imaginary data
			} else {
				((VDT *) data[j])[offset+(index&1)] = vec[2*j];				// Real data
				((VDT *) data[j])[offset+(index&1)+2] = vec[2*j+1];			// Imaginary data
			}
		}
	}
}

#elif !defined (X86_64)				// x87 FFTs only supported in 32-bit mode

bool read_line_dbl (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options, bool post_process_monics)
{
	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
//GW: check if GW_FFT1 should be used instead
		vec[0] = 1.0;
		vec[1] = 0.0;
		vec += 2;
	}

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
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			vec[2*j] = data[j][0];
//bug? return 0.0 instead of -11 ---- fix when we support FFTing invec
			vec[2*j+1] = 0.0; //data[j][-11];
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
			if (gwdata->PASS2_SIZE == 0) vec[2*j] = data[j][1];
			else vec[2*j] = data[j][2];
//bug? return 0.0 instead of -10.  How are we going to implement FFTed poly read?  Its own routine? It's own FFT state? Compressed?
			vec[2*j+1] = 0.0; //data[j][-10];
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 4) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[2*j], 0, 2*sizeof(VDT)); continue; }
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

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
//GW: check if GW_FFT1 should be used instead
		vec[2*size] = 1.0;
		vec[2*size+1] = 0.0;
	}

	// Return success
	return (TRUE);
}

void write_line_dbl (gwhandle *gwdata, gwnum *data, int size, int index, VDT *vec, int options)
{
	// If either input was a monic RLP, back up the output vector pointer to the newly created entry by monic_line_adjustments
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
		unsigned long offset = cache_line_offset (gwdata, index / 4) / sizeof (VDT);
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

#endif


/*--------------------------------------------------------------------------
|	     Generalized routines (vector data type independent)
+-------------------------------------------------------------------------*/

#if defined (AVX512) || defined (FMA) || defined (AVX) || defined (SSE2) || !defined (X86_64)

void brute_line (
	VDT	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	VDT	*invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial
	VDT	*outvec,		// Output poly
	int	outvec_size,		// Either invec1_size+invec2_size-1 or size of gw_outvec if circular convolution
	int	options,		// Options passed to polymult
	gwnum	*gw_outvec)		// Gwnum array passed to polymult (so we can detect outputs that do not need to be computed)
{
	// AVX implementation.  Loop through output indices.
	for (int out = 0; out < invec1_size + invec2_size - 1; out++) {
		// Adjust output index for circular convolutions
		int adjusted_out = out;
		bool first_time = TRUE;
		if (out >= outvec_size) {
			adjusted_out = out - outvec_size;
			first_time = FALSE;
		}
		// See if this output needs to be computed
		if ((options & POLYMULT_CIRCULAR) && gw_outvec[adjusted_out] == NULL) continue;			//GW: Catch more cases such as RLP * RLP
		// Get index into input vectors
		int in1 = out < invec2_size ? 0 : out - invec2_size + 1;	// Index into invec1
		int in2 = out - in1;						// Index into invec2
		// Loop through all the pairs to accumulate in this output
		for ( ; in1 < invec1_size && in2 >= 0; in1++, in2--) {
			VDT real_x = invec1[2*in1];
			VDT imag_x = invec1[2*in1+1];
			VDT real_y = invec2[2*in2];
			VDT imag_y = invec2[2*in2+1];
			VDT real_tmp = fmsubpd (real_x, real_y, mulpd (imag_x, imag_y));
			VDT imag_tmp = fmaddpd (real_x, imag_y, mulpd (imag_x, real_y));
			if (first_time) {
				outvec[2*adjusted_out] = real_tmp;
				outvec[2*adjusted_out+1] = imag_tmp;
				first_time = FALSE;
			} else {
				outvec[2*adjusted_out] = addpd (outvec[2*adjusted_out], real_tmp);
				outvec[2*adjusted_out+1] = addpd (outvec[2*adjusted_out+1], imag_tmp);
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
	VDT	*invec1,		// Coefficients of first input poly
	int	invec1_size,		// Size of the first input polynomial
	VDT	*invec2,		// Coefficients of second input poly
	int	invec2_size,		// Size of the second input polynomial
	VDT	*tmp,			// Temporary space
	VDT	*outvec)		// Output poly coefficients
{
	int	i, a_size, b_size, c_size, d_size, a_plus_b_size, c_plus_d_size, bd_size, ac_size, a_plus_b_times_c_plus_d_size;
	VDT	*a, *b, *c, *d, *ac, *bd, *a_plus_b, *c_plus_d, *a_plus_b_times_c_plus_d;

	// Handle end of recursion
	if (invec1_size == 1 || invec2_size == 1 || invec1_size + invec2_size < pmdata->KARAT_BREAK) {
		brute_line (invec1, invec1_size, invec2, invec2_size, outvec, invec1_size + invec2_size - 1, 0, NULL);
		return;
	}

	// To simplify code, always make invec2 the larger poly.
	if (invec1_size > invec2_size) {
		int tmpsize = invec1_size; invec1_size = invec2_size; invec2_size = tmpsize;
		VDT *tmpvec = invec1; invec1 = invec2; invec2 = tmpvec;
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
		VDT *bc = tmp;
		int bc_size = invec1_size + c_size - 1;
		tmp += 2*bc_size;
		karatsuba_line (pmdata, invec1, invec1_size, c, c_size, tmp, bc);
		// add/copy bc to output
		for (i = 0; i < 2*(invec1_size-1); i++) outvec[i + 2*d_size] = addpd (outvec[i + 2*d_size], bc[i]);
		for ( ; i < 2*bc_size; i++) outvec[i + 2*d_size] = bc[i];
		return;
	}

	// Split smaller vector, the size of the two lower halves must be the same.  For readability, create pointers to two halves.
	b_size = d_size;
	a_size = invec1_size - b_size;
	a = invec1 + 2*b_size, b = invec1;

	// Compute (a+b)
	a_plus_b = outvec;
	for (i = 0; i < 2*a_size && i < 2*b_size; i++) a_plus_b[i] = addpd (a[i], b[i]);
//	for ( ; i < 2*a_size; i++) a_plus_b[i] = a[i];		// In case a_size > b_size
	for ( ; i < 2*b_size; i++) a_plus_b[i] = b[i];		// In case b_size > a_size
	a_plus_b_size = i / 2;
	// Compute (c+d)
	c_plus_d = &outvec[i];
	for (i = 0; i < 2*c_size && i < 2*d_size; i++) c_plus_d[i] = addpd (c[i], d[i]);
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
	memset (&outvec[2*bd_size], 0, 2*sizeof(VDT));
	// Compute ac
	ac = &outvec[2*(bd_size+1)];
	karatsuba_line (pmdata, a, a_size, c, c_size, tmp, ac);
	ac_size = a_size + c_size - 1;

	// subtract ac from middle
	for (i = 0; i < 2*ac_size; i++) a_plus_b_times_c_plus_d[i] = subpd (a_plus_b_times_c_plus_d[i], ac[i]);
	// subtract bd from middle
	for (i = 0; i < 2*bd_size; i++) a_plus_b_times_c_plus_d[i] = subpd (a_plus_b_times_c_plus_d[i], bd[i]);
	// add middle to output
	for (i = 0; i < 2*a_plus_b_times_c_plus_d_size; i++) outvec[i + 2*b_size] = addpd (outvec[i + 2*b_size], a_plus_b_times_c_plus_d[i]);
}

void fft_line (
	pmhandle *pmdata,		// Handle for polymult library
	VDT	*vec,			// Array of complex values
	unsigned int fft_size)		// The fft size
{
	unsigned int size = fft_size;		// The current group size to process
	unsigned int stride;			// Stride through FFT data
	unsigned int twiddle_stride;		// Stride through the twiddles data

	generate_sincos (pmdata, fft_size);

	// Radix-3 FFT
	// A 3-complex FFT is:
	// Res1:  (R1+R2+R3) + (I1+I2+I3)i
	// Res2:  (R1-.5R2-.5R3-.866I2+.866I3) + (I1-.5I2-.5I3+.866R2-.866R3)i
	// Res3:  (R1-.5R2-.5R3+.866I2-.866I3) + (I1-.5I2-.5I3-.866R2+.866R3)i
	twiddle_stride = 1;
	for ( ; size % 3 == 0; size /= 3, twiddle_stride *= 3) {
		stride = size / 3;
		VDT half = broadcastsd (0.5);
		VDT _866 = broadcastsd (0.86602540378443864676372317075294);
		for (unsigned int group = 0; group < fft_size; group += stride * 3) {
			for (unsigned int group_member = 0, twiddle_idx = 0; group_member < stride; group_member++, twiddle_idx += twiddle_stride) {
				VDT real_a = vec[2*(group + group_member)];
				VDT imag_a = vec[2*(group + group_member)+1];
				VDT real_b = vec[2*(group + group_member + stride)];
				VDT imag_b = vec[2*(group + group_member + stride)+1];
				VDT real_c = vec[2*(group + group_member + 2*stride)];
				VDT imag_c = vec[2*(group + group_member + 2*stride)+1];

				VDT real_b_plus_c = addpd (real_b, real_c);				// b + c
				VDT imag_b_plus_c = addpd (imag_b, imag_c);

				VDT real_tmp23a = fnmaddpd (half, real_b_plus_c, real_a);		// a - 0.5 * (b + c)
				VDT imag_tmp23a = fnmaddpd (half, imag_b_plus_c, imag_a);
//GW: delay mul by 866 on true FMA CPUs
				VDT real_tmp23b = mulpd (_866, subpd (imag_b, imag_c));	// 0.866 * (b - c)
				VDT imag_tmp23b = mulpd (_866, subpd (real_b, real_c));

				real_a = addpd (real_a, real_b_plus_c);					// a + b + c
				imag_a = addpd (imag_a, imag_b_plus_c);
				real_b = subpd (real_tmp23a, real_tmp23b);				// combine 23a and 23b
				imag_b = addpd (imag_tmp23a, imag_tmp23b);
				real_c = addpd (real_tmp23a, real_tmp23b);
				imag_c = subpd (imag_tmp23a, imag_tmp23b);

				vec[2*(group + group_member)]   = real_a;
				vec[2*(group + group_member)+1] = imag_a;
				if (group_member == 0) {	// First group member's twiddles are 1+0i, skip twiddle multiply
					vec[2*(group + group_member + stride)]     = real_b;
					vec[2*(group + group_member + stride)+1]   = imag_b;
					vec[2*(group + group_member + 2*stride)]   = real_c;
					vec[2*(group + group_member + 2*stride)+1] = imag_c;
				}
				else {			// Apply twiddles
					VDT real_twiddle1 = broadcastsd (pmdata->twiddles1[2*twiddle_idx]);
					VDT imag_twiddle1 = broadcastsd (pmdata->twiddles1[2*twiddle_idx+1]);
					vec[2*(group + group_member + stride)]     = fmsubpd (real_b, real_twiddle1, mulpd (imag_b, imag_twiddle1));
					vec[2*(group + group_member + stride)+1]   = fmaddpd (imag_b, real_twiddle1, mulpd (real_b, imag_twiddle1));
					vec[2*(group + group_member + 2*stride)]   = fmaddpd (real_c, real_twiddle1, mulpd (imag_c, imag_twiddle1));
					vec[2*(group + group_member + 2*stride)+1] = fmsubpd (imag_c, real_twiddle1, mulpd (real_c, imag_twiddle1));
				}
			}
		}
	}

	// Radix-5 FFT
	// R1= r1     +(r2+r5)     +(r3+r4)
	// R2= r1 +.309(r2+r5) -.809(r3+r4)    -.951(i2-i5) -.588(i3-i4)
	// R5= r1 +.309(r2+r5) -.809(r3+r4)    +.951(i2-i5) +.588(i3-i4)
	// R3= r1 -.809(r2+r5) +.309(r3+r4)    -.588(i2-i5) +.951(i3-i4)
	// R4= r1 -.809(r2+r5) +.309(r3+r4)    +.588(i2-i5) -.951(i3-i4)
	// I1= i1     +(i2+i5)     +(i3+i4)
	// I2= i1 +.309(i2+i5) -.809(i3+i4)    +.951(r2-r5) +.588(r3-r4)
	// I5= i1 +.309(i2+i5) -.809(i3+i4)    -.951(r2-r5) -.588(r3-r4)
	// I3= i1 -.809(i2+i5) +.309(i3+i4)    +.588(r2-r5) -.951(r3-r4)
	// I4= i1 -.809(i2+i5) +.309(i3+i4)    -.588(r2-r5) +.951(r3-r4)
	twiddle_stride = 1;
	for ( ; size % 5 == 0; size /= 5, twiddle_stride *= 5) {
		VDT _309 = broadcastsd (0.30901699437494742410229341718282);
		VDT _809 = broadcastsd (0.80901699437494742410229341718282);
		VDT _951 = broadcastsd (0.95105651629515357211643933337938);
		VDT _588_951 = broadcastsd (0.61803398874989484820458683436564);
		stride = size / 5;
		for (unsigned int group = 0; group < fft_size; group += stride * 5) {
			for (unsigned int group_member = 0, twiddle_idx = 0; group_member < stride; group_member++, twiddle_idx += 2*twiddle_stride) {
				VDT real_a = vec[2*(group + group_member)];
				VDT imag_a = vec[2*(group + group_member)+1];
				VDT real_b = vec[2*(group + group_member + stride)];
				VDT imag_b = vec[2*(group + group_member + stride)+1];
				VDT real_c = vec[2*(group + group_member + 2*stride)];
				VDT imag_c = vec[2*(group + group_member + 2*stride)+1];
				VDT real_d = vec[2*(group + group_member + 3*stride)];
				VDT imag_d = vec[2*(group + group_member + 3*stride)+1];
				VDT real_e = vec[2*(group + group_member + 4*stride)];
				VDT imag_e = vec[2*(group + group_member + 4*stride)+1];

				VDT real_b_plus_e = addpd (real_b, real_e);				// b + e
				VDT imag_b_plus_e = addpd (imag_b, imag_e);
				VDT real_c_plus_d = addpd (real_c, real_d);				// c + d
				VDT imag_c_plus_d = addpd (imag_c, imag_d);
				VDT real_b_minus_e = subpd (real_b, real_e);				// b - e
				VDT imag_b_minus_e = subpd (imag_b, imag_e);
				VDT real_c_minus_d = subpd (real_c, real_d);				// c - d
				VDT imag_c_minus_d = subpd (imag_c, imag_d);

				VDT real_tmp25a = fmaddpd (_309, real_b_plus_e, fnmaddpd (_809, real_c_plus_d, real_a)); // a +.309*(b+e) -.809*(c+d)
				VDT imag_tmp25a = fmaddpd (_309, imag_b_plus_e, fnmaddpd (_809, imag_c_plus_d, imag_a));
				VDT real_tmp34a = fnmaddpd (_809, real_b_plus_e, fmaddpd (_309, real_c_plus_d, real_a)); // a -.809*(b+e) +.309*(c+d)
				VDT imag_tmp34a = fnmaddpd (_809, imag_b_plus_e, fmaddpd (_309, imag_c_plus_d, imag_a));
//GW: delay mul by .951 on true FMA CPUs
				VDT real_tmp25b = mulpd (_951, fmaddpd (_588_951, imag_c_minus_d, imag_b_minus_e));	// .951*(b-e) + 0.588*(c-d)
				VDT imag_tmp25b = mulpd (_951, fmaddpd (_588_951, real_c_minus_d, real_b_minus_e));
				VDT real_tmp34b = mulpd (_951, fmsubpd (_588_951, imag_b_minus_e, imag_c_minus_d));	// .588*(b-e) - 0.951*(c-d)
				VDT imag_tmp34b = mulpd (_951, fmsubpd (_588_951, real_b_minus_e, real_c_minus_d));

				real_a = addpd (addpd (real_a, real_b_plus_e), real_c_plus_d);	// a + b + c + d + e
				imag_a = addpd (addpd (imag_a, imag_b_plus_e), imag_c_plus_d);
				real_b = subpd (real_tmp25a, real_tmp25b);				// combine 25a and 25b
				imag_b = addpd (imag_tmp25a, imag_tmp25b);
				real_e = addpd (real_tmp25a, real_tmp25b);
				imag_e = subpd (imag_tmp25a, imag_tmp25b);
				real_c = subpd (real_tmp34a, real_tmp34b);				// combine 34a and 34b
				imag_c = addpd (imag_tmp34a, imag_tmp34b);
				real_d = addpd (real_tmp34a, real_tmp34b);
				imag_d = subpd (imag_tmp34a, imag_tmp34b);

				vec[2*(group + group_member)]   = real_a;
				vec[2*(group + group_member)+1] = imag_a;
				if (group_member == 0) {	// First group member's twiddles are 1+0i, skip twiddle multiply
					vec[2*(group + group_member + stride)]     = real_b;
					vec[2*(group + group_member + stride)+1]   = imag_b;
					vec[2*(group + group_member + 2*stride)]   = real_c;
					vec[2*(group + group_member + 2*stride)+1] = imag_c;
					vec[2*(group + group_member + 3*stride)]   = real_d;
					vec[2*(group + group_member + 3*stride)+1] = imag_d;
					vec[2*(group + group_member + 4*stride)]   = real_e;
					vec[2*(group + group_member + 4*stride)+1] = imag_e;
				}
				else {			// Apply twiddles
					VDT real_twiddle1 = broadcastsd (pmdata->twiddles2[2*twiddle_idx]);
					VDT imag_twiddle1 = broadcastsd (pmdata->twiddles2[2*twiddle_idx+1]);
					VDT real_twiddle2 = broadcastsd (pmdata->twiddles2[2*(twiddle_idx+1)]);
					VDT imag_twiddle2 = broadcastsd (pmdata->twiddles2[2*(twiddle_idx+1)+1]);
					vec[2*(group + group_member + stride)]     = fmsubpd (real_b, real_twiddle1, mulpd (imag_b, imag_twiddle1));
					vec[2*(group + group_member + stride)+1]   = fmaddpd (imag_b, real_twiddle1, mulpd (real_b, imag_twiddle1));
					vec[2*(group + group_member + 2*stride)]   = fmsubpd (real_c, real_twiddle2, mulpd (imag_c, imag_twiddle2));
					vec[2*(group + group_member + 2*stride)+1] = fmaddpd (imag_c, real_twiddle2, mulpd (real_c, imag_twiddle2));
					vec[2*(group + group_member + 3*stride)]   = fmaddpd (real_d, real_twiddle2, mulpd (imag_d, imag_twiddle2));
					vec[2*(group + group_member + 3*stride)+1] = fmsubpd (imag_d, real_twiddle2, mulpd (real_d, imag_twiddle2));
					vec[2*(group + group_member + 4*stride)]   = fmaddpd (real_e, real_twiddle1, mulpd (imag_e, imag_twiddle1));
					vec[2*(group + group_member + 4*stride)+1] = fmsubpd (imag_e, real_twiddle1, mulpd (real_e, imag_twiddle1));
				}
			}
		}
	}

	// Radix-4 FFT, one fewer complex multiply than radix-2, does more work while data is in registers
	for ( ; size % 4 == 0; size /= 4, twiddle_stride *= 4) {
		stride = size / 4;
		for (unsigned int group = 0; group < fft_size; group += stride * 4) {
			for (unsigned int group_member = 0, twiddle_idx = 0; group_member < stride; group_member++, twiddle_idx += 2*twiddle_stride) {
				VDT real_a = vec[2*(group + group_member)];
				VDT imag_a = vec[2*(group + group_member)+1];
				VDT real_b = vec[2*(group + group_member + stride)];
				VDT imag_b = vec[2*(group + group_member + stride)+1];
				VDT real_c = vec[2*(group + group_member + 2*stride)];
				VDT imag_c = vec[2*(group + group_member + 2*stride)+1];
				VDT real_d = vec[2*(group + group_member + 3*stride)];
				VDT imag_d = vec[2*(group + group_member + 3*stride)+1];

				VDT real_new_a = addpd (real_a, real_c);			// a + c
				VDT imag_new_a = addpd (imag_a, imag_c);
				VDT real_new_c = subpd (real_a, real_c);			// a - c
				VDT imag_new_c = subpd (imag_a, imag_c);
				VDT real_new_b = addpd (real_b, real_d);			// b + d
				VDT imag_new_b = addpd (imag_b, imag_d);
				VDT real_new_d = subpd (real_b, real_d);			// b - d
				VDT imag_new_d = subpd (imag_b, imag_d);

				real_a = addpd (real_new_a, real_new_b);				// a + b
				imag_a = addpd (imag_new_a, imag_new_b);
				real_b = subpd (real_new_a, real_new_b);				// a - b
				imag_b = subpd (imag_new_a, imag_new_b);
				real_c = subpd (real_new_c, imag_new_d);				// c + id
				imag_c = addpd (imag_new_c, real_new_d);
				real_d = addpd (real_new_c, imag_new_d);				// c - id
				imag_d = subpd (imag_new_c, real_new_d);

				vec[2*(group + group_member)]   = real_a;
				vec[2*(group + group_member)+1] = imag_a;
				if (group_member == 0) {	// First group member's twiddles are 1+0i, skip twiddle multiply
					vec[2*(group + group_member + stride)]     = real_b;
					vec[2*(group + group_member + stride)+1]   = imag_b;
					vec[2*(group + group_member + 2*stride)]   = real_c;
					vec[2*(group + group_member + 2*stride)+1] = imag_c;
					vec[2*(group + group_member + 3*stride)]   = real_d;
					vec[2*(group + group_member + 3*stride)+1] = imag_d;
				}
				else {			// Apply twiddles
					VDT real_twiddle1 = broadcastsd (pmdata->twiddles2[2*twiddle_idx]);
					VDT imag_twiddle1 = broadcastsd (pmdata->twiddles2[2*twiddle_idx+1]);
					VDT real_twiddle2 = broadcastsd (pmdata->twiddles2[2*(twiddle_idx+1)]);
					VDT imag_twiddle2 = broadcastsd (pmdata->twiddles2[2*(twiddle_idx+1)+1]);
					vec[2*(group + group_member + stride)]     = fmsubpd (real_b, real_twiddle2, mulpd (imag_b, imag_twiddle2));
					vec[2*(group + group_member + stride)+1]   = fmaddpd (imag_b, real_twiddle2, mulpd (real_b, imag_twiddle2));
					vec[2*(group + group_member + 2*stride)]   = fmsubpd (real_c, real_twiddle1, mulpd (imag_c, imag_twiddle1));
					vec[2*(group + group_member + 2*stride)+1] = fmaddpd (imag_c, real_twiddle1, mulpd (real_c, imag_twiddle1));
					vec[2*(group + group_member + 3*stride)]   = fmaddpd (real_d, real_twiddle1, mulpd (imag_d, imag_twiddle1));
					vec[2*(group + group_member + 3*stride)+1] = fmsubpd (imag_d, real_twiddle1, mulpd (real_d, imag_twiddle1));
				}
			}
		}
	}

	// Radix-2 FFT
	if (size == 2) {
		stride = 1;		// Stride for this radix-2 step
		for (unsigned int group = 0; group < fft_size; group += stride * 2) {
			for (unsigned int group_member = 0; group_member < stride; group_member++) {
				VDT real_a = vec[2*(group + group_member)];
				VDT imag_a = vec[2*(group + group_member)+1];
				VDT real_b = vec[2*(group + group_member + stride)];
				VDT imag_b = vec[2*(group + group_member + stride)+1];
				vec[2*(group + group_member)] = addpd (real_a, real_b);
				vec[2*(group + group_member)+1] = addpd (imag_a, imag_b);
				vec[2*(group + group_member + stride)] = subpd (real_a, real_b);
				vec[2*(group + group_member + stride)+1] = subpd (imag_a, imag_b);
			}
		}
	}
}

void invfft_line (
	pmhandle *pmdata,		// Handle for polymult library
	VDT	*vec,			// Array of complex values
	unsigned int fft_size)		// The fft size
{
	unsigned int stride;			// Stride through FFT data
	unsigned int twiddle_stride;		// Stride through the twiddles data
	unsigned int powers_of_two;

	generate_sincos (pmdata, fft_size);

	// Radix-2 FFT
	stride = 1;										// First FFT data stride
	for (twiddle_stride = fft_size;	twiddle_stride % 3 == 0; twiddle_stride /= 3);		// First twiddle stride
	for (powers_of_two = 0; (fft_size & (1 << powers_of_two)) == 0; powers_of_two++);
	if (powers_of_two & 1) {		// If odd number of powers of two do a radix-2 step
		for (unsigned int group = 0; group < fft_size; group += stride * 2) {
			for (unsigned int group_member = 0; group_member < stride; group_member++) {
				VDT real_a = vec[2*(group + group_member)];
				VDT imag_a = vec[2*(group + group_member)+1];
				VDT real_b = vec[2*(group + group_member + stride)];
				VDT imag_b = vec[2*(group + group_member + stride)+1];
				vec[2*(group + group_member)] = addpd (real_a, real_b);
				vec[2*(group + group_member)+1] = addpd (imag_a, imag_b);
				vec[2*(group + group_member + stride)] = subpd (real_a, real_b);
				vec[2*(group + group_member + stride)+1] = subpd (imag_a, imag_b);
			}
		}
		stride *= 2;			// Next stride amount
		twiddle_stride /= 2;		// Next twiddle stride amount
		powers_of_two--;
	}

	// Radix-4 FFT, one fewer complex multiply than radix-2, does more work while data is in registers
	for ( ; powers_of_two; powers_of_two -= 2, stride *= 4) {
		twiddle_stride /= 4;
		for (unsigned int group = 0; group < fft_size; group += stride * 4) {
			for (unsigned int group_member = 0, twiddle_idx = 0; group_member < stride; group_member++, twiddle_idx += 2*twiddle_stride) {
				VDT real_a = vec[2*(group + group_member)];
				VDT imag_a = vec[2*(group + group_member)+1];
				VDT real_b = vec[2*(group + group_member + stride)];
				VDT imag_b = vec[2*(group + group_member + stride)+1];
				VDT real_c = vec[2*(group + group_member + 2*stride)];
				VDT imag_c = vec[2*(group + group_member + 2*stride)+1];
				VDT real_d = vec[2*(group + group_member + 3*stride)];
				VDT imag_d = vec[2*(group + group_member + 3*stride)+1];

				// Special case first group member, no twiddles necessary
				VDT real_twiddled_b, imag_twiddled_b, real_twiddled_c, imag_twiddled_c, real_twiddled_d, imag_twiddled_d;
				if (group_member == 0) {
					real_twiddled_b = real_b;
					imag_twiddled_b = imag_b;
					real_twiddled_c = real_c;
					imag_twiddled_c = imag_c;
					real_twiddled_d = real_d;
					imag_twiddled_d = imag_d;
				}
				// Apply the CONJUGATE of the twiddle
				else {
					VDT real_twiddle1 = broadcastsd (pmdata->twiddles2[2*twiddle_idx]);
					VDT imag_twiddle1 = broadcastsd (pmdata->twiddles2[2*twiddle_idx+1]);
					VDT real_twiddle2 = broadcastsd (pmdata->twiddles2[2*(twiddle_idx+1)]);
					VDT imag_twiddle2 = broadcastsd (pmdata->twiddles2[2*(twiddle_idx+1)+1]);
					real_twiddled_b = fmaddpd (real_b, real_twiddle2, mulpd (imag_b, imag_twiddle2));
					imag_twiddled_b = fmsubpd (imag_b, real_twiddle2, mulpd (real_b, imag_twiddle2));
					real_twiddled_c = fmaddpd (real_c, real_twiddle1, mulpd (imag_c, imag_twiddle1));
					imag_twiddled_c = fmsubpd (imag_c, real_twiddle1, mulpd (real_c, imag_twiddle1));
					real_twiddled_d = fmsubpd (real_d, real_twiddle1, mulpd (imag_d, imag_twiddle1));
					imag_twiddled_d = fmaddpd (imag_d, real_twiddle1, mulpd (real_d, imag_twiddle1));
				}

				VDT real_new_a = addpd (real_a, real_twiddled_b);			// a + b
				VDT imag_new_a = addpd (imag_a, imag_twiddled_b);
				VDT real_new_b = subpd (real_a, real_twiddled_b);			// a - b
				VDT imag_new_b = subpd (imag_a, imag_twiddled_b);
				VDT real_new_c = addpd (real_twiddled_c, real_twiddled_d);		// c + d
				VDT imag_new_c = addpd (imag_twiddled_c, imag_twiddled_d);
				VDT real_new_d = subpd (real_twiddled_c, real_twiddled_d);		// c - d
				VDT imag_new_d = subpd (imag_twiddled_c, imag_twiddled_d);

				real_a = addpd (real_new_a, real_new_c);					// a + c
				imag_a = addpd (imag_new_a, imag_new_c);
				real_c = subpd (real_new_a, real_new_c);					// a - c
				imag_c = subpd (imag_new_a, imag_new_c);
				real_b = addpd (real_new_b, imag_new_d);					// b - id
				imag_b = subpd (imag_new_b, real_new_d);
				real_d = subpd (real_new_b, imag_new_d);					// b + id
				imag_d = addpd (imag_new_b, real_new_d);

				// Save results
				vec[2*(group + group_member)]              = real_a;
				vec[2*(group + group_member)+1]            = imag_a;
				vec[2*(group + group_member + stride)]     = real_b;
				vec[2*(group + group_member + stride)+1]   = imag_b;
				vec[2*(group + group_member + 2*stride)]   = real_c;
				vec[2*(group + group_member + 2*stride)+1] = imag_c;
				vec[2*(group + group_member + 3*stride)]   = real_d;
				vec[2*(group + group_member + 3*stride)+1] = imag_d;
			}
		}
	}

	// Radix-5 FFT
	// R1= r1     +(r2+r5)     +(r3+r4)
	// R2= r1 +.309(r2+r5) -.809(r3+r4)    +.951(i2-i5) +.588(i3-i4)
	// R5= r1 +.309(r2+r5) -.809(r3+r4)    -.951(i2-i5) -.588(i3-i4)
	// R3= r1 -.809(r2+r5) +.309(r3+r4)    +.588(i2-i5) -.951(i3-i4)
	// R4= r1 -.809(r2+r5) +.309(r3+r4)    -.588(i2-i5) +.951(i3-i4)
	// I1= i1     +(i2+i5)     +(i3+i4)
	// I2= i1 +.309(i2+i5) -.809(i3+i4)    -.951(r2-r5) -.588(r3-r4)
	// I5= i1 +.309(i2+i5) -.809(i3+i4)    +.951(r2-r5) +.588(r3-r4)
	// I3= i1 -.809(i2+i5) +.309(i3+i4)    -.588(r2-r5) +.951(r3-r4)
	// I4= i1 -.809(i2+i5) +.309(i3+i4)    +.588(r2-r5) -.951(r3-r4)
	for ( ; (fft_size / stride) % 5 == 0; stride *= 5) {
		twiddle_stride /= 5;
		VDT _309 = broadcastsd (0.30901699437494742410229341718282);
		VDT _809 = broadcastsd (0.80901699437494742410229341718282);
		VDT _951 = broadcastsd (0.95105651629515357211643933337938);
		VDT _588_951 = broadcastsd (0.61803398874989484820458683436564);
		for (unsigned int group = 0; group < fft_size; group += stride * 5) {
			for (unsigned int group_member = 0, twiddle_idx = 0; group_member < stride; group_member++, twiddle_idx += 2*twiddle_stride) {
				VDT real_a = vec[2*(group + group_member)];
				VDT imag_a = vec[2*(group + group_member)+1];
				VDT real_b = vec[2*(group + group_member + stride)];
				VDT imag_b = vec[2*(group + group_member + stride)+1];
				VDT real_c = vec[2*(group + group_member + 2*stride)];
				VDT imag_c = vec[2*(group + group_member + 2*stride)+1];
				VDT real_d = vec[2*(group + group_member + 3*stride)];
				VDT imag_d = vec[2*(group + group_member + 3*stride)+1];
				VDT real_e = vec[2*(group + group_member + 4*stride)];
				VDT imag_e = vec[2*(group + group_member + 4*stride)+1];

				// Special case first group member, no twiddles necessary
				VDT real_twiddled_b, imag_twiddled_b, real_twiddled_c, imag_twiddled_c;
				VDT real_twiddled_d, imag_twiddled_d, real_twiddled_e, imag_twiddled_e;
				if (group_member == 0) {
					real_twiddled_b = real_b;
					imag_twiddled_b = imag_b;
					real_twiddled_c = real_c;
					imag_twiddled_c = imag_c;
					real_twiddled_d = real_d;
					imag_twiddled_d = imag_d;
					real_twiddled_e = real_e;
					imag_twiddled_e = imag_e;
				}
				// Apply the CONJUGATE of the twiddle
				else {
					VDT real_twiddle1 = broadcastsd (pmdata->twiddles2[2*twiddle_idx]);
					VDT imag_twiddle1 = broadcastsd (pmdata->twiddles2[2*twiddle_idx+1]);
					VDT real_twiddle2 = broadcastsd (pmdata->twiddles2[2*(twiddle_idx+1)]);
					VDT imag_twiddle2 = broadcastsd (pmdata->twiddles2[2*(twiddle_idx+1)+1]);
					real_twiddled_b = fmaddpd (real_b, real_twiddle1, mulpd (imag_b, imag_twiddle1));
					imag_twiddled_b = fmsubpd (imag_b, real_twiddle1, mulpd (real_b, imag_twiddle1));
					real_twiddled_c = fmaddpd (real_c, real_twiddle2, mulpd (imag_c, imag_twiddle2));
					imag_twiddled_c = fmsubpd (imag_c, real_twiddle2, mulpd (real_c, imag_twiddle2));
					real_twiddled_d = fmsubpd (real_d, real_twiddle2, mulpd (imag_d, imag_twiddle2));
					imag_twiddled_d = fmaddpd (imag_d, real_twiddle2, mulpd (real_d, imag_twiddle2));
					real_twiddled_e = fmsubpd (real_e, real_twiddle1, mulpd (imag_e, imag_twiddle1));
					imag_twiddled_e = fmaddpd (imag_e, real_twiddle1, mulpd (real_e, imag_twiddle1));
				}

				VDT real_b_plus_e = addpd (real_twiddled_b, real_twiddled_e);				// b + e
				VDT imag_b_plus_e = addpd (imag_twiddled_b, imag_twiddled_e);
				VDT real_c_plus_d = addpd (real_twiddled_c, real_twiddled_d);				// c + d
				VDT imag_c_plus_d = addpd (imag_twiddled_c, imag_twiddled_d);
				VDT real_b_minus_e = subpd (real_twiddled_b, real_twiddled_e);			// b - e
				VDT imag_b_minus_e = subpd (imag_twiddled_b, imag_twiddled_e);
				VDT real_c_minus_d = subpd (real_twiddled_c, real_twiddled_d);			// c - d
				VDT imag_c_minus_d = subpd (imag_twiddled_c, imag_twiddled_d);

				VDT real_tmp25a = fmaddpd (_309, real_b_plus_e, fnmaddpd (_809, real_c_plus_d, real_a)); // a +.309*(b+e) -.809*(c+d)
				VDT imag_tmp25a = fmaddpd (_309, imag_b_plus_e, fnmaddpd (_809, imag_c_plus_d, imag_a));
				VDT real_tmp34a = fnmaddpd (_809, real_b_plus_e, fmaddpd (_309, real_c_plus_d, real_a)); // a -.809*(b+e) +.309*(c+d)
				VDT imag_tmp34a = fnmaddpd (_809, imag_b_plus_e, fmaddpd (_309, imag_c_plus_d, imag_a));
//GW: delay mul by .951 on true FMA CPUs
				VDT real_tmp25b = mulpd (_951, fmaddpd (_588_951, imag_c_minus_d, imag_b_minus_e));	// .951*(b-e) + 0.588*(c-d)
				VDT imag_tmp25b = mulpd (_951, fmaddpd (_588_951, real_c_minus_d, real_b_minus_e));
				VDT real_tmp34b = mulpd (_951, fmsubpd (_588_951, imag_b_minus_e, imag_c_minus_d));	// .588*(b-e) - 0.951*(c-d)
				VDT imag_tmp34b = mulpd (_951, fmsubpd (_588_951, real_b_minus_e, real_c_minus_d));

				real_a = addpd (addpd (real_a, real_b_plus_e), real_c_plus_d);	// a + b + c + d + e
				imag_a = addpd (addpd (imag_a, imag_b_plus_e), imag_c_plus_d);
				real_b = addpd (real_tmp25a, real_tmp25b);				// combine 25a and 25b
				imag_b = subpd (imag_tmp25a, imag_tmp25b);
				real_e = subpd (real_tmp25a, real_tmp25b);
				imag_e = addpd (imag_tmp25a, imag_tmp25b);
				real_c = addpd (real_tmp34a, real_tmp34b);				// combine 34a and 34b
				imag_c = subpd (imag_tmp34a, imag_tmp34b);
				real_d = subpd (real_tmp34a, real_tmp34b);
				imag_d = addpd (imag_tmp34a, imag_tmp34b);

				// Save results
				vec[2*(group + group_member)]              = real_a;
				vec[2*(group + group_member)+1]            = imag_a;
				vec[2*(group + group_member + stride)]     = real_b;
				vec[2*(group + group_member + stride)+1]   = imag_b;
				vec[2*(group + group_member + 2*stride)]   = real_c;
				vec[2*(group + group_member + 2*stride)+1] = imag_c;
				vec[2*(group + group_member + 3*stride)]   = real_d;
				vec[2*(group + group_member + 3*stride)+1] = imag_d;
				vec[2*(group + group_member + 4*stride)]   = real_e;
				vec[2*(group + group_member + 4*stride)+1] = imag_e;
			}
		}
	}

	// Radix-3 FFT
	// A 3-complex inverse FFT is:
	// Res1: (R1+R2+R3) + (I1+I2+I3)i
	// Res2:  (R1-.5R2-.5R3+.866I2-.866I3) + (I1-.5I2-.5I3-.866R2+.866R3)i
	// Res3:  (R1-.5R2-.5R3-.866I2+.866I3) + (I1-.5I2-.5I3+.866R2-.866R3)i
	twiddle_stride = fft_size / stride;		// First radix-3 twiddle stride
	for ( ; stride < fft_size; stride *= 3) {
		twiddle_stride /= 3;
		VDT half = broadcastsd (0.5);
		VDT _866 = broadcastsd (0.86602540378443864676372317075294);
		for (unsigned int group = 0; group < fft_size; group += stride * 3) {
			for (unsigned int group_member = 0, twiddle_idx = 0; group_member < stride; group_member++, twiddle_idx += twiddle_stride) {
				VDT real_a = vec[2*(group + group_member)];
				VDT imag_a = vec[2*(group + group_member)+1];
				VDT real_b = vec[2*(group + group_member + stride)];
				VDT imag_b = vec[2*(group + group_member + stride)+1];
				VDT real_c = vec[2*(group + group_member + 2*stride)];
				VDT imag_c = vec[2*(group + group_member + 2*stride)+1];

				// Special case first group member, no twiddles necessary
				VDT real_twiddled_b, imag_twiddled_b, real_twiddled_c, imag_twiddled_c;
				if (group_member == 0) {
					real_twiddled_b = real_b;
					imag_twiddled_b = imag_b;
					real_twiddled_c = real_c;
					imag_twiddled_c = imag_c;
				}
				// Apply the CONJUGATE of the twiddle
				else {
					VDT real_twiddle1 = broadcastsd (pmdata->twiddles1[2*twiddle_idx]);
					VDT imag_twiddle1 = broadcastsd (pmdata->twiddles1[2*twiddle_idx+1]);
					real_twiddled_b = addpd (mulpd (real_b, real_twiddle1), mulpd (imag_b, imag_twiddle1));
					imag_twiddled_b = subpd (mulpd (imag_b, real_twiddle1), mulpd (real_b, imag_twiddle1));
					real_twiddled_c = subpd (mulpd (real_c, real_twiddle1), mulpd (imag_c, imag_twiddle1));
					imag_twiddled_c = addpd (mulpd (imag_c, real_twiddle1), mulpd (real_c, imag_twiddle1));
				}

				VDT real_b_plus_c = addpd (real_twiddled_b, real_twiddled_c);				// b + c
				VDT imag_b_plus_c = addpd (imag_twiddled_b, imag_twiddled_c);

				VDT real_tmp23a = fnmaddpd (half, real_b_plus_c, real_a);				// a - 0.5 * (b + c)
				VDT imag_tmp23a = fnmaddpd (half, imag_b_plus_c, imag_a);
//GW: delay mul by 866 on true FMA CPUs
				VDT real_tmp23b = mulpd (_866, subpd (imag_twiddled_b, imag_twiddled_c));	// 0.866 * (b - c)
				VDT imag_tmp23b = mulpd (_866, subpd (real_twiddled_b, real_twiddled_c));

				real_a = addpd (real_a, real_b_plus_c);							// a + b + c
				imag_a = addpd (imag_a, imag_b_plus_c);
				real_b = addpd (real_tmp23a, real_tmp23b);						// combine 23a and 23b
				imag_b = subpd (imag_tmp23a, imag_tmp23b);
				real_c = subpd (real_tmp23a, real_tmp23b);
				imag_c = addpd (imag_tmp23a, imag_tmp23b);

				vec[2*(group + group_member)] = real_a;
				vec[2*(group + group_member)+1] = imag_a;
				vec[2*(group + group_member + stride)] = real_b;
				vec[2*(group + group_member + stride)+1] = imag_b;
				vec[2*(group + group_member + 2*stride)] = real_c;
				vec[2*(group + group_member + 2*stride)+1] = imag_c;
			}
		}
	}
}

// out = in1 * in2 pointwise also divide by FFT length
void pointmul_line (VDT *invec1, VDT *invec2, VDT *outvec, int size)
{
	VDT inv_fft_size = broadcastsd (1.0 / (double) size);
	for (int j = 0; j < size; j++) {
		VDT real_x = invec1[2*j];
		VDT imag_x = invec1[2*j+1];
		VDT real_y = invec2[2*j];
		VDT imag_y = invec2[2*j+1];
		outvec[2*j] = mulpd (fmsubpd (real_x, real_y, mulpd (imag_x, imag_y)), inv_fft_size);
		outvec[2*j+1] = mulpd (fmaddpd (real_x, imag_y, mulpd (imag_x, real_y)), inv_fft_size);
	}
}

// outvec += input vector shifted the proper amount to create a monic polymult result
void monic_line_add_hi (VDT *invec, int invec_size, VDT *outvec, int outvec_size)
{
	// Add words starting from the highest word
	for (invec += 2*(invec_size-1), outvec += 2*(outvec_size-1); invec_size; invec -= 2, outvec -= 2, invec_size--) {
		outvec[0] = addpd (outvec[0], invec[0]);
		outvec[1] = addpd (outvec[1], invec[1]);
	}
}

// outvec += input vector shifted to the low part of outvec
void monic_line_add_lo (VDT *invec, int invec_size, VDT *outvec, int outvec_size)
{
	for ( ; invec_size; invec += 2, outvec += 2, invec_size--) {
		outvec[0] = addpd (outvec[0], invec[0]);
		outvec[1] = addpd (outvec[1], invec[1]);
	}
}

// Make adjustments to outvec for any monic input vectors
void monic_line_adjustments (VDT *invec1, int invec1_size, VDT *invec2, int invec2_size, VDT *outvec, int outvec_size, int options)
{
	// If neither input is monic, no adjustments are necessary
	if (! (options & POLYMULT_INVEC1_MONIC) && ! (options & POLYMULT_INVEC2_MONIC)) return;

	// If either input is a monic RLP, create a new entry before the output vector calculated thusfar
	if (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) ||
	    ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP))) {
		outvec -= 2;
		outvec_size++;
		memset (&outvec[0], 0, 2*sizeof(VDT));
	}

	// Create a new entry after the output vector calculated thusfar so that we can add to it
	outvec_size++;
	memset (&outvec[2*(outvec_size-1)], 0, 2*sizeof(VDT));

	// Adjust for an implied high coefficient of one -- add in the other poly in the high words
	if (options & POLYMULT_INVEC1_MONIC) monic_line_add_hi (invec2, invec2_size, outvec, outvec_size);
	if (options & POLYMULT_INVEC2_MONIC) monic_line_add_hi (invec1, invec1_size, outvec, outvec_size);

	// Adjust for an implied low coefficient of one -- add in the other poly to the low words
	if ((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) monic_line_add_lo (invec2, invec2_size, outvec, outvec_size);
	if ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP)) monic_line_add_lo (invec1, invec1_size, outvec, outvec_size);
}

// Select from various polymult algorithms
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
	int	adjusted_invec1_size, adjusted_invec2_size, adjusted_outvec_size, fft_size;
	bool	must_brute, must_fft, post_process_monics;
	VDT	*invec1, *invec2, *outvec, *tmpvec;

	// Adjust sizes due to RLPs.  Calculate the output size assuming no monic inputs.
	adjusted_invec1_size = (options & POLYMULT_INVEC1_RLP) ? 2 * invec1_size - 1 : invec1_size;
	adjusted_invec2_size = (options & POLYMULT_INVEC2_RLP) ? 2 * invec2_size - 1 : invec2_size;
	adjusted_outvec_size = adjusted_invec1_size + adjusted_invec2_size - 1;

	// Sometimes we must (or it is beneficial) handle monic input polynomials as part of invec1/invec2 as opposed to a post-processing call to
	// monic_line_adjustments.  Find those cases and adjust vector sizes.
	if (options & POLYMULT_CIRCULAR) {
		ASSERTG (!(options & POLYMULT_INVEC1_MONIC) || !(options & POLYMULT_INVEC2_MONIC));
//GW: Find more cases where monic processing would be better done this way
		if (options & POLYMULT_INVEC1_MONIC) adjusted_invec1_size += (options & POLYMULT_INVEC1_RLP) ? 2 : 1;
		if (options & POLYMULT_INVEC2_MONIC) adjusted_invec2_size += (options & POLYMULT_INVEC2_RLP) ? 2 : 1;
		adjusted_outvec_size = outvec_size;
		post_process_monics = FALSE;
	} else
		post_process_monics = TRUE;

	// Polymult using AVX instructions.  Allocate memory.  For the output vector allocate two extra in case there are any monic RLP inputs.
	must_brute = (options & POLYMULT_CIRCULAR) && outvec_size != polymult_fft_size (outvec_size);
	must_fft = !must_brute && (options & (POLYMULT_FORWARD_FFT | POLYMULT_INVEC1_FFTED | POLYMULT_INVEC2_FFTED));
	if (must_brute || (!must_fft && adjusted_outvec_size < pmdata->KARAT_BREAK)) {
		invec1 = (VDT *) aligned_malloc (2 * adjusted_invec1_size * sizeof (VDT), 64);
		invec2 = (VDT *) aligned_malloc (2 * adjusted_invec2_size * sizeof (VDT), 64);
		outvec = (VDT *) aligned_malloc (2 * (adjusted_outvec_size + 2) * sizeof (VDT), 64) + 2;
		tmpvec = NULL;
	}
	else if (!must_fft && !(options & POLYMULT_CIRCULAR) && adjusted_outvec_size < pmdata->FFT_BREAK) {
		invec1 = (VDT *) aligned_malloc (2 * adjusted_invec1_size * sizeof (VDT), 64);
		invec2 = (VDT *) aligned_malloc (2 * adjusted_invec2_size * sizeof (VDT), 64);
		outvec = (VDT *) aligned_malloc (2 * (adjusted_outvec_size + 2) * sizeof (VDT), 64) + 2;
		tmpvec = (VDT *) aligned_malloc (2 * (adjusted_outvec_size + 32) * sizeof (VDT), 64);	// Could need one extra for each recursion
	}
	else {
		fft_size = (options & (POLYMULT_FORWARD_FFT | POLYMULT_CIRCULAR)) ? outvec_size :
			   (options & POLYMULT_INVEC1_FFTED) ? invec1_size :
			   (options & POLYMULT_INVEC2_FFTED) ? invec2_size : polymult_fft_size (adjusted_outvec_size);
		tmpvec = (VDT *) aligned_malloc (2 * fft_size * sizeof (VDT), 64);		// FFT of invec1 goes here
		outvec = (VDT *) aligned_malloc (2 * (fft_size + 2) * sizeof (VDT), 64) + 2;	// FFT of invec2 goes here as does outvec
		if (post_process_monics && (options & POLYMULT_INVEC2_MONIC)) invec1 = (VDT *) aligned_malloc (2 * fft_size * sizeof (VDT), 64);
		else invec1 = tmpvec;
		if (post_process_monics && (options & POLYMULT_INVEC1_MONIC)) invec2 = (VDT *) aligned_malloc (2 * fft_size * sizeof (VDT), 64);
		else invec2 = outvec;
	}

	// Grab four complex values from each gwnum coefficient.  Work that data set.  Then move onto remaining complex vals in gwnums.
	for ( ; ; ) {
		int line = next_line (pmdata);
		if (!read_line (pmdata->gwdata, gw_invec1, invec1_size, line, invec1, options & INVEC1_OPTIONS, post_process_monics)) break;
		read_line (pmdata->gwdata, gw_invec2, invec2_size, line, invec2, options & INVEC2_OPTIONS, post_process_monics);

		if (must_brute || (!must_fft && adjusted_outvec_size < pmdata->KARAT_BREAK)) {
			brute_line (invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, outvec, adjusted_outvec_size, options, gw_outvec);
		}
		else if (!must_fft && !(options & POLYMULT_CIRCULAR) && adjusted_outvec_size < pmdata->FFT_BREAK) {
			karatsuba_line (pmdata, invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, tmpvec, outvec);
		}
		else {
			if (!(options & POLYMULT_INVEC1_FFTED)) {
				if (invec1 != tmpvec) memcpy (tmpvec, invec1, 2*adjusted_invec1_size*sizeof(VDT));
				memset (tmpvec + 2*adjusted_invec1_size, 0, 2*(fft_size-adjusted_invec1_size)*sizeof(VDT));
				fft_line (pmdata, tmpvec, fft_size);
			}
			if (options & POLYMULT_FORWARD_FFT) {
				write_line (pmdata->gwdata, gw_outvec, fft_size, line, tmpvec, 0);
				continue;
			}
			if (!(options & POLYMULT_INVEC2_FFTED)) {
				if (invec2 != outvec) memcpy (outvec, invec2, 2*adjusted_invec2_size*sizeof(VDT));
				memset (outvec + 2*adjusted_invec2_size, 0, 2*(fft_size-adjusted_invec2_size)*sizeof(VDT));
				fft_line (pmdata, outvec, fft_size);
			}
			pointmul_line (tmpvec, outvec, outvec, fft_size);
			invfft_line (pmdata, outvec, fft_size);
		}
		if (post_process_monics) monic_line_adjustments (invec1, adjusted_invec1_size, invec2, adjusted_invec2_size, outvec, adjusted_outvec_size, options);
		write_line (pmdata->gwdata, gw_outvec, outvec_size, line, outvec, options & POLYMULT_CIRCULAR ? 0 : options);
	}

	// Free memory
	if (invec1 != tmpvec) aligned_free (invec1);
	if (invec2 != outvec) aligned_free (invec2);
	aligned_free (outvec - 2);
	aligned_free (tmpvec);
}

#endif


/*--------------------------------------------------------------------------
|			Polymult helper thread
+-------------------------------------------------------------------------*/

#if !defined (SSE2) && !defined (AVX) && !defined (FMA) && !defined (AVX512)

/* Entry point for helper thread to do part of the polymult */
void polymult_thread (void *arg)
{
	pmhandle *pmdata = (pmhandle *) arg;

/* Generate a unique helper thread number (so that callback routine to uniquely identify helper threads) */
/* Call optional user provided callback routine so that the caller can set the thread's priority and affinity */

	gwmutex_lock (&pmdata->poly_mutex);
	int thread_num = ++pmdata->next_thread_num;
	gwmutex_unlock (&pmdata->poly_mutex);
	if (pmdata->gwdata->thread_callback != NULL)
		(*pmdata->gwdata->thread_callback) (thread_num, 20, pmdata->gwdata->thread_callback_data);

/* Loop doing polymult work */

	while (!pmdata->termination_in_progress) {
		// Optimized routines for AVX512 instructions
		if (pmdata->gwdata->cpu_flags & CPU_AVX512F)
			polymult_line_avx512 (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					      pmdata->outvec, pmdata->outvec_size, pmdata->options);
		// Optimized routines for FMA3 instructions
		else if (pmdata->gwdata->cpu_flags & CPU_FMA3)
			polymult_line_fma (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					   pmdata->outvec, pmdata->outvec_size, pmdata->options);
		// Optimized routines for AVX instructions
		else if (pmdata->gwdata->cpu_flags & CPU_AVX)
			polymult_line_avx (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					   pmdata->outvec, pmdata->outvec_size, pmdata->options);
		// Optimized routines for SSE2 instructions
		else if (pmdata->gwdata->cpu_flags & CPU_SSE2)
			polymult_line_sse2 (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					    pmdata->outvec, pmdata->outvec_size, pmdata->options);
		// Default routines using no vector instructions
#ifndef X86_64
		else
			polymult_line_dbl (pmdata, pmdata->invec1, pmdata->invec1_size, pmdata->invec2, pmdata->invec2_size,
					   pmdata->outvec, pmdata->outvec_size, pmdata->options);
#endif

		// When all helpers are finished, signal main thread
		gwmutex_lock (&pmdata->poly_mutex);
		int num_active = --pmdata->helpers_active;
		gwmutex_unlock (&pmdata->poly_mutex);
		if (num_active == 0) gwevent_signal (&pmdata->helpers_done);

/* Wait for more work (or termination) */

		gwevent_wait (&pmdata->work_to_do, 0);
	}

/* Call optional user provided callback routine so that the caller can do any necessary cleanup */

	if (pmdata->gwdata->thread_callback != NULL)
		(*pmdata->gwdata->thread_callback) (thread_num, 21, pmdata->gwdata->thread_callback_data);
}

/*--------------------------------------------------------------------------
|			Polymult main routine
+-------------------------------------------------------------------------*/

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
	int	outvec_size,		// Size of the output polynomial (or fft_size if POLYMULT_FORWARD_FFT or POLYMULT_CIRCULAR)
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
	ASSERTG ((options & (POLYMULT_FORWARD_FFT | POLYMULT_CIRCULAR | POLYMULT_INVEC1_FFTED | POLYMULT_INVEC2_FFTED) || outvec_size == outsz2));
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
	pmdata->invec1 = invec1;
	pmdata->invec1_size = invec1_size;
	pmdata->invec2 = invec2;
	pmdata->invec2_size = invec2_size;
	pmdata->outvec = outvec;
	pmdata->outvec_size = outvec_size;
	pmdata->options = options;
	pmdata->next_line = 0;

	// Fire up each helper thread
	if (pmdata->num_threads > 1) {
		pmdata->helpers_active = pmdata->num_threads - 1;
		if (pmdata->thread_ids == NULL) {
			pmdata->thread_ids = (gwthread *) malloc (pmdata->num_threads * sizeof (gwthread));
			gwmutex_init (&pmdata->poly_mutex);
			gwevent_init (&pmdata->work_to_do);
			gwevent_init (&pmdata->helpers_done);
			pmdata->next_thread_num = 0;
			for (int i = 1; i < pmdata->num_threads; i++) gwthread_create_waitable (&pmdata->thread_ids[i], &polymult_thread, (void *) pmdata);
		} else {
			gwmutex_lock (&pmdata->poly_mutex);
			gwevent_reset (&pmdata->helpers_done);
			gwevent_signal (&pmdata->work_to_do);		// Start all helper threads
			gwevent_reset (&pmdata->work_to_do);		// Reset event immediately
			gwmutex_unlock (&pmdata->poly_mutex);
		}
	}

	// This thread works on the polymult too.  Select which optimized routines to use.
	if (pmdata->gwdata->cpu_flags & CPU_AVX512F)
		polymult_line_avx512 (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
	else if (pmdata->gwdata->cpu_flags & CPU_FMA3)
		polymult_line_fma (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
	else if (pmdata->gwdata->cpu_flags & CPU_AVX)
		polymult_line_avx (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
	else if (pmdata->gwdata->cpu_flags & CPU_SSE2)
		polymult_line_sse2 (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
#ifndef X86_64
	else
		polymult_line_dbl (pmdata, invec1, invec1_size, invec2, invec2_size, outvec, outvec_size, options);
#endif

	// Wait for all helper threads to complete
	if (pmdata->num_threads > 1) gwevent_wait (&pmdata->helpers_done, 0);

	// Figure out which output entries (if any) need to get 1.0 added to the result because a monic poly is multiplied by a monic RLP poly.
	// We used to add the FFT(1) in monic_adjustments but this led to excessive round off error as FFT(1) is NOT a random signal.  It only
	// affects the least significant bits and is apt to get drowned out by the many multiplied signals.

	int	addin1 = -1;
	int	addin2 = -1;
	if ((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP)) addin1 = outvec_size - (2 * invec2_size);
	if ((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) addin2 = outvec_size - (2 * invec1_size);
	if (options & POLYMULT_INVEC2_MONIC_TIMES_MONIC_RLP_OK) addin2 = invec2_size;
			
	// Unfft all output coefficients
	if (!(options & POLYMULT_FORWARD_FFT)) {
		for (int j = 0; j < outvec_size; j++) {
			if (outvec[j] == NULL) continue;
			FFT_state (outvec[j]) = FULLY_FFTed;
			gwunfft2 (pmdata->gwdata, outvec[j], outvec[j], (j != addin1 && j != addin2 && (options & POLYMULT_STARTNEXTFFT)) ? GWMUL_STARTNEXTFFT : 0);
			if (j == addin1) gwaddsmall (pmdata->gwdata, outvec[j], 1);
			if (j == addin2) gwaddsmall (pmdata->gwdata, outvec[j], 1);
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

