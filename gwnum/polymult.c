/*----------------------------------------------------------------------
| polymult.c
|
| This file contains the C routines for fast polynomial multiplication
| 
|  Copyright 2021-2022 Mersenne Research, Inc.  All rights reserved.
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
// 1) Support non-power-of-two FFTs below 32?  Support FFT size requiring a radix-2 step?
// 9) Is a fast convert-gwnum-to-larger-FFT feasible? necessary?
// 10) Tune for Brute / Karatsuba / FFT crossovers.
// 12) Tune routine that returns needed safety margin needed for given poly sizes
// 15) Add polymult to gwtest or as a standalone program
// 16) Handle out-of-memory gracefully
// 17) Right collection of PFA FFT sizes?
// 19) For large ffts, use Pavel's idea of copying the complex values to a gwnum and use an unweighted FFT with a special no-carry-propagation norm routine.
//     Perhaps only do this for most common architectures, such as AVX-512 and FMA 64-bit.
// 22) Can we add asserts that polymult isn't violating EXTRA_BITS?
// 25) Optimize.  Read/write several cache lines at a time (fill up L2 cache?).  Prefetch future cache lines.
// 28) Optimize brute_force to do less compute when one or both invec entries are RLPs.
// 32) If monic input and polymult is doing an FFT and there is zero padding going on, then pad with one (or GW_FFT1) rather than use monic_add.  Be careful
//	that this does not cause any unwanted wrapping of data!
// 33) Are there cases with brute/Karatsuba would benefit from using a bigger poly rather than using monic_add?
// 34) Support GWMUL_BY_CONST, GWMUL_ADD_CONST option in gwunfft2.  No, use the no_UNFFT2 option (requires user to know where 1*1 took place)??
// 35) A monic 1024 times monic 1024 using length 2048 FFT could apply the two ones in the FFT vectors and patch outvec[0] (possibly with gwsetaddin(-1))??how?
//	A monic 1023 times monic 1023 using length 2048 FFT could apply the two ones in the FFT vectors with no patching rather than monic  post-processing.
//	That should be a little faster.  Other monic combinations have similar optimizations.
// 38)  For poly #2, is there a quick block compression we can do where the extra read/write cost will be more than compensated by the few more extra
//	gwnums we can add to poly #2?
// 39)  Think more on unbalanced poly mults.  Is a 100 x 200 mult better done as two 100 x 100 multiplies?  Answer: Never for brute force, maybe for FFT.
//	If we can save the smaller RLP in FFTed state using smaller FFT size that might influence our decision.
// 41) Add callback routine to check for ^C/ESC
// 42) Transparent / or easy way to store gwnums on disk as part of input and/or output vectors.  Picture using 12GB memory and 50GB disk for really huge polys.
// 43) Read_line/write_line that reads/writes data with stride=1?
// 44) Cache several sincos tables for different FFT sizes.
// 45) If twiddle_stride gets large, build a smaller twiddle table for smaller strides.
// 46) Look for more FMA opportunities.
// 48) Test if hyperthreading is beneficial
// 51) Scratch area (like prime95) to avoid large strides,  Read_partial_line could populate the scratch area (makes monic_adjustment difficult)?
// 52) Half-sized DCT for RLPs?  The DJB FFTs may make reordering DCT outputs for an inverse FFT tricky.
// 53) Can P-1's circular monic-poly times monic-RLP be safely done without monic_adjustments allowing us to save the FFT'ed RLP that is half-sized due to DCT?
// 54) Does monolithic two-forward-FFTs-pointmul-and-one-inverse-FFT routine suggest different preferred pass sizes (use half the L2 cache)?
// 56) Option to not mul by 1/FFTLEN?
// 57) Option to save memory by processing AVX-512 data using FMA3.  (reduces sizeof(CVDT)).
// 58) Should we multi-thread FFTing the polymult input arguments?  Should we multithread polymult_preprocess?
// 60) In polymult_line_preprocess, most cases do not need an invec1 allocation.  We can write directly into the output buffer.  I think only RLPs without
//	POLYMULT_FFT are an issue.  Perhaps fix readline not to duplicate RLP data for when preprocessing without POLYMULT_FFT, or allocate a little more
//	space in the output buffer.
// 61) Can loop unrolling help in radix-3/4/5 sections?  Seems unlikely.
// 62) Linux stage 2 timings vary widely from 6.7 to 7.0 sec.  What explains the 5% slower runs?  4KB/64KB alignments?
// 63) Multithread large poly FFTs to reduce memory usage
// 64) Save memory in polymult_line.  For large FFTs with big buffer requirements, offer to re-read invec1 and/or invec2 for post-monic addin (much the same
//	way FMAvec is read in without requiring a buffer).  This is awkward if the vector being re-read has been pre-processed.
// 65) An fmavec of NULL could mean all zeros.  Combined with FNMADD would turn into negate outputs.  Code handles this, ASSERTS do not.
// 66) Offer polymult option where caller guarantees that invec1 and/or invec2 have all been FFTed??  Too minor an optimization?
// 67) Offer polymult option where plan is saved??  Too minor an optimization?
// 69) When compressing with POLYMULT_PRE_COMPRESS, figure out a way to malloc a smaller buffer and grow if necessary rather than malloc a large one and
//	shrink at end.

// Forward declarations required for the SSE2, AVX, FMA, AVX512 compiles
unsigned long cache_line_offset (gwhandle *gwdata, unsigned long i);
void generate_sincos (pmhandle *pmdata, const unsigned int size);
void zpad7_fft (gwhandle *gwdata, gwnum g);
int compress_vec (char *buf, int size);
void decompress_vec (char *inbuf, char *outbuf, int size);
gwarray polymult_line_preprocess_sse2 (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, int invec2_size, int outvec_size, int options);
gwarray polymult_line_preprocess_avx (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, int invec2_size, int outvec_size, int options);
gwarray polymult_line_preprocess_fma (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, int invec2_size, int outvec_size, int options);
gwarray polymult_line_preprocess_avx512 (pmhandle *pmdata, gwnum *gw_invec1, int invec1_size, int invec2_size, int outvec_size, int options);
void polymult_line_sse2 (pmhandle *pmdata);
void polymult_line_avx (pmhandle *pmdata);
void polymult_line_fma (pmhandle *pmdata);
void polymult_line_avx512 (pmhandle *pmdata);

#define POLYMULT_IMPL_BRUTE		0		// Brute force polymult
#define POLYMULT_IMPL_KARATSUBA		1		// Karasuba polymult
#define POLYMULT_IMPL_FFT		2		// FFT based polymult

#define INVEC1_OPTIONS	(POLYMULT_INVEC1_MONIC | POLYMULT_INVEC1_RLP | POLYMULT_INVEC1_NEGATE)
#define INVEC2_OPTIONS	(POLYMULT_INVEC2_MONIC | POLYMULT_INVEC2_RLP | POLYMULT_INVEC2_NEGATE)

#if defined (AVX512)
#define VLEN			512			// Vector is 256 bits wide
#define VDT			__m512d			// Vector data type
#define VFMA			1			// FMA is supported
#define	read_line		read_line_avx512	// Naming extension for this option
#define	read_preprocess_line	read_preprocess_line_avx512
#define	write_line		write_line_avx512
#define	brute_line		brute_line_avx512
#define	karatsuba_line		karatsuba_line_avx512
#define	fft_line		fft_line_avx512
#define pick_pass_sizes		pick_pass_sizes_avx512
#define	monic_line_add_hi	monic_line_add_hi_avx512
#define	monic_line_add_lo	monic_line_add_lo_avx512
#define	monic_line_adjustments	monic_line_adjustments_avx512
#define	polymult_line_preprocess polymult_line_preprocess_avx512
#define	polymult_line		polymult_line_avx512
#define broadcastsd		_mm512_set1_pd
#define addpd			_mm512_add_pd
#define subpd			_mm512_sub_pd
#define mulpd			_mm512_mul_pd
#define divpd			_mm512_div_pd
#define fmaddpd			_mm512_fmadd_pd
#define fmsubpd			_mm512_fmsub_pd
#define fnmaddpd		_mm512_fnmadd_pd
#elif defined (FMA)
#define VLEN			256			// Vector is 256 bits wide
#define VDT			__m256d			// Vector data type
#define VFMA			1			// FMA is supported
#define	read_line		read_line_fma		// Naming extension for this option
#define	read_preprocess_line	read_preprocess_line_fma
#define	write_line		write_line_fma
#define	brute_line		brute_line_fma
#define	karatsuba_line		karatsuba_line_fma
#define	fft_line		fft_line_fma
#define pick_pass_sizes		pick_pass_sizes_fma
#define	monic_line_add_hi	monic_line_add_hi_fma
#define	monic_line_add_lo	monic_line_add_lo_fma
#define	monic_line_adjustments	monic_line_adjustments_fma
#define	polymult_line_preprocess polymult_line_preprocess_fma
#define	polymult_line		polymult_line_fma
#define broadcastsd		_mm256_set1_pd
#define addpd			_mm256_add_pd
#define subpd			_mm256_sub_pd
#define mulpd			_mm256_mul_pd
#define divpd			_mm256_div_pd
#define fmaddpd			_mm256_fmadd_pd
#define fmsubpd			_mm256_fmsub_pd
#define fnmaddpd		_mm256_fnmadd_pd
#elif defined (AVX)
#define VLEN			256			// Vector is 256 bits wide
#define VDT			__m256d			// Vector data type
#define VFMA			0			// FMA not supported
#define	read_line		read_line_avx		// Naming extension for this option
#define	read_preprocess_line	read_preprocess_line_avx
#define	write_line		write_line_avx
#define	brute_line		brute_line_avx
#define	karatsuba_line		karatsuba_line_avx
#define	fft_line		fft_line_avx
#define pick_pass_sizes		pick_pass_sizes_avx
#define	monic_line_add_hi	monic_line_add_hi_avx
#define	monic_line_add_lo	monic_line_add_lo_avx
#define	monic_line_adjustments	monic_line_adjustments_avx
#define	polymult_line_preprocess polymult_line_preprocess_avx
#define	polymult_line		polymult_line_avx
#define broadcastsd		_mm256_set1_pd
#define addpd			_mm256_add_pd
#define subpd			_mm256_sub_pd
#define mulpd			_mm256_mul_pd
#define divpd			_mm256_div_pd
#elif defined (SSE2)
#define VLEN			128			// Vector is 256 bits wide
#define VDT			__m128d			// Vector data type
#define VFMA			0			// FMA not supported
#define	read_line		read_line_sse2		// Naming extension for this option
#define	read_preprocess_line	read_preprocess_line_sse2
#define	write_line		write_line_sse2
#define	brute_line		brute_line_sse2
#define	karatsuba_line		karatsuba_line_sse2
#define	fft_line		fft_line_sse2
#define pick_pass_sizes		pick_pass_sizes_sse2
#define	monic_line_add_hi	monic_line_add_hi_sse2
#define	monic_line_add_lo	monic_line_add_lo_sse2
#define	monic_line_adjustments	monic_line_adjustments_sse2
#define	polymult_line_preprocess polymult_line_preprocess_sse2
#define	polymult_line		polymult_line_sse2
#define broadcastsd		_mm_set1_pd
#define addpd			_mm_add_pd
#define subpd			_mm_sub_pd
#define mulpd			_mm_mul_pd
#define divpd			_mm_div_pd
#else
#define VLEN			64			// Vector is 256 bits wide
#define VDT			double			// Vector data type
#define VFMA			0			// FMA not supported
#define	read_line		read_line_dbl		// Naming extension for this option
#define	read_preprocess_line	read_preprocess_line_dbl
#define	write_line		write_line_dbl
#define	brute_line		brute_line_dbl
#define	karatsuba_line		karatsuba_line_dbl
#define	fft_line		fft_line_dbl
#define pick_pass_sizes		pick_pass_sizes_dbl
#define	monic_line_add_hi	monic_line_add_hi_dbl
#define	monic_line_add_lo	monic_line_add_lo_dbl
#define	monic_line_adjustments	monic_line_adjustments_dbl
#define	polymult_line_preprocess polymult_line_preprocess_dbl
#define	polymult_line		polymult_line_dbl
#define broadcastsd		_mm64_set1_pd
#define addpd			_mm64_add_pd
#define subpd			_mm64_sub_pd
#define mulpd			_mm64_mul_pd
#define divpd			_mm64_div_pd
// Emulate vector length 1 using Intel-style instrinsic naming conventions
#define __m64d			double
#define _mm64_set1_pd(a)	(a)
#define _mm64_add_pd(a,b)	((a)+(b))
#define _mm64_sub_pd(a,b)	((a)-(b))
#define _mm64_mul_pd(a,b)	((a)*(b))
#define _mm64_div_pd(a,b)	((a)/(b))
#endif

// Emulate FMA on CPUs that do not have an FMA instruction
#if !VFMA
#define fmaddpd(a,b,c)	addpd(mulpd(a,b),c)
#define fmsubpd(a,b,c)	subpd(mulpd(a,b),c)
#define fnmaddpd(a,b,c)	subpd(c,mulpd(a,b))
#endif

// Structure for complex VDTs
typedef struct {
	VDT	real;
	VDT	imag;
} CVDT;
#define cvload(a,b,c)		(a).real = broadcastsd (b), (a).imag = broadcastsd (c)
#define cvrmul(a,const_d,c)	(a).real = mulpd (const_d, (c).real), (a).imag = mulpd (const_d, (c).imag)				// a = const_d*c
#define cvadd(a,b,c)		(a).real = addpd ((b).real, (c).real), (a).imag = addpd ((b).imag, (c).imag)				// a = b + c
#define cvsub(a,b,c)		(a).real = subpd ((b).real, (c).real), (a).imag = subpd ((b).imag, (c).imag)				// a = b - c
#define cviadd(a,b,c)		(a).real = subpd ((b).real, (c).imag), (a).imag = addpd ((b).imag, (c).real)				// a = b + i*c
#define cvisub(a,b,c)		(a).real = addpd ((b).real, (c).imag), (a).imag = subpd ((b).imag, (c).real)				// a = b - i*c
#define cvaddfm(a,b,const_d,c)	(a).real = fmaddpd (const_d, (c).real, (b).real), (a).imag = fmaddpd (const_d, (c).imag, (b).imag)	// a = b + const_d*c
#define cvsubfm(a,b,const_d,c)	(a).real = fnmaddpd (const_d, (c).real, (b).real), (a).imag = fnmaddpd (const_d, (c).imag, (b).imag)	// a = b - const_d*c
#define cviaddfm(a,b,const_d,c)	(a).real = fnmaddpd (const_d, (c).imag, (b).real), (a).imag = fmaddpd (const_d, (c).real, (b).imag)	// a = b + i*const_d*c
#define cvisubfm(a,b,const_d,c)	(a).real = fmaddpd (const_d, (c).imag, (b).real), (a).imag = fnmaddpd (const_d, (c).real, (b).imag)	// a = b - i*const_d*c
//#define cvmul(a,b,c)		(a).real = fmsubpd ((b).real, (c).real, mulpd ((b).imag, (c).imag)), \
//				(a).imag = fmaddpd ((b).imag, (c).real, mulpd ((b).real, (c).imag))					// a = b * c
#define cvmul(a,b,c)		{ VDT __tmpproduct__ = mulpd ((b).imag, (c).imag); \
				  (a).imag = fmaddpd ((b).imag, (c).real, mulpd ((b).real, (c).imag)); \
				  (a).real = fmsubpd ((b).real, (c).real, __tmpproduct__); }						// a = b * c
#define cvconjmul(a,b,c)	(a).real = fmaddpd ((b).real, (c).real, mulpd ((b).imag, (c).imag)), \
				(a).imag = fmsubpd ((b).imag, (c).real, mulpd ((b).real, (c).imag))					// a = b * conjugate(c)
//#define cvmulsafe(a,b,c)	{ CVDT tmp; cvmul (tmp, b, c); a = tmp; }								// a = b * c safely
#define cvmulsafe(a,b,c)	cvmul (a, b, c);											// a = b * c safely

// Structure for twiddle vector data type.
// Twiddles aren't stored as real & imag.  Instead they are stored as cos & cos/sin (a.k.a. real & real/imag).  This provides more FMA opportunities.
typedef struct {
	VDT	cos_over_sin;
	VDT	sin;
} TVDT;
#define twidload(t,b)		(t).cos_over_sin = broadcastsd ((b)[0]), (t).sin = broadcastsd ((b)[1])
#define twidmul(a,b,t)		(a).real = mulpd (fmsubpd ((b).real, (t).cos_over_sin, (b).imag), (t).sin), \
				(a).imag = mulpd (fmaddpd ((b).imag, (t).cos_over_sin, (b).real), (t).sin)		// a = b * t
#define twidconjmul(a,b,t)	(a).real = mulpd (fmaddpd ((b).real, (t).cos_over_sin, (b).imag), (t).sin), \
				(a).imag = mulpd (fmsubpd ((b).imag, (t).cos_over_sin, (b).real), (t).sin)		// a = b * conjugate(t)
#define twidmuldelay(a,b,t)	(a).real = fmsubpd ((b).real, (t).cos_over_sin, (b).imag), \
				(a).imag = fmaddpd ((b).imag, (t).cos_over_sin, (b).real)				// a = b * t, delay mul by sin
#define twidconjmuldelay(a,b,t)	(a).real = fmaddpd ((b).real, (t).cos_over_sin, (b).imag), \
				(a).imag = fmsubpd ((b).imag, (t).cos_over_sin, (b).real)				// a = b * conjugate(t), delay mul by sin

#if !defined (SSE2) && !defined (AVX) && !defined (FMA) && !defined (AVX512)

// Get the needed safety_margin required for an invec1_size by invec2_size polymult
double polymult_safety_margin (int invec1_size, int invec2_size) {
	// The smaller poly determines the number of partial products that are added together
	int	n = (invec1_size < invec2_size) ? invec1_size : invec2_size;
	// The larger poly approximates the FFT size of the polymult (assumes POLYMULT_CIRCULAR)
	int	fft_size = (invec1_size > invec2_size) ? invec1_size : invec2_size;
	// According to roe_gwnum.cpp, the number of EXTRA_BITS (output FFT bits) consumed by a polymult ~= log2(n)/2.  However, gwinit's safety margin
	// parameter refers to input FFT bits.  Thus, divide by another 2.  However, in practice the formula fails for large polymults.  P-1 with AVX-512 FFT
	// size 4608, and an RLP poly1 of size 320000 times poly2 of size 1.3M we get roundoff errors.  My theory is that while roe_gwnum.cpp accurately estimates
	// the number of most significant bits needed to safely add 640000 partial products, it is not accurately estimating the rounding errors introduced
	// to the least significant bits computing a 1.3M size FFT in fft_line.  Until we get more data points (the example got a 0.48+ roundoff where 0.4 should
	// be the target maximum acceptable error -- that's log2(0.49/0.4)=0.3 extra output bits), adjust the needed expected bits by 0.3*log2(FFTsize)/log2(1.3M).
	// 2022-02-05: AVX-512 FFT size 4608 on M76333 fails with poly1-size:322560, poly2-size:1339033.  Increasing to 0.4.
	// 2022-03-02: FMA FFT size 640 on M10303 fails with poly1-size:9123840, poly2-size:37280030, roundoff: 0.4453.  Increasing to 0.45.
	return ((log2 ((double) n) / 2.0 + 0.45 * log2 ((double) fft_size) / log2 (1300000.0)) / 2.0);
}

// Get the FFT size that will be used for an n = invec1_size + invec2_size polymult
int polymult_fft_size (int n)
{
//GW: These FFT sizes haven't been timed.  Some may be slower than a larger FFT or more radix-3/5 options might be profitable
	int	pfas[] = {40, 45, 48, 50, 54, 60, 64, 72, 75};		// Provide a variety of radix-3 and radix-5 FFT sizes

	// Power of two when n is small -- brute force is probably faster anyway
	if (n <= 32) for (int i = 4; ; i *= 2) if (i >= n) return (i);

	// Increase n by a power of two until we find the smallest fft size larger than n
	for (int two_multiplier = 1; ; two_multiplier *= 2) {
		if (n > two_multiplier * 75) continue;
		for (int i = 0; i < sizeof (pfas) / sizeof (int); i++) {
			int fft_size = two_multiplier * pfas[i];
			if (fft_size >= n && (fft_size % 4 == 0 || fft_size % 2 == 1)) return (fft_size);  // Weed out fft sizes that would require a radix-2 step
		}
	}
}

// Get the memory (in bytes) required for an FFT based polymult.  Use the information to ensure over-allocating memory does not occur.
uint64_t polymult_mem_required (
	int	invec1_size,		// Size of poly #1
	int	invec2_size,		// Size of poly #2
	int	options,		// Polymult options
	int	cpu_flags,		// CPU flags.  Vector size is determined by AVX-512, AVX, SSE2 settings.
	int	num_threads)		// Each polymult thread allocates memory
{
	int	adjusted_invec1_size, adjusted_invec2_size, adjusted_outvec_size, pass1_size, vector_size;
	bool	post_process_monics;
	uint64_t fft_size, memory;

	ASSERTG (num_threads >= 1);

	// Determine the complex vector size
	vector_size = (cpu_flags & CPU_AVX512F) ? 2*64 : (cpu_flags & CPU_AVX) ? 2*32 : (cpu_flags & CPU_SSE2) ? 2*16 : 2*8;
	
	// Adjust sizes due to RLPs.  Calculate the output size assuming no monic inputs.
	adjusted_invec1_size = (options & POLYMULT_INVEC1_RLP) ? 2 * invec1_size - 1 : invec1_size;
	adjusted_invec2_size = (options & POLYMULT_INVEC2_RLP) ? 2 * invec2_size - 1 : invec2_size;
	adjusted_outvec_size = adjusted_invec1_size + adjusted_invec2_size - 1;
	if (options & POLYMULT_CIRCULAR) {
		if (options & POLYMULT_INVEC1_MONIC) adjusted_invec1_size += (options & POLYMULT_INVEC1_RLP) ? 2 : 1;
		if (options & POLYMULT_INVEC2_MONIC) adjusted_invec2_size += (options & POLYMULT_INVEC2_RLP) ? 2 : 1;
		// Assume circular size will be the larger of the input sizes rounded up to the next poly FFT size
		adjusted_outvec_size = (adjusted_invec1_size > adjusted_invec2_size ? adjusted_invec1_size : adjusted_invec2_size);
		post_process_monics = FALSE;
	} else
		post_process_monics = TRUE;

	// Compute FFT size and memory usage
	fft_size = polymult_fft_size (adjusted_outvec_size);
	pass1_size = (int) sqrt ((double) fft_size);		// Assume pass1 size roughly equals pass2 size
	memory = fft_size * vector_size;			// FFT of invec1
	memory += (fft_size + 2) * vector_size;			// FFT of invec2/outvec
	if (post_process_monics && (options & POLYMULT_INVEC2_MONIC)) memory += fft_size * vector_size;
	if (post_process_monics && (options & POLYMULT_INVEC1_MONIC)) memory += fft_size * vector_size;
	memory += pass1_size * vector_size;

	// Each thread requires memory
	memory *= num_threads;

	// Also account for sin/cos tables
	return (memory + (fft_size % 3 == 0 ? fft_size * 2 / 3 : fft_size) * sizeof (double));
}

// Initialize a polymult handle
void polymult_init (
	pmhandle *pmdata,		// Handle for polymult library that needs initializing
	gwhandle *gwdata)		// Handle for gwnum FFT library
{
	memset (pmdata, 0, sizeof (pmhandle));
	pmdata->gwdata = gwdata;

	// Initialize default number of threads
	pmdata->max_num_threads = gwget_num_threads (gwdata);
	pmdata->num_threads = pmdata->max_num_threads;

	// Initialize default breakeven points
	if (pmdata->gwdata->cpu_flags & CPU_AVX) {
		pmdata->KARAT_BREAK = 32;			//GW: Fix me
		pmdata->FFT_BREAK = 64;				//GW: Fix me
	}
	else {
		pmdata->KARAT_BREAK = 32;			//GW: Fix me
		pmdata->FFT_BREAK = 64;				//GW: Fix me
	}

	// Default enables caching twiddles
	pmdata->cached_twiddles_enabled = TRUE;
	pmdata->twiddle_cache_additions_disabled = FALSE;
			
	// Set default cache size to optimize FFTs for
	pmdata->L2_CACHE_SIZE = 256;

	// Init FFT(1) if needed
	if (pmdata->gwdata->k != 1.0) gwuser_init_FFT1 (pmdata->gwdata);
	if (pmdata->gwdata->ZERO_PADDED_FFT) zpad7_fft (pmdata->gwdata, pmdata->gwdata->GW_FFT1);
}

// Terminate use of a polymult handle.  Free up memory.
void polymult_done (
	pmhandle *pmdata)		// Handle for polymult library
{
	if (pmdata->thread_ids != NULL) {
		pmdata->helpers_must_exit = TRUE;
		gwevent_signal (&pmdata->work_to_do);
		atomic_set (pmdata->alt_work_to_do, 1);
		for (int i = 1; i < pmdata->max_num_threads; i++) gwthread_wait_for_exit (&pmdata->thread_ids[i]);
		free (pmdata->thread_ids);
		pmdata->thread_ids = NULL;
		gwmutex_destroy (&pmdata->poly_mutex);
		gwevent_destroy (&pmdata->work_to_do);
		gwevent_destroy (&pmdata->all_helpers_done);
	}
	if (pmdata->cached_twiddles_enabled) {
		for (int i = 0; i < pmdata->cached_twiddles_count; i++) {
			aligned_free (pmdata->cached_twiddles[i].twiddles1);
			aligned_free (pmdata->cached_twiddles[i].twiddles2);
		}
		pmdata->cached_twiddles_count = 0;
	}
	if (pmdata->twiddles_initialized && !pmdata->twiddles_are_from_cache) {
		aligned_free (pmdata->twiddles1);
		aligned_free (pmdata->twiddles2);
	}
	pmdata->twiddles_initialized = 0;
	pmdata->twiddles1 = NULL;
	pmdata->twiddles2 = NULL;
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

//
// Process the 7 values in the header of zero-padded FFT gwnums.  These seven values cannot be multiplied pointwise as is required by polymult.
// The assembly code (see zmult7) multiplies the 7 values keeping the high half of the result.  An zero-padded 13-reals FFT of the 7 values allows
// us to do pointwise multiplies.  A zero-padded FFT produces 1 real value and 6 complex values for a total of 13 doubles.  The extra 6 doubles need
// to be stored somewhere in the gwnum header or a pad line in the gwnum.
//

// Roots of unity for radix-13 FFT
const long double w1c = 0.8854560256532098959003755220151;	// cos(2*pi/13*1)
const long double w1s = 0.4647231720437685456560153351331;	// sin(2*pi/13*1)
const long double w2c = 0.56806474673115580251180755912752;	// cos(2*pi/13*2)
const long double w2s = 0.82298386589365639457961742343938;	// sin(2*pi/13*2)
const long double w3c = 0.12053668025532305334906768745254;	// cos(2*pi/13*3)
const long double w3s = 0.99270887409805399280075164949252;	// sin(2*pi/13*3)
const long double w4c = -0.35460488704253562596963789260002;	// cos(2*pi/13*4)
const long double w4s = 0.93501624268541482343978459983783;	// sin(2*pi/13*4)
const long double w5c = -0.74851074817110109863463059970135;	// cos(2*pi/13*5)
const long double w5s = 0.66312265824079520237678549266677;	// sin(2*pi/13*5)
const long double w6c = -0.97094181742605202715698227629379;	// cos(2*pi/13*6)
const long double w6s = 0.23931566428755776714875372626021;	// sin(2*pi/13*6)
const long double w7c = -0.97094181742605202715698227629379;	// cos(2*pi/13*7)
const long double w7s = -0.23931566428755776714875372626021;	// sin(2*pi/13*7)
const long double w8c = -0.74851074817110109863463059970135;	// cos(2*pi/13*8)
const long double w8s = -0.66312265824079520237678549266677;	// sin(2*pi/13*8)
const long double w9c = -0.35460488704253562596963789260002;	// cos(2*pi/13*9)
const long double w9s = -0.93501624268541482343978459983783;	// sin(2*pi/13*9)
const long double w10c = 0.12053668025532305334906768745254;	// cos(2*pi/13*10)
const long double w10s = -0.99270887409805399280075164949252;	// sin(2*pi/13*10)
const long double w11c = 0.56806474673115580251180755912752;	// cos(2*pi/13*11)
const long double w11s = -0.82298386589365639457961742343938;	// sin(2*pi/13*11)
const long double w12c = 0.8854560256532098959003755220151;	// cos(2*pi/13*12)
const long double w12s = -0.4647231720437685456560153351331;	// sin(2*pi/13*12)

// FFT the 7 zpad values
void zpad7_fft (
	gwhandle *gwdata,		// Handle for gwnum library
	gwnum	g)			// gwnum to process the 7 zpad values in the header
{
	double r1, r2, r3, r4, r5, r6, r7;		// The 7 reals from the header (r7=MSW, r1=LSW)

	// Load the zpad doubles from the header
	if (gwdata->cpu_flags & CPU_AVX512F) {
		r1 = g[-11];		// Load input FFT word at halfway - 3
		r2 = g[-10];		// Load input FFT word at halfway - 2
		r3 = g[-9];		// Load input FFT word at halfway - 1
		r4 = g[-8];		// Load input FFT word at halfway + 0
		r5 = g[-7];		// Load input FFT word at halfway + 1
		r6 = g[-6];		// Load input FFT word at halfway + 2
		r7 = g[-5];		// Load input FFT word at halfway + 3
	} else {
		r1 = g[-11];		// Load input FFT word at halfway - 3
		r2 = g[-10];		// Load input FFT word at halfway - 2
		r3 = g[-9];		// Load input FFT word at halfway - 1
		r4 = g[-5];		// Load input FFT word at halfway + 0
		r5 = g[-6];		// Load input FFT word at halfway + 1
		r6 = g[-7];		// Load input FFT word at halfway + 2
		r7 = g[-8];		// Load input FFT word at halfway + 3
	}

	// Brute force a zero padded 13-point FFT.  Hermetian symetry allows this to be stored in 7 complex values.
	g[-12] = r1 + r2 + r3 + r4 + r5 + r6 + r7;
	g[-13] = r1 + w1c*r2 +  w2c*r3 +  w3c*r4 +  w4c*r5 +  w5c*r6 +  w6c*r7;
	g[-14] = r1 + w2c*r2 +  w4c*r3 +  w6c*r4 +  w8c*r5 + w10c*r6 + w12c*r7;
	g[-15] = r1 + w3c*r2 +  w6c*r3 +  w9c*r4 + w12c*r5 +  w2c*r6 +  w5c*r7;
	g[-16] = r1 + w4c*r2 +  w8c*r3 + w12c*r4 +  w3c*r5 +  w7c*r6 + w11c*r7;
	g[-17] = r1 + w5c*r2 + w10c*r3 +  w2c*r4 +  w7c*r5 + w12c*r6 +  w4c*r7;
	g[-18] = r1 + w6c*r2 + w12c*r3 +  w5c*r4 + w11c*r5 +  w4c*r6 + w10c*r7;
	g[-19] =      w1s*r2 +  w2s*r3 +  w3s*r4 +  w4s*r5 +  w5s*r6 +  w6s*r7;
	g[-20] =      w2s*r2 +  w4s*r3 +  w6s*r4 +  w8s*r5 + w10s*r6 + w12s*r7;
	g[-21] =      w3s*r2 +  w6s*r3 +  w9s*r4 + w12s*r5 +  w2s*r6 +  w5s*r7;
	g[-22] =      w4s*r2 +  w8s*r3 + w12s*r4 +  w3s*r5 +  w7s*r6 + w11s*r7;
	g[-23] =      w5s*r2 + w10s*r3 +  w2s*r4 +  w7s*r5 + w12s*r6 +  w4s*r7;
	g[-24] =      w6s*r2 + w12s*r3 +  w5s*r4 + w11s*r5 +  w4s*r6 + w10s*r7;
}

// Inverse FFT to reconstruct the 7 zpad values
void zpad7_unfft (
	gwhandle *gwdata,		// Handle for gwnum library
	gwnum	g)			// gwnum to process the 7 zpad values in the header
{
const	double INV13 = 0.07692307692307692307692307692308;	// 1/13
	double r1, r2, r3, r4, r5, r6, r7;			// The 7 FFTed reals
	double i2, i3, i4, i5, i6, i7;				// The 7 FFTed imaginaries (first imaginary is zero)
	double o[7];						// The 13 output reals (we don't need first 6 reals)

	// Load FFT of the zpad doubles
	r1 = g[-12];
	r2 = g[-13];
	r3 = g[-14];
	r4 = g[-15];
	r5 = g[-16];
	r6 = g[-17];
	r7 = g[-18];
	i2 = g[-19];
	i3 = g[-20];
	i4 = g[-21];
	i5 = g[-22];
	i6 = g[-23];
	i7 = g[-24];

	// Brute force the padded 13-point inverse FFT discarding the lower 6 values
	o[0] = r1 + 2.0 * ( w6c*r2 + w12c*r3 +  w5c*r4 + w11c*r5 +  w4c*r6 + w10c*r7 +  w6s*i2 + w12s*i3 +  w5s*i4 + w11s*i5 +  w4s*i6 + w10s*i7);
	o[1] = r1 + 2.0 * ( w7c*r2 +  w1c*r3 +  w8c*r4 +  w2c*r5 +  w9c*r6 +  w3c*r7 +  w7s*i2 +  w1s*i3 +  w8s*i4 +  w2s*i5 +  w9s*i6 +  w3s*i7);
	o[2] = r1 + 2.0 * ( w8c*r2 +  w3c*r3 + w11c*r4 +  w6c*r5 +  w1c*r6 +  w9c*r7 +  w8s*i2 +  w3s*i3 + w11s*i4 +  w6s*i5 +  w1s*i6 +  w9s*i7);
	o[3] = r1 + 2.0 * ( w9c*r2 +  w5c*r3 +  w1c*r4 + w10c*r5 +  w6c*r6 +  w2c*r7 +  w9s*i2 +  w5s*i3 +  w1s*i4 + w10s*i5 +  w6s*i6 +  w2s*i7);
	o[4] = r1 + 2.0 * (w10c*r2 +  w7c*r3 +  w4c*r4 +  w1c*r5 + w11c*r6 +  w8c*r7 + w10s*i2 +  w7s*i3 +  w4s*i4 +  w1s*i5 + w11s*i6 +  w8s*i7);
	o[5] = r1 + 2.0 * (w11c*r2 +  w9c*r3 +  w7c*r4 +  w5c*r5 +  w3c*r6 +  w1c*r7 + w11s*i2 +  w9s*i3 +  w7s*i4 +  w5s*i5 +  w3s*i6 +  w1s*i7);
	o[6] = r1 + 2.0 * (w12c*r2 + w11c*r3 + w10c*r4 +  w9c*r5 +  w8c*r6 +  w7c*r7 + w12s*i2 + w11s*i3 + w10s*i4 +  w9s*i5 +  w8s*i6 +  w7s*i7);

	// Store the zpad doubles in the header
	if (gwdata->cpu_flags & CPU_AVX512F) {
		g[-11] = o[0] * INV13;	// Store output FFT word at halfway - 3
		g[-10] = o[1] * INV13;	// Store output FFT word at halfway - 2
		g[-9]  = o[2] * INV13;	// Store output FFT word at halfway - 1
		g[-8]  = o[3] * INV13;	// Store output FFT word at halfway + 0
		g[-7]  = o[4] * INV13;	// Store output FFT word at halfway + 1
		g[-6]  = o[5] * INV13;	// Store output FFT word at halfway + 2
		g[-5]  = o[6] * INV13;	// Store output FFT word at halfway + 3
	} else {
		g[-11] = o[0] * INV13;	// Store output FFT word at halfway - 3
		g[-10] = o[1] * INV13;	// Store output FFT word at halfway - 2
		g[-9]  = o[2] * INV13;	// Store output FFT word at halfway - 1
		g[-5]  = o[3] * INV13;	// Store output FFT word at halfway + 0
		g[-6]  = o[4] * INV13;	// Store output FFT word at halfway + 1
		g[-7]  = o[5] * INV13;	// Store output FFT word at halfway + 2
		g[-8]  = o[6] * INV13;	// Store output FFT word at halfway + 3
	}
}

// Helper routines for FFT code

// Generate radix-3/4/5 twiddle factors
void generate_sincos (
	pmhandle *pmdata,		// Handle for polymult library
	const unsigned int fft_size)
{
	const long double twopi = 2.0L * 3.1415926535897932384626433L;
#define gen_twid(t,a)	{long double c,s; c=cosl(a); s=sinl(a)+1e-80; (t)[0]=(double)(c/s); (t)[1]=(double)s; }

	// Return quickly if twiddles have already been initialized
	if (pmdata->twiddles_initialized && pmdata->twiddles_initialized == fft_size) return;

	// Make sure only one thread initializes the twiddles
	if (pmdata->num_threads > 1) gwmutex_lock (&pmdata->poly_mutex);

	// If some other thread beat us to the initialization, then return
	if (pmdata->twiddles_initialized && pmdata->twiddles_initialized == fft_size) {
		if (pmdata->num_threads > 1) gwmutex_unlock (&pmdata->poly_mutex);
		return;
	}

//GW: accept fft_size < pmdata->twiddles_initialized (at least for power of 2)??

	// Free the previous twiddles unless they are cached
	if (pmdata->twiddles_initialized && !pmdata->twiddles_are_from_cache) {
		aligned_free (pmdata->twiddles1);
		aligned_free (pmdata->twiddles2);
	}

	// Look for desired twiddles in the cache
	if (pmdata->cached_twiddles_enabled) {
		for (int i = 0; i < pmdata->cached_twiddles_count; i++) {
			if (fft_size == pmdata->cached_twiddles[i].size) {
				pmdata->twiddles_initialized = fft_size;
				pmdata->twiddles1 = pmdata->cached_twiddles[i].twiddles1;
				pmdata->twiddles2 = pmdata->cached_twiddles[i].twiddles2;
				pmdata->twiddles_are_from_cache = TRUE;
				if (pmdata->num_threads > 1) gwmutex_unlock (&pmdata->poly_mutex);
				return;
			}
		}
	}

	// Clear previous state
	pmdata->twiddles_initialized = 0;
	pmdata->twiddles1 = NULL;
	pmdata->twiddles2 = NULL;

	// Generate twiddles for radix-3 butterflies
	unsigned int size = fft_size;
	if (size % 3 == 0) {
		long double twopi_over_size = twopi / (long double) size;
		pmdata->twiddles1 = (double *) aligned_malloc (size/3 * 2 * sizeof (double), 64);
		for (unsigned int i = 0; i < size/3; i++) {
			long double angle = twopi_over_size * (long double) i;
			gen_twid (&pmdata->twiddles1[2*i], angle);
		}
		while (size % 3 == 0) size /= 3;
	}

	// Generate twiddles for radix-4/5 butterflies
	long double twopi_over_size = twopi / (long double) size;
	pmdata->twiddles2 = (double *) aligned_malloc (size/4 * 2 * 2 * sizeof (double), 64);
	for (unsigned int i = 0; i < size/4; i++) {
		double long angle = twopi_over_size * (long double) i;
		gen_twid (&pmdata->twiddles2[2*(2*i)], angle);
		gen_twid (&pmdata->twiddles2[2*(2*i+1)], 2.0L*angle);
	}

	// Set size
	pmdata->twiddles_initialized = fft_size;
	pmdata->twiddles_are_from_cache = FALSE;

	// Add twiddles to the caches
	if (pmdata->cached_twiddles_enabled && !pmdata->twiddle_cache_additions_disabled) {
		int i;
		if (pmdata->cached_twiddles_count < 40) i = pmdata->cached_twiddles_count++;
		else {
			i = fft_size % 37;		// Delete random entry
			aligned_free (pmdata->cached_twiddles[i].twiddles1);
			aligned_free (pmdata->cached_twiddles[i].twiddles2);
		}
		pmdata->cached_twiddles[i].size = fft_size;
		pmdata->cached_twiddles[i].twiddles1 = pmdata->twiddles1;
		pmdata->cached_twiddles[i].twiddles2 = pmdata->twiddles2;
		pmdata->twiddles_are_from_cache = TRUE;
	}

	// Unlock and return
	if (pmdata->num_threads > 1) gwmutex_unlock (&pmdata->poly_mutex);
}

// Compress a vector of doubles
// The compression algorithm is pretty simple.  It relies on the observation that exponents in the double are clustered around a peak, dropping rapidly
// in frequency as you move away from the peak.
// The three exponents that constitute the peak are output with 2-bit values.
// Three more exponents are output with 4-bit values (either 1 above the peak and 2 below the peak or 3 below the peak).
// The next exponent below the peak is output with a 5-bit value.
// The remaining exponents are output with a 5-bit escape code followed by the 11-bit exponent (many 11-bit exponents are "impossible" and we can
// likely reduce this to a 5-bit escape code followed by a 6-bit exponent).
int compress_vec (		// Returns size of the compressed buffer
	char	*buf,		// Array of doubles to compress
	int	size)		// Size of buffer in bytes
{
	short	counts[2048];
	int	peak;			// The highest exponent in the set of three most common exponents
	bool	above_peak_flag;	// TRUE if exponent above the peak is compressed

//char *orig_buf = malloc (size);
//memcpy(orig_buf, buf, size);

	// Init the output queues
	unsigned char *outptr = (unsigned char *) buf;
	int	outval = 0;		// Compressed exponent and sign bits output in multiples of 8 bits
	int	outlen = 0;		// Number of bits in outval
	unsigned char mbuf[48];		// Mantissa buffer
	int	mlen = 0;		// Number of chars in mbuf
	int	last_m48len = 0;	// Number of chars in the last mantissa48 (either 0 or 6)
#define put_bits(v,len)	outval = (outval << (len)) + (int)(v), outlen += len
#define put_m48(m)	{ union {uint64_t a; char b[8];} x; x.a = m; memcpy (mbuf + mlen, x.b, 6); mlen += 6; last_m48len = 6; }
#define emptyQ()	{ if (outlen>=8) outlen-=8, *outptr++ = (unsigned char)(outval>>outlen); \
			  if (mlen > last_m48len) mlen-=last_m48len, memcpy(outptr,mbuf,mlen), outptr+=mlen, memcpy(mbuf,mbuf+mlen,last_m48len), mlen=last_m48len; \
			  while (outlen>=8) outlen-=8, *outptr++ = (unsigned char)(outval>>outlen); \
			  if (outlen == 0 && mlen) memcpy(outptr,mbuf,mlen), outptr+=mlen, mlen=0; \
			  last_m48len = 0; }

	// Sample 1000 doubles to find the most common exponents in the doubles
//GW: count only expos between 960 an 1080  (a smaller array)
	memset (&counts, 0, sizeof (counts));
	for (int i = 0; i < size && i < 8000; i += 8) {
		uint64_t expo;
		expo = * (uint64_t *) (buf + i);
		expo = (expo << 1) >> 53;	// strip off sign bit and mantissa
		counts[expo]++;
	}

	// Find the set of three exponents that are most common
	peak = 960;
	for (int i = 961; i < 1080; i++)
		if (counts[i-2] + counts[i-1]  + counts[i] > counts[peak-2] + counts[peak-1] + counts[peak]) peak = i;
	ASSERTG (peak != 960);

	// Determine expo just above the peak will be compressed (or output raw)
	above_peak_flag = (counts[peak+1] > counts[peak-5]);

	// Output header bits
	put_bits (peak - 960, 7);
	put_bits (above_peak_flag, 1);

	// Loop compressing doubles
	for (int i = 0; i < size; i += 8) {
		uint64_t expo, sign, nibble, mantissa48;

		// Parse the double
		expo = * (uint64_t *) (buf + i);
		sign = expo >> 63;
		nibble = (expo << 12) >> 60;
		mantissa48 = expo & 0xFFFFFFFFFFFFULL;
		expo = (expo << 1) >> 53;	// strip off sign bit and mantissa

		// Empty the output queues (as much as possible)
		emptyQ ();

		// Rare case of compressing FFTed data that contained lots of zeros.  The FFTed data may contain lots of really small exponent values (noise)
		// that would greatly dimish our compression.  Turn the noise into zeros.
		if (expo < peak - 50) expo = 0;

		// Compress the exponent
		if (expo == peak) put_bits (0, 2);
		else if (expo == peak-1) put_bits (1, 2);
		else if (expo == peak-2) put_bits (2, 2);
		else if (expo == peak-3) put_bits (12, 4);
		else if (expo == peak-4) put_bits (13, 4);
		else if (above_peak_flag && expo == peak+1) put_bits (14, 4);
		else if (above_peak_flag && expo == peak-5) put_bits (30, 5);
		else if (!above_peak_flag && expo == peak-5) put_bits (14, 4);
		else if (!above_peak_flag && expo == peak-6) put_bits (30, 5);
		else {
			put_bits (31, 5);
			ASSERTG (expo == 0 || (expo >= peak-50 && expo <= peak-6) || (expo >= peak+1 && expo <= peak+18));
			if (expo > peak) expo = expo - (peak+1) + 1;		// Handle up to 18 values above peak coded as 1 through 18
			else if (expo > 0) expo = expo - (peak-50) + 19;	// Handle up to 45 values below peak coded as 19 through 63
			ASSERTG (expo <= 63);
			put_bits (expo, 6);
		}

		// Output the sign bit, high nibble of mantissa, and rest of the mantissa
		if (expo != 0) { put_bits (sign, 1); put_bits (nibble, 4); put_m48 (mantissa48); }
	}

	// Empty the output queues
	if (outlen & 7) put_bits (0, 8 - (outlen & 7));
	emptyQ ();

//char *xxx = malloc (size);
//decompress_vec (buf, xxx, size);
//int foo;
//foo=memcmp(orig_buf, xxx, size);
//free(xxx);
//free(orig_buf);

	// Return size of the compressed vector
	return ((int) ((char *) outptr - buf));

#undef put_bits
#undef put_m48
#undef emptyQ
}

// Decompress a vector of doubles
void decompress_vec (
	char	*inbuf,		// Compressed array of doubles
	char	*outbuf,	// Decompressed array of doubles
	int	size)		// Size of output buffer in bytes
{
	unsigned int inval, inbits;	// Queue to extract bits from inbuf
	int	peak;			// The highest exponent in the set of three most common exponents
	bool	above_peak_flag;	// TRUE if exponent above the peak is compressed

	// Init input queue
	unsigned char *inptr = (unsigned char *) inbuf;
	inval = inbits = 0;
#define get_bits(v,len)	{ if (inbits<len) inval=(inval<<8)+(*inptr++),inbits+=8; inbits-=len; v=(inval>>inbits)&((1<<(len))-1); }
#define get_m48(v)	{ union {uint64_t a; char b[8];} x; x.a = 0; memcpy (x.b, inptr, 6); inptr += 6; v = x.a; }

	// Read peak exponent and above-peak-flag
	get_bits (peak, 7);  peak += 960;
	get_bits (above_peak_flag, 1);

#ifdef STRAIGHTFORWARD_IMPL
	while (size) {
		uint64_t expo, sign, nibble, mantissa48;
		int	code;

		// Decompress the exponent
		get_bits (code, 2);
		if (code <= 2) expo = peak - code;
		else {
			get_bits (code, 2);
			if (code <= 1) expo = peak - 3 - code;
			else if (code == 2) expo = above_peak_flag ? peak + 1 : peak - 5;
			else {
				get_bits (code, 1);
				if (code == 0) expo = above_peak_flag ? peak - 5 : peak - 6;
				else {
					get_bits (code, 6);
					if (code >= 19) expo = code - 19 + (peak-50);	// Handle up to 45 values below peak coded as 19 through 63
					else if (code >= 1) expo = code - 1 + (peak+1);	// Handle up to 18 values above peak coded as 1 through 18
					else expo = 0;
				}
			}
		}

		// Get the sign bit and mantissa
		if (expo != 0) { get_bits (sign, 1); get_bits (nibble, 4); get_m48 (mantissa48); }
		else sign = nibble = mantissa48 = 0;

		// Output the uncompressed double, move onto next double
		* (uint64_t *) outbuf = (sign << 63) + (expo << 52) + (nibble << 48) + mantissa48;
		outbuf += sizeof (double);
		size -= sizeof (double);
	}

#else		// This version is less readable but a little faster

	for ( ; size; outbuf += sizeof (double), size -= sizeof (double)) {
		int	code, more_code, expo, sign, nibble;

		// We need to get at least 7 bits, cheaper to get them all at once.
		// The goal is to decompress ~75% of the expos using just one get_bits call.
		get_bits (code, 7);
		
		// Decompress the exponent
		if (code < 0b1100000) expo = peak - (code >> 5);
		else if (code < 0b1111000) {
			if (code < 0b1110000) expo = peak - 3 - ((code & 0b0011000) >> 3);
			else expo = above_peak_flag ? peak + 1 : peak - 5;
			get_bits (more_code, 2);
			code = (code << 2) + more_code;
		} else if (code < 0b1111100) {
			expo = above_peak_flag ? peak - 5 : peak - 6;
			get_bits (more_code, 3);
			code = (code << 3) + more_code;
		} else {
			get_bits (more_code, 4);
			code = ((code & 0b0000011) << 4) + more_code;
			if (code >= 19) expo = code - 19 + (peak-50);	// Handle up to 45 values below peak coded as 19 through 63
			else if (code >= 1) expo = code - 1 + (peak+1);	// Handle up to 18 values above peak coded as 1 through 18
			else { memset (outbuf, 0, sizeof (double)); continue; }
			get_bits (code, 5);
		}

		// Get the sign bit and high nibble of mantissa
		sign = code >> 4;
		nibble = code & 0xF;

		// Copy 6 mantissa bytes
		memcpy (outbuf, inptr, 6), inptr += 6;

		// Output the top two bytes of the double
		uint16_t top2 = (uint16_t) ((sign << 15) + (expo << 4) + nibble);
		memcpy (outbuf+6, &top2, 2);
	}
#endif

#undef get_bits
#undef get_m48
}

#endif

//--------------------------------------------------------------------
//		Read/write cache lines from/to a gwnum
//    Requires intimate knowledge of layout of gwnums in memory
//--------------------------------------------------------------------

#if defined (AVX512)

void read_line_avx512 (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics)
{
	int	orig_index = index;			// Save for a possible recursive call to read_line

	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
		if (gwdata->k == 1.0) {
			vec[0].real = broadcastsd (1.0);
			vec[0].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, orig_index, vec, 0, FALSE);
		}
		vec++;
	}

	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Check if we've read all the gwnum data
	ASSERTG (index * 16 < (int) gwdata->FFTLEN);

	// AVX512 implementation
	if (index == -1) {			// First real-only value and optionally seven zero-pad values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			// Return the first real-only value
			if (!gwdata->ZERO_PADDED_FFT) {
				vec[j].real = _mm512_set_pd (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, data[j][0]);
				vec[j].imag = broadcastsd (0.0);
			}
			// Zero padded FFTs also return the 7 FFTed values for the 7 real values stored in the header
			else {
				vec[j].real = _mm512_set_pd (data[j][-18], data[j][-17], data[j][-16], data[j][-15],
							     data[j][-14], data[j][-13], data[j][-12], data[j][0]);
				vec[j].imag = _mm512_set_pd (data[j][-24], data[j][-23], data[j][-22], data[j][-21],
							     data[j][-20], data[j][-19], 0.0, 0.0);
			}
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {	// Second real-only value and 7 complex values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			vec[j].real = _mm512_set_pd (data[j][7], data[j][6], data[j][5], data[j][4], data[j][3], data[j][2], data[j][1], data[j][8]);
			vec[j].imag = _mm512_set_pd (data[j][15], data[j][14], data[j][13], data[j][12], data[j][11], data[j][10], data[j][9], 0.0);
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index * 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof(CVDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			vec[j].real = ((VDT *) data[j])[offset];
			vec[j].imag = ((VDT *) data[j])[offset+1];
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		VDT zero = broadcastsd (0.0);
		for (int j = 0; j < size; j++) {
			vec[j].real = subpd (zero, vec[j].real);
			vec[j].imag = subpd (zero, vec[j].imag);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) vec[-j] = vec[j];
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
		if (gwdata->k == 1.0) {
			vec[size].real = broadcastsd (1.0);
			vec[size].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, orig_index, &vec[size], 0, FALSE);
		}
	}
}

void write_line_avx512 (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, gwnum *fma_data, int options)
{
	CVDT	fmaval, oneval;

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// AVX512 implementation
	union {
		VDT	a;
		double	b[8];
	} x;
	if (index == -1) {			// First real-only value and optionally seven zero-pad values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL)
					fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else
					fmaval.real = _mm512_set_pd (fma_data[j][-18], fma_data[j][-17], fma_data[j][-16], fma_data[j][-15],
								     fma_data[j][-14], fma_data[j][-13], fma_data[j][-12], fma_data[j][0]),
					fmaval.imag = _mm512_set_pd (fma_data[j][-24], fma_data[j][-23], fma_data[j][-22], fma_data[j][-21],
								     fma_data[j][-20], fma_data[j][-19], 0.0, 0.0);
				if (gwdata->k != 1.0) {
					oneval.real = _mm512_set_pd (gwdata->GW_FFT1[-18], gwdata->GW_FFT1[-17], gwdata->GW_FFT1[-16], gwdata->GW_FFT1[-15],
								     gwdata->GW_FFT1[-14], gwdata->GW_FFT1[-13], gwdata->GW_FFT1[-12], gwdata->GW_FFT1[0]);
					oneval.imag = _mm512_set_pd (gwdata->GW_FFT1[-24], gwdata->GW_FFT1[-23], gwdata->GW_FFT1[-22], gwdata->GW_FFT1[-21],
								     gwdata->GW_FFT1[-20], gwdata->GW_FFT1[-19], 0.0, 0.0);
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			// Return the first real-only value
			x.a = vec[j].real;
			data[j][0] = x.b[0];
			// Zero padded FFTs also return the 7 real values stored in the header
			if (gwdata->ZERO_PADDED_FFT) {
				data[j][-18] = x.b[7];
				data[j][-17] = x.b[6];
				data[j][-16] = x.b[5];
				data[j][-15] = x.b[4];
				data[j][-14] = x.b[3];
				data[j][-13] = x.b[2];
				data[j][-12] = x.b[1];
				x.a = vec[j].imag;
				data[j][-24] = x.b[7];
				data[j][-23] = x.b[6];
				data[j][-22] = x.b[5];
				data[j][-21] = x.b[4];
				data[j][-20] = x.b[3];
				data[j][-19] = x.b[2];
			}
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {	// Second real-only value and 7 complex values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL)
					fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else
					fmaval.real = _mm512_set_pd (fma_data[j][7], fma_data[j][6], fma_data[j][5], fma_data[j][4],
								     fma_data[j][3], fma_data[j][2], fma_data[j][1], fma_data[j][8]),
					fmaval.imag = _mm512_set_pd (fma_data[j][15], fma_data[j][14], fma_data[j][13], fma_data[j][12],
								     fma_data[j][11], fma_data[j][10], fma_data[j][9], 0.0);
				if (gwdata->k != 1.0) {
					oneval.real = _mm512_set_pd (gwdata->GW_FFT1[7], gwdata->GW_FFT1[6], gwdata->GW_FFT1[5], gwdata->GW_FFT1[4],
								     gwdata->GW_FFT1[3], gwdata->GW_FFT1[2], gwdata->GW_FFT1[1], gwdata->GW_FFT1[8]);
					oneval.imag = _mm512_set_pd (gwdata->GW_FFT1[15], gwdata->GW_FFT1[14], gwdata->GW_FFT1[13], gwdata->GW_FFT1[12],
								     gwdata->GW_FFT1[11], gwdata->GW_FFT1[10], gwdata->GW_FFT1[9], 0.0);
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			x.a = vec[j].real, data[j][7] = x.b[7], data[j][6] = x.b[6], data[j][5] = x.b[5], data[j][4] = x.b[4],
					   data[j][3] = x.b[3], data[j][2] = x.b[2], data[j][1] = x.b[1], data[j][8] = x.b[0];
			x.a = vec[j].imag, data[j][15] = x.b[7], data[j][14] = x.b[6], data[j][13] = x.b[5], data[j][12] = x.b[4],
					   data[j][11] = x.b[3], data[j][10] = x.b[2], data[j][9] = x.b[1];
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index * 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL) fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else fmaval.real = ((VDT *) fma_data[j])[offset], fmaval.imag = ((VDT *) fma_data[j])[offset+1];
				if (gwdata->k != 1.0) {
					oneval.real = ((VDT *) gwdata->GW_FFT1)[offset], oneval.imag = ((VDT *) gwdata->GW_FFT1)[offset+1];
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			((VDT *) data[j])[offset] = vec[j].real;
			((VDT *) data[j])[offset+1] = vec[j].imag;
		}
	}
}

#elif defined (FMA)				// FMA and AVX FFTs use the same memory layout

#define read_line_fma		read_line_avx
#define read_preprocess_line_fma read_preprocess_line_avx
#define write_line_fma		write_line_avx
#define pick_pass_sizes_fma	pick_pass_sizes_avx
void read_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics);
void read_preprocess_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics);
void write_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, gwnum *fma_data, int options);
void pick_pass_sizes_avx (pmhandle *pmdata, unsigned int fftsize, unsigned int *p1size, unsigned int *p2size);

#elif defined (AVX)

void read_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics)
{
	int	orig_index = index;			// Save for a possible recursive call to read_line

	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
		if (gwdata->k == 1.0) {
			vec[0].real = broadcastsd (1.0);
			vec[0].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, orig_index, vec, 0, FALSE);
		}
		vec++;
	}

	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs return another cache line for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index--;

	// Check if we've read all the gwnum data
	ASSERTG (index * 8 < (int) gwdata->FFTLEN);

	// AVX implementation
	if (index == -2) {			// Four of the seven FFTed zero pad header values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			// Zero padded FFTs, return 4 of the FFTed values for the 7 real values stored in the header
			vec[j].real = _mm256_set_pd (data[j][-18], data[j][-17], data[j][-16], data[j][-15]);
			vec[j].imag = _mm256_set_pd (data[j][-24], data[j][-23], data[j][-22], data[j][-21]);
		}
	}
	else if (index == -1) {			// First real-only value and optionally three of the seven FFTed zero-pad values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			// Return the first real-only value
			if (!gwdata->ZERO_PADDED_FFT) {
				vec[j].real = _mm256_set_pd (0.0, 0.0, 0.0, data[j][0]);
				vec[j].imag = broadcastsd (0.0);
			}
			// Zero padded FFTs also return 3 of the FFTed values for the 7 real values stored in the header
			else {
				vec[j].real = _mm256_set_pd (data[j][-14], data[j][-13], data[j][-12], data[j][0]);
				vec[j].imag = _mm256_set_pd (data[j][-20], data[j][-19], 0.0, 0.0);
			}
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {	// Second real-only value and 3 complex values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			vec[j].real = _mm256_set_pd (data[j][3], data[j][2], data[j][1], data[j][4]);
			vec[j].imag = _mm256_set_pd (data[j][7], data[j][6], data[j][5], 0.0);
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			vec[j].real = ((VDT *) data[j])[offset];
			vec[j].imag = ((VDT *) data[j])[offset+1];
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		VDT zero = broadcastsd (0.0);
		for (int j = 0; j < size; j++) {
			vec[j].real = subpd (zero, vec[j].real);
			vec[j].imag = subpd (zero, vec[j].imag);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) vec[-j] = vec[j];
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
		if (gwdata->k == 1.0) {
			vec[size].real = broadcastsd (1.0);
			vec[size].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, orig_index, &vec[size], 0, FALSE);
		}
	}
}

void write_line_avx (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, gwnum *fma_data, int options)
{
	CVDT	fmaval, oneval;

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs return another cache line for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index--;

	// AVX implementation
	union {
		VDT	a;
		double	b[4];
	} x;
	if (index == -2) {			// Four of the seven zero pad value in the header
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL)
					fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else
					fmaval.real = _mm256_set_pd (fma_data[j][-18], fma_data[j][-17], fma_data[j][-16], fma_data[j][-15]),
					fmaval.imag = _mm256_set_pd (fma_data[j][-24], fma_data[j][-23], fma_data[j][-22], fma_data[j][-21]);
				if (gwdata->k != 1.0) {
					oneval.real = _mm256_set_pd (gwdata->GW_FFT1[-18], gwdata->GW_FFT1[-17], gwdata->GW_FFT1[-16], gwdata->GW_FFT1[-15]),
					oneval.imag = _mm256_set_pd (gwdata->GW_FFT1[-24], gwdata->GW_FFT1[-23], gwdata->GW_FFT1[-22], gwdata->GW_FFT1[-21]);
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			x.a = vec[j].real;
			data[j][-18] = x.b[3];
			data[j][-17] = x.b[2];
			data[j][-16] = x.b[1];
			data[j][-15] = x.b[0];
			x.a = vec[j].imag;
			data[j][-24] = x.b[3];
			data[j][-23] = x.b[2];
			data[j][-22] = x.b[1];
			data[j][-21] = x.b[0];
		}
	}
	else if (index == -1) {			// First real-only value and optionally three of the seven zero-pad values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL)
					fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else
					fmaval.real = _mm256_set_pd (fma_data[j][-14], fma_data[j][-13], fma_data[j][-12], fma_data[j][0]),
					fmaval.imag = _mm256_set_pd (fma_data[j][-20], fma_data[j][-19], 0.0, 0.0);
				if (gwdata->k != 1.0) {
					oneval.real = _mm256_set_pd (gwdata->GW_FFT1[-14], gwdata->GW_FFT1[-13], gwdata->GW_FFT1[-12], gwdata->GW_FFT1[0]),
					oneval.imag = _mm256_set_pd (gwdata->GW_FFT1[-20], gwdata->GW_FFT1[-19], 0.0, 0.0);
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			// Return the first real-only value
			x.a = vec[j].real;
			data[j][0] = x.b[0];
			// Return 3 of the 7 real values stored in the header
			if (gwdata->ZERO_PADDED_FFT) {
				data[j][-14] = x.b[3];
				data[j][-13] = x.b[2];
				data[j][-12] = x.b[1];
				x.a = vec[j].imag;
				data[j][-20] = x.b[3];
				data[j][-19] = x.b[2];
			}
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL)
					fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else
					fmaval.real = _mm256_set_pd (fma_data[j][3], fma_data[j][2], fma_data[j][1], fma_data[j][4]),
					fmaval.imag = _mm256_set_pd (fma_data[j][7], fma_data[j][6], fma_data[j][5], 0.0);
				if (gwdata->k != 1.0) {
					oneval.real = _mm256_set_pd (gwdata->GW_FFT1[3], gwdata->GW_FFT1[2], gwdata->GW_FFT1[1], gwdata->GW_FFT1[4]),
					oneval.imag = _mm256_set_pd (gwdata->GW_FFT1[7], gwdata->GW_FFT1[6], gwdata->GW_FFT1[5], 0.0);
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			x.a = vec[j].real, data[j][3] = x.b[3], data[j][2] = x.b[2], data[j][1] = x.b[1], data[j][4] = x.b[0];
			x.a = vec[j].imag, data[j][7] = x.b[3], data[j][6] = x.b[2], data[j][5] = x.b[1];
		}
	}
	else {				// Process complex values where imaginary part is in the high half of the cache line
		unsigned long offset = cache_line_offset (gwdata, index) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL) fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else fmaval.real = ((VDT *) fma_data[j])[offset], fmaval.imag = ((VDT *) fma_data[j])[offset+1];
				if (gwdata->k != 1.0) {
					oneval.real = ((VDT *) gwdata->GW_FFT1)[offset], oneval.imag = ((VDT *) gwdata->GW_FFT1)[offset+1];
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			((VDT *) data[j])[offset] = vec[j].real;
			((VDT *) data[j])[offset+1] = vec[j].imag;
		}
	}
}

#elif defined (SSE2)

void read_line_sse2 (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics)
{
	int	orig_index = index;			// Save for a possible recursive call to read_line

	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
		if (gwdata->k == 1.0) {
			vec[0].real = broadcastsd (1.0);
			vec[0].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, orig_index, vec, 0, FALSE);
		}
		vec++;
	}

	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs return another 3 cache lines for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 3;

	// Check if we've read all the gwnum data
	ASSERTG (index * 4 < (int) gwdata->FFTLEN);

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

	// SSE2 implementation
	if (index <= -2) {			// Six of the seven zero pad header values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			// Zero padded FFTs, return 1 of the FFTed values for the 7 real values stored in the header
			vec[j].real = _mm_set_pd (data[j][index*2-10], data[j][index*2-9]);
			vec[j].imag = _mm_set_pd (data[j][index*2-16], data[j][index*2-15]);
		}
	}
	else if (index == -1) {			// First real-only valueand optionally one of the seven zero-pad values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			// Return the first real-only value
			if (!gwdata->ZERO_PADDED_FFT) {
				vec[j].real = _mm_set_pd (0.0, data[j][0]);
				vec[j].imag = broadcastsd (0.0);
			}
			// Zero padded FFTs also return 1 of the FFTed values for the 7 real values stored in the header
			else {
				vec[j].real = _mm_set_pd (data[j][-12], data[j][0]);
				vec[j].imag = broadcastsd (0.0);
			}
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {	// Second real-only value and one complex value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			if (imag_dist == 8) {
				vec[j].real = _mm_set_pd (data[j][2], data[j][1]);
				vec[j].imag = _mm_set_pd (data[j][3], 0.0);
			} else if (imag_dist == 16) {
				vec[j].real = _mm_set_pd (data[j][1], data[j][2]);
				vec[j].imag = _mm_set_pd (data[j][3], 0.0);
			} else {
				vec[j].real = _mm_set_pd (data[j][1], data[j][4]);
				vec[j].imag = _mm_set_pd (data[j][5], 0.0);
			}
		}
	}
	else {				// Process complex values where imaginary part is 8, 16, or 32 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
			if (imag_dist == 8) {
				vec[j].real = _mm_set_pd (data[j][index*4+2], data[j][index*4]);
				vec[j].imag = _mm_set_pd (data[j][index*4+3], data[j][index*4+1]);
			} else if (imag_dist == 16) {
				vec[j].real = ((VDT *) data[j])[offset+(index&1)*2];
				vec[j].imag = ((VDT *) data[j])[offset+(index&1)*2+1];
			} else {
				vec[j].real = ((VDT *) data[j])[offset+(index&1)];
				vec[j].imag = ((VDT *) data[j])[offset+(index&1)+2];
			}
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		VDT zero = broadcastsd (0.0);
		for (int j = 0; j < size; j++) {
			vec[j].real = subpd (zero, vec[j].real);
			vec[j].imag = subpd (zero, vec[j].imag);
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) vec[-j] = vec[j];
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
		if (gwdata->k == 1.0) {
			vec[size].real = broadcastsd (1.0);
			vec[size].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, orig_index, &vec[size], 0, FALSE);
		}
	}
}

void write_line_sse2 (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, gwnum *fma_data, int options)
{
	CVDT	fmaval, oneval;

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs also return 3 extra cache lines for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 3;

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
	if (index <= -2) {		// Six of the seven real values in the header
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL)
					fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else
					fmaval.real = _mm_set_pd (fma_data[j][index*2-10], fma_data[j][index*2-9]),
					fmaval.imag = _mm_set_pd (fma_data[j][index*2-16], fma_data[j][index*2-15]);
				if (gwdata->k != 1.0) {
					oneval.real = _mm_set_pd (gwdata->GW_FFT1[index*2-10], gwdata->GW_FFT1[index*2-9]),
					oneval.imag = _mm_set_pd (gwdata->GW_FFT1[index*2-16], gwdata->GW_FFT1[index*2-15]);
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			x.a = vec[j].real;
			data[j][index*2-10] = x.b[1];
			data[j][index*2-9] = x.b[0];
			x.a = vec[j].imag;
			data[j][index*2-16] = x.b[1];
			data[j][index*2-15] = x.b[0];
		}
	}
	else if (index == -1) {		// First real-only value and optionally one of the seven zero-pad values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL) fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else fmaval.real = _mm_set_pd (fma_data[j][-12], fma_data[j][0]);
				if (gwdata->k != 1.0) {
					oneval.real = _mm_set_pd (gwdata->GW_FFT1[-12], gwdata->GW_FFT1[0]);
					fmaval.real = mulpd (fmaval.real, oneval.real);
				}
				if (options & POLYMULT_FMADD) vec[j].real = addpd (vec[j].real, fmaval.real);
				else if (options & POLYMULT_FMSUB) vec[j].real = subpd (vec[j].real, fmaval.real);
				else vec[j].real = subpd (fmaval.real, vec[j].real);
			}
			// Return the first real-only value
			x.a = vec[j].real;
			data[j][0] = x.b[0];
			// Return 1 of the 7 real values stored in the header
			if (gwdata->ZERO_PADDED_FFT) data[j][-12] = x.b[1];
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {		// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL)
					fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else if (imag_dist == 8)
					fmaval.real = _mm_set_pd (fma_data[j][2], fma_data[j][1]), fmaval.imag = _mm_set_pd (fma_data[j][3], 0.0);
				else if (imag_dist == 16)
					fmaval.real = _mm_set_pd (fma_data[j][1], fma_data[j][2]), fmaval.imag = _mm_set_pd (fma_data[j][3], 0.0);
				else
					fmaval.real = _mm_set_pd (fma_data[j][1], fma_data[j][4]), fmaval.imag = _mm_set_pd (fma_data[j][5], 0.0);
				if (gwdata->k != 1.0) {
					if (imag_dist == 8)
						oneval.real = _mm_set_pd (gwdata->GW_FFT1[2], gwdata->GW_FFT1[1]), oneval.imag = _mm_set_pd (gwdata->GW_FFT1[3], 0.0);
					else if (imag_dist == 16)
						oneval.real = _mm_set_pd (gwdata->GW_FFT1[1], gwdata->GW_FFT1[2]), oneval.imag = _mm_set_pd (gwdata->GW_FFT1[3], 0.0);
					else
						oneval.real = _mm_set_pd (gwdata->GW_FFT1[1], gwdata->GW_FFT1[4]), oneval.imag = _mm_set_pd (gwdata->GW_FFT1[5], 0.0);
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			if (imag_dist == 8) {
				x.a = vec[j].real, data[j][2] = x.b[1], data[j][1] = x.b[0];
				x.a = vec[j].imag, data[j][3] = x.b[1];
			} else if (imag_dist == 16) {
				x.a = vec[j].real, data[j][1] = x.b[1], data[j][2] = x.b[0];
				x.a = vec[j].imag, data[j][3] = x.b[1];
			} else {
				x.a = vec[j].real, data[j][1] = x.b[1], data[j][4] = x.b[0];
				x.a = vec[j].imag, data[j][5] = x.b[1];
			}
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 2) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL)
					fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else if (imag_dist == 8)
					fmaval.real = _mm_set_pd (fma_data[j][index*4+2], fma_data[j][index*4]),
					fmaval.imag = _mm_set_pd (fma_data[j][index*4+3], fma_data[j][index*4+1]);
				else if (imag_dist == 16)
					fmaval.real = ((VDT *) fma_data[j])[offset+(index&1)*2],
					fmaval.imag = ((VDT *) fma_data[j])[offset+(index&1)*2+1];
				else 
					fmaval.real = ((VDT *) fma_data[j])[offset+(index&1)],
					fmaval.imag = ((VDT *) fma_data[j])[offset+(index&1)+2];
				if (gwdata->k != 1.0) {
					if (imag_dist == 8)
						oneval.real = _mm_set_pd (gwdata->GW_FFT1[index*4+2], gwdata->GW_FFT1[index*4]),
						oneval.imag = _mm_set_pd (gwdata->GW_FFT1[index*4+3], gwdata->GW_FFT1[index*4+1]);
					else if (imag_dist == 16)
						oneval.real = ((VDT *) gwdata->GW_FFT1)[offset+(index&1)*2],
						oneval.imag = ((VDT *) gwdata->GW_FFT1)[offset+(index&1)*2+1];
					else 
						oneval.real = ((VDT *) gwdata->GW_FFT1)[offset+(index&1)],
						oneval.imag = ((VDT *) gwdata->GW_FFT1)[offset+(index&1)+2];
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			if (imag_dist == 8) {
				x.a = vec[j].real, data[j][index*4+2] = x.b[1], data[j][index*4] = x.b[0];
				x.a = vec[j].imag, data[j][index*4+3] = x.b[1], data[j][index*4+1] = x.b[0];
			} else if (imag_dist == 16) {
				((VDT *) data[j])[offset+(index&1)*2] = vec[j].real;
				((VDT *) data[j])[offset+(index&1)*2+1] = vec[j].imag;
			} else {
				((VDT *) data[j])[offset+(index&1)] = vec[j].real;
				((VDT *) data[j])[offset+(index&1)+2] = vec[j].imag;
			}
		}
	}
}

#elif !defined (X86_64)				// x87 FFTs only supported in 32-bit mode

void read_line_dbl (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics)
{
	int	orig_index = index;			// Save for a possible recursive call to read_line

	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
		if (gwdata->k == 1.0) {
			vec[0].real = broadcastsd (1.0);
			vec[0].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, orig_index, vec, 0, FALSE);
		}
		vec++;
	}

	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += (size-1);

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;

	// Check if we've read all the gwnum data
	ASSERTG (index * 2 < (int) gwdata->FFTLEN);

	// No special instructions implementation
	if (index <= -2) {			// Seven zero pad header values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			// Zero padded FFTs return 1 of the FFTed values for the 7 real values stored in the header
			vec[j].real = data[j][index-10];
			vec[j].imag = data[j][index-17];
		}
	}
	else if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			vec[j].real = data[j][0];
			vec[j].imag = 0.0;
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {			// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			if (gwdata->PASS2_SIZE == 0) vec[j].real = data[j][1];
			else vec[j].real = data[j][2];
			vec[j].imag = 0.0;
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 4) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) { memset (&vec[j], 0, sizeof (CVDT)); continue; }
			ASSERTG (FFT_state (data[j]) == FULLY_FFTed);
// bug???  Index assuming what?  skip every other set of 4?  Examine ancient FFT sources!
			if (gwdata->PASS2_SIZE == 0 && index < 8 && !gwdata->ALL_COMPLEX_FFT) {
				vec[j].real = data[j][offset+(index&3)*2];
				vec[j].imag = data[j][offset+(index&3)*2+1];
			} else {
				vec[j].real = data[j][offset+(index&2)*2+(index&1)];
				vec[j].imag = data[j][offset+(index&2)*2+(index&1)+2];
			}
		}
	}

	// Negate coefficients if requested
	if (options & (POLYMULT_INVEC1_NEGATE | POLYMULT_INVEC2_NEGATE)) {
		for (int j = 0; j < size; j++) {
			vec[j].real = -vec[j].real;
			vec[j].imag = -vec[j].imag;
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) vec[-j] = vec[j];
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
		if (gwdata->k == 1.0) {
			vec[size].real = broadcastsd (1.0);
			vec[size].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, orig_index, &vec[size], 0, FALSE);
		}
	}
}

void write_line_dbl (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, gwnum *fma_data, int options)
{
	CVDT	fmaval, oneval;

	// Transformed real data ends up with 2 real values and FFTLEN/2-1 complex values.  Transformed complex data has FFTLEN/2 complex values.
	// Adjust index to handle the extra data point that must be returned for transformed real data.
	if (!gwdata->ALL_COMPLEX_FFT) index--;

	// Zero padded FFTs must also return for processing the 7 real values stored in the header
	if (gwdata->ZERO_PADDED_FFT) index -= 7;

	// No special instructions implementation
	if (index <= -2) {			// Seven zero pad header values
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL) fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else fmaval.real = fma_data[j][index-10], fmaval.imag = fma_data[j][index-17];
				if (gwdata->k != 1.0) {
					oneval.real = gwdata->GW_FFT1[index-10], oneval.imag = gwdata->GW_FFT1[index-17];
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			data[j][index-10] = vec[j].real;
			data[j][index-17] = vec[j].imag;
		}
	}
	else if (index == -1) {			// First real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL) fmaval.real = broadcastsd (0.0);
				else fmaval.real = fma_data[j][0];
				if (gwdata->k != 1.0) {
					oneval.real = gwdata->GW_FFT1[0];
					fmaval.real = mulpd (fmaval.real, oneval.real);
				}
				if (options & POLYMULT_FMADD) vec[j].real = vec[j].real + fmaval.real;
				else if (options & POLYMULT_FMSUB) vec[j].real = vec[j].real - fmaval.real;
				else vec[j].real = fmaval.real - vec[j].real;
			}
			data[j][0] = vec[j].real;
		}
	}
	else if (index == 0 && !gwdata->ALL_COMPLEX_FFT) {		// Second real-only value
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL) fmaval.real = broadcastsd (0.0);
				else if (gwdata->PASS2_SIZE == 0) fmaval.real = fma_data[j][1];
				else fmaval.real = fma_data[j][2];
				if (gwdata->k != 1.0) {
					oneval.real = gwdata->GW_FFT1[2];
					fmaval.real = mulpd (fmaval.real, oneval.real);
				}
				if (options & POLYMULT_FMADD) vec[j].real = vec[j].real + fmaval.real;
				else if (options & POLYMULT_FMSUB) vec[j].real = vec[j].real - fmaval.real;
				else vec[j].real = fmaval.real - vec[j].real;
			}
			if (gwdata->PASS2_SIZE == 0) data[j][1] = vec[j].real;
			else data[j][2] = vec[j].real;
		}
	}
	else {				// Process complex values where imaginary part is 8 or 16 bytes above the real part
		unsigned long offset = cache_line_offset (gwdata, index / 4) / sizeof (VDT);
		for (int j = 0; j < size; j++) {
			if (data[j] == NULL) continue;
			if (options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) {
				if (fma_data == NULL || fma_data[j] == NULL) fmaval.real = broadcastsd (0.0), fmaval.imag = broadcastsd (0.0);
				else if (gwdata->PASS2_SIZE == 0 && index < 8 && !gwdata->ALL_COMPLEX_FFT)
					fmaval.real = fma_data[j][offset+(index&3)*2], fmaval.imag = fma_data[j][offset+(index&3)*2+1];
				else
					fmaval.real = fma_data[j][offset+(index&2)*2+(index&1)], fmaval.imag = fma_data[j][offset+(index&2)*2+(index&1)+2];
				if (gwdata->k != 1.0) {
					if (gwdata->PASS2_SIZE == 0 && index < 8 && !gwdata->ALL_COMPLEX_FFT)
						oneval.real = gwdata->GW_FFT1[offset+(index&3)*2], oneval.imag = gwdata->GW_FFT1[offset+(index&3)*2+1];
					else
						oneval.real = gwdata->GW_FFT1[offset+(index&2)*2+(index&1)], oneval.imag = gwdata->GW_FFT1[offset+(index&2)*2+(index&1)+2];
					cvmulsafe (fmaval, fmaval, oneval);
				}
				if (options & POLYMULT_FMADD) cvadd (vec[j], vec[j], fmaval);
				else if (options & POLYMULT_FMSUB) cvsub (vec[j], vec[j], fmaval);
				else cvsub (vec[j], fmaval, vec[j]);
			}
			if (gwdata->PASS2_SIZE == 0 && index < 8 && !gwdata->ALL_COMPLEX_FFT) {
				data[j][offset+(index&3)*2] = vec[j].real;
				data[j][offset+(index&3)*2+1] = vec[j].imag;
			} else {
				data[j][offset+(index&2)*2+(index&1)] = vec[j].real;
				data[j][offset+(index&2)*2+(index&1)+2] = vec[j].imag;
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
	CVDT	*invec1,		// First input poly (may be a monic with ones stripped)
	int	invec1_size,		// Size of the first input polynomial
	CVDT	*invec2,		// Second input poly (may be a monic with ones stripped)
	int	invec2_size,		// Size of the second input polynomial
	CVDT	*outvec,		// Output poly
	int	circular_size,		// Output poly size (may be less than invec1_size + invec2_size - 1)
	gwnum	*gw_outvec,		// Gwnum array passed to polymult (so we can detect outputs that do not need to be computed)
	int	gw_outvec_size,		// Size of the gwnum output vector (may be less than output poly size)
	int	adjusted_shift,		// Difference between base of outvec and base of true_outvec.  Can be non-zero when monic RLPs are multiplied.
	int	LSWs_skipped)		// Least significant coefficients from true_outvec_size that do not need to be returned
{
	// Loop through output indices
	for (int out = 0; out < invec1_size + invec2_size - 1; out++) {
		// Adjust output index for circular convolutions
		int adjusted_out = out + adjusted_shift;
		bool first_time = TRUE;
		if (adjusted_out >= circular_size) {
			adjusted_out = adjusted_out % circular_size;
			first_time = FALSE;
		}
		// See if this output needs to be computed
		if (gw_outvec != NULL &&
		    (adjusted_out < LSWs_skipped || adjusted_out >= LSWs_skipped + gw_outvec_size || gw_outvec[adjusted_out - LSWs_skipped] == NULL)) continue;
		// Get index into input vectors
		int in1 = out < invec2_size ? 0 : out - invec2_size + 1;	// Index into invec1
		int in2 = out - in1;						// Index into invec2
		// Loop through all the pairs to accumulate in this output
		for ( ; in1 < invec1_size && in2 >= 0; in1++, in2--) {
			CVDT tmp;
			cvmul (tmp, invec1[in1], invec2[in2]);
			if (first_time) {
				outvec[adjusted_out] = tmp;
				first_time = FALSE;
			} else {
				cvadd (outvec[adjusted_out], outvec[adjusted_out], tmp);
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
	CVDT	*invec1,		// Coefficients of first input poly
	int	invec1_size,		// Size of the first input polynomial
	CVDT	*invec2,		// Coefficients of second input poly
	int	invec2_size,		// Size of the second input polynomial
	CVDT	*tmp,			// Temporary space
	CVDT	*outvec)		// Output poly coefficients
{
	int	i, a_size, b_size, c_size, d_size, a_plus_b_size, c_plus_d_size, bd_size, ac_size, a_plus_b_times_c_plus_d_size;
	CVDT	*a, *b, *c, *d, *ac, *bd, *a_plus_b, *c_plus_d, *a_plus_b_times_c_plus_d;

	// Handle end of recursion
	if (invec1_size == 1 || invec2_size == 1 || invec1_size + invec2_size < pmdata->KARAT_BREAK) {
		brute_line (invec1, invec1_size, invec2, invec2_size, outvec, invec1_size + invec2_size - 1, NULL, 0, 0, 0);
		return;
	}

	// To simplify code, always make invec2 the larger poly.
	if (invec1_size > invec2_size) {
		int tmpsize = invec1_size; invec1_size = invec2_size; invec2_size = tmpsize;
		CVDT *tmpvec = invec1; invec1 = invec2; invec2 = tmpvec;
	}

	// Split larger vector roughly in half.  For readability, create pointers to two halves.
	d_size = (invec2_size + 1) / 2;
	c_size = invec2_size - d_size;
	c = invec2 + d_size, d = invec2;

	// Handle case where vectors are of wildly unequal lengths
	if (d_size >= invec1_size) {
		// Compute bd
		bd = outvec;
		karatsuba_line (pmdata, invec1, invec1_size, d, d_size, tmp, bd);
		bd_size = invec1_size + d_size - 1;
		// Compute bc and store it in tmp
		CVDT *bc = tmp;
		int bc_size = invec1_size + c_size - 1;
		tmp += bc_size;
		karatsuba_line (pmdata, invec1, invec1_size, c, c_size, tmp, bc);
		// add/copy bc to output
		for (i = 0; i < invec1_size - 1; i++) cvadd (outvec[i + d_size], outvec[i + d_size], bc[i]);
		for ( ; i < bc_size; i++) outvec[i + d_size] = bc[i];
		return;
	}

	// Split smaller vector, the size of the two lower halves must be the same.  For readability, create pointers to two halves.
	b_size = d_size;
	a_size = invec1_size - b_size;
	a = invec1 + b_size, b = invec1;

	// Compute (a+b)
	a_plus_b = outvec;
	for (i = 0; i < a_size && i < b_size; i++) cvadd (a_plus_b[i], a[i], b[i]);
//	for ( ; i < a_size; i++) a_plus_b[i] = a[i];		// In case a_size > b_size
	for ( ; i < b_size; i++) a_plus_b[i] = b[i];		// In case b_size > a_size
	a_plus_b_size = i;
	// Compute (c+d)
	c_plus_d = &outvec[i];
	for (i = 0; i < c_size && i < d_size; i++) cvadd (c_plus_d[i], c[i], d[i]);
//	for ( ; i < c_size; i++) c_plus_d[i] = c[i];		// In case c_size > d_size
	for ( ; i < d_size; i++) c_plus_d[i] = d[i];		// In case d_size > c_size
	c_plus_d_size = i;
	// Compute (a+b) * (c+d) and store it in tmp
	a_plus_b_times_c_plus_d = tmp;
	a_plus_b_times_c_plus_d_size = a_plus_b_size + c_plus_d_size - 1;
	tmp += a_plus_b_times_c_plus_d_size;
	karatsuba_line (pmdata, a_plus_b, a_plus_b_size, c_plus_d, c_plus_d_size, tmp, a_plus_b_times_c_plus_d);

	// Compute bd
	bd = outvec;
	karatsuba_line (pmdata, b, b_size, d, d_size, tmp, bd);
	bd_size = b_size + d_size - 1;
	// Output the complex value 0+0i
	memset (&outvec[bd_size], 0, sizeof (CVDT));
	// Compute ac
	ac = &outvec[bd_size+1];
	karatsuba_line (pmdata, a, a_size, c, c_size, tmp, ac);
	ac_size = a_size + c_size - 1;

	// subtract ac from middle
	for (i = 0; i < ac_size; i++) cvsub (a_plus_b_times_c_plus_d[i], a_plus_b_times_c_plus_d[i], ac[i]);
	// subtract bd from middle
	for (i = 0; i < bd_size; i++) cvsub (a_plus_b_times_c_plus_d[i], a_plus_b_times_c_plus_d[i], bd[i]);
	// add middle to output
	for (i = 0; i < a_plus_b_times_c_plus_d_size; i++) cvadd (outvec[i + b_size], outvec[i + b_size], a_plus_b_times_c_plus_d[i]);
}

// Pick FFT pass sizes
#ifndef FMA
void pick_pass_sizes (
	pmhandle *pmdata,		// Handle for polymult library
	unsigned int fftsize,		// FFT size
	unsigned int *p1size,		// Returned pass 1 size
	unsigned int *p2size)		// Returned pass 2 size
{
	unsigned int powers_of_two;
	double	max_p2size, log_fftsize;

	// Figure out what power-of-two FFT size will fit in the L2 cache.  One complex number fits in one CVDT.
	max_p2size = log2 ((double) pmdata->L2_CACHE_SIZE * 1024.0 / (double) sizeof (CVDT));

	// If the FFT can be done in one pass, choose that option
	log_fftsize = log2 ((double) fftsize);
	if (log_fftsize <= max_p2size) { *p1size = 1; *p2size = fftsize; return; }

	// Make pass sizes roughly equal
	*p2size = (unsigned int) log_fftsize / 2;
	if (*p2size > (unsigned int) max_p2size) *p2size = (unsigned int) max_p2size;

	// Make sure any odd power of two is done in pass 2
	for (powers_of_two = 0; (fftsize & (1 << powers_of_two)) == 0; powers_of_two++);
	if (*p2size > powers_of_two) *p2size = powers_of_two;
	if ((powers_of_two & 1) == 0) {		// no radix-2/8 step
		if ((*p2size & 1) == 1) *p2size -= 1;
	} else {
		if ((*p2size & 1) == 0) *p2size -= 1;
	}
	*p1size = fftsize >> *p2size;
	*p2size = 1 << *p2size;
}
#endif

// Perform forward FFT, pointmul, inverse FFT in as few passes over memory as possible.
// This routine can be called several different ways:
#define	FORWARD_FFT_INVEC1	0x1	// Forward FFT invec1 required
#define	FORWARD_FFT_INVEC2	0x2	// Forward FFT invec2 required
#define FORWARD_FFT_ONLY	0x4	// No pointwise multiply and inverse FFT

void fft_line (
	pmhandle *pmdata,		// Handle for polymult library
	CVDT	*invec1,		// Array of complex values
	CVDT	*invec2,		// Array of complex values
	CVDT	*outvec,		// Array of complex values
	CVDT	*scratch,		// Scratch area for two-pass FFTs       
	unsigned int fft_size,		// The fft size
	int	options)
{
	unsigned int pass1_size, pass2_size;	// Parameters for creating a two-pass FFT
	unsigned int twiddle3_stride = 1;	// Stride through the radix-3 twiddles data
	unsigned int twiddle45_stride = 1;	// Stride through the radix-4/5 twiddles data
	unsigned int inv_twiddle3_stride;	// Starting stride for the inverse FFT radix-3 twiddles data
	unsigned int inv_twiddle45_stride;	// Starting stride for the inverse FFT radix-4/5 twiddles data
	unsigned int powers_of_two;
	VDT	inv_fft_size = broadcastsd (1.0 / (double) fft_size);
	VDT	inv_fft_size_707 = broadcastsd (0.70710678118654752440084436210485 / (double) fft_size);
	bool	inv_fft_applied;		// Set to TRUE when all inverse FFT data has been multiplied by inv_fft_size

	generate_sincos (pmdata, fft_size);

	// Init some needed items
	for (powers_of_two = 0; (fft_size & (1 << powers_of_two)) == 0; powers_of_two++);		// Inverse FFT check for a radix-8 step
	for (inv_twiddle3_stride = 1, inv_twiddle45_stride = fft_size; inv_twiddle45_stride % 3 == 0; ) {// Calc first inverse FFT twiddle strides
		inv_twiddle3_stride *= 3;
		inv_twiddle45_stride /= 3;
	}

	// When the data set gets large it overflows the L2 cache.  This causes every radix-step to read & write data from main memory.  To combat
	// this we use two passes to complete the FFT, doing as much work as possible while FFT data is in the L2 cache.
	pick_pass_sizes (pmdata, fft_size, &pass1_size, &pass2_size);

	// Do two-pass forward FFT merged with a pointmul and a two-pass inverse FFT.
	// Pass 1 is first half of the forward FFT.  Pass 2 is second half of the forward FFT, pointwise multiply, and half of the inverse FFT.
	// Pass 3 is the remaining half of the inverse FFT.
	for (unsigned int pass = 1; pass <= 3; pass++) {
	    unsigned int pass_size;			// The size of this pass
	    unsigned int pass_stride;			// Distance between smallest elements in the pass
	    unsigned int num_interior_blks;		// Num blocks smaller than pass_stride (for 2-pass FFTs this equals pass2_size or 1)
	    unsigned int num_exterior_blks;		// Num blocks larger than pass_stride (for 2-pass FFTs this equals 1 or pass1_size)
	    unsigned int saved_twiddle3_stride = twiddle3_stride;
	    unsigned int saved_twiddle45_stride = twiddle45_stride;

	    // Skip pass 1 and 3 for small FFTs (perform a one-pass FFT)
	    if (pass1_size == 1 && pass != 2) continue;

	    // Init some pass-specific variables
	    if (pass == 1 || pass == 3) {
		pass_size = pass1_size;
		pass_stride = pass2_size;
		num_interior_blks = pass2_size;
		num_exterior_blks = 1;
	    } else {
		pass_size = pass2_size;
		pass_stride = 1;
		num_interior_blks = 1;
		num_exterior_blks = pass1_size;
	    }

	    // Do pass an exterior a block at a time, then do pass an interior a block at a time (we're prepared for future three-pass FFTs)
	    for (unsigned int blke = 0; blke < num_exterior_blks; blke++) {
	    for (unsigned int blki = 0; blki < num_interior_blks; blki++) {
		CVDT *src;				// FFT source data
		CVDT *dest;				// FFT destination data
		unsigned int size;			// The current group size to process
		unsigned int stride;			// Logical stride through the FFT (often the same as input and output stride through memory)
		unsigned int instride1;			// Unit stride through input FFT data in memory
		unsigned int outstride1;		// Unit stride through output FFT data in memory

// Forward FFT.  Loop over both input vectors.  Both may need a forward FFT.
		
		if (pass <= 2) {
		    for (int srcarg = 1; srcarg <= 2; srcarg++) {

			// Skip input vectors that do not need a forward FFT
			if (srcarg == 1 && !(options & FORWARD_FFT_INVEC1)) continue;
			if (srcarg == 2 && !(options & FORWARD_FFT_INVEC2)) continue;
//GW: If neither are to be FFTed, need to skip to pointmul

			// Pass 1 reads from an input vector into a scratch area to reduce large strides which CPU caches hate.  At end, result written
			// back to input vector.  Pass 2 reads/writes directly from/to an input vector.
			if (pass == 1) {
			    src = (srcarg == 1 ? invec1 : invec2) + blki;
			    dest = scratch;
			    instride1 = pass2_size;
			    outstride1 = 1;
			}
			else {
			    src = dest = (srcarg == 1 ? invec1 : invec2) + blke * pass2_size;
			    instride1 = outstride1 = 1;
			}

			// Forward FFT.  Start working with group size as the whole pass size, then work on progressively smaller group sizes.
			size = pass_size;

			// Reset twiddle strides each time we execute the pass
			twiddle3_stride = saved_twiddle3_stride;
			twiddle45_stride = saved_twiddle45_stride;

//GW: Assume all radix-3 and radix-5 work is in the first (or only) pass???

			// Radix-3 FFT
			// A 3-complex FFT is:
			// Res1:  (R1+R2+R3) + (I1+I2+I3)i
			// Res2:  (R1-.5R2-.5R3-.866I2+.866I3) + (I1-.5I2-.5I3+.866R2-.866R3)i
			// Res3:  (R1-.5R2-.5R3+.866I2-.866I3) + (I1-.5I2-.5I3-.866R2+.866R3)i
			for ( ; size % 3 == 0; size /= 3, twiddle3_stride *= 3, src = dest, instride1 = 1) {
			    stride = size / 3;
			    if (pass == 1 && stride == 1) dest = (srcarg == 1 ? invec1 : invec2) + blki, outstride1 = pass2_size;
			    VDT half = broadcastsd (0.5);
			    VDT _866 = broadcastsd (0.86602540378443864676372317075294);
			    for (unsigned int group = 0; group < pass_size; group += stride * 3) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle3_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];

					CVDT	b_plus_c, b_minus_c, tmp23;
					cvadd (b_plus_c, b, c);			// b + c
					cvsub (b_minus_c, b, c);		// b - c

					cvsubfm (tmp23, a, half, b_plus_c);			// tmp23 = a - 0.5 * (b + c)

					if (twiddle_idx == 0) {	// First group member's twiddles are 1+0i, skip twiddle multiply
#if !VFMA
						cvrmul (b_minus_c, _866, b_minus_c);
						cviadd (b, tmp23, b_minus_c);			// b = tmp23 + i * .866 * (b - c)
						cvisub (c, tmp23, b_minus_c);			// c = tmp23 - i * .866 * (b - c)
#else
						cviaddfm (b, tmp23, _866, b_minus_c);		// b = tmp23 + i * .866 * (b - c)
						cvisubfm (c, tmp23, _866, b_minus_c);		// c = tmp23 - i * .866 * (b - c)
#endif
						cvadd (a, a, b_plus_c);				// a + b + c
						dest[outstride1*(group + group_member)] = a;
						dest[outstride1*(group + group_member + stride)]   = b;
						dest[outstride1*(group + group_member + 2*stride)] = c;
					}
					else {			// Apply twiddles
						TVDT twiddle1; twidload (twiddle1, &pmdata->twiddles1[2*twiddle_idx]);
#if !VFMA
						cvrmul (b_minus_c, _866, b_minus_c);
						cvadd (a, a, b_plus_c);				// a + b + c
						cviadd (b, tmp23, b_minus_c);			// b = tmp23 + i * .866 * (b - c)
						cvisub (c, tmp23, b_minus_c);			// c = tmp23 - i * .866 * (b - c)
						dest[outstride1*(group + group_member)] = a;
						twidmul     (dest[outstride1*(group + group_member + stride)], b, twiddle1);
						twidconjmul (dest[outstride1*(group + group_member + 2*stride)], c, twiddle1);
#else
						VDT sin866 = mulpd (_866, twiddle1.sin);
						cvrmul (tmp23, twiddle1.sin, tmp23);
						cvadd (a, a, b_plus_c);				// a + b + c
						cviaddfm (b, tmp23, sin866, b_minus_c);		// b = tmp23 + i * .866 * (b - c)
						cvisubfm (c, tmp23, sin866, b_minus_c);		// c = tmp23 - i * .866 * (b - c)
						dest[outstride1*(group + group_member)] = a;
						twidmuldelay     (dest[outstride1*(group + group_member + stride)], b, twiddle1);
						twidconjmuldelay (dest[outstride1*(group + group_member + 2*stride)], c, twiddle1);
#endif
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
			for ( ; size % 5 == 0; size /= 5, twiddle45_stride *= 5, src = dest, instride1 = 1) {
			    VDT _309 = broadcastsd (0.30901699437494742410229341718282);
			    VDT _809 = broadcastsd (0.80901699437494742410229341718282);
			    VDT _951 = broadcastsd (0.95105651629515357211643933337938);
			    VDT _588_951 = broadcastsd (0.61803398874989484820458683436564);
			    stride = size / 5;
			    if (pass == 1 && stride == 1) dest = (srcarg == 1 ? invec1 : invec2) + blki, outstride1 = pass2_size;
			    for (unsigned int group = 0; group < pass_size; group += stride * 5) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle45_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];
					CVDT d = src[instride1*(group + group_member + 3*stride)];
					CVDT e = src[instride1*(group + group_member + 4*stride)];

					CVDT b_plus_e, b_minus_e, c_plus_d, c_minus_d, tmp25a, tmp25b, tmp34a, tmp34b;

					cvadd (b_plus_e, b, e);					// b + e
					cvsub (b_minus_e, b, e);				// b - e
					cvadd (c_plus_d, c, d);					// c + d
					cvsub (c_minus_d, c, d);				// c - d

					cvaddfm (tmp25a, a, _309, b_plus_e);
					cvsubfm (tmp34a, a, _809, b_plus_e);
					cvadd   (a, a, b_plus_e);
					cvaddfm (tmp25b, b_minus_e, _588_951, c_minus_d);	// tmp25b = delayed .951*(b-e) + 0.588*(c-d)
					cvsubfm (tmp34b, c_minus_d, _588_951, b_minus_e);	// tmp34b = delayed -.588*(b-e) + 0.951*(c-d)
					cvsubfm (tmp25a, tmp25a, _809, c_plus_d);		// tmp25a = a +.309*(b+e) -.809*(c+d)
					cvaddfm (tmp34a, tmp34a, _309, c_plus_d);		// tmp34a = a -.809*(b+e) +.309*(c+d)
					cvadd   (a, a, c_plus_d);				// a + b + c + d + e

					if (twiddle_idx == 0) {	// First group member's twiddles are 1+0i, skip twiddle multiply
#if !VFMA
						cvrmul (tmp25b, _951, tmp25b);			// tmp25b *= .951
						cvrmul (tmp34b, _951, tmp34b);			// tmp34b *= .951
						cviadd (b, tmp25a, tmp25b);			// b = tmp25a + i*tmp25b
						cvisub (e, tmp25a, tmp25b);			// e = tmp25a - i*tmp25b
						cvisub (c, tmp34a, tmp34b);			// c = tmp34a + i*tmp34b
						cviadd (d, tmp34a, tmp34b);			// d = tmp34a - i*tmp34b
#else
						cviaddfm (b, tmp25a, _951, tmp25b);		// b = tmp25a + i*.951*tmp25b
						cvisubfm (e, tmp25a, _951, tmp25b);		// e = tmp25a - i*.951*tmp25b
						cvisubfm (c, tmp34a, _951, tmp34b);		// c = tmp34a + i*.951*tmp34b
						cviaddfm (d, tmp34a, _951, tmp34b);		// d = tmp34a - i*.951*tmp34b
#endif
						dest[outstride1*(group + group_member)] = a;
						dest[outstride1*(group + group_member + stride)]   = b;
						dest[outstride1*(group + group_member + 2*stride)] = c;
						dest[outstride1*(group + group_member + 3*stride)] = d;
						dest[outstride1*(group + group_member + 4*stride)] = e;
					}
					else {			// Apply twiddles
						TVDT twiddle1; twidload (twiddle1, &pmdata->twiddles2[4*twiddle_idx]);
						TVDT twiddle2; twidload (twiddle2, &pmdata->twiddles2[4*twiddle_idx+2]);
#if !VFMA
						cvrmul (tmp25b, _951, tmp25b);			// tmp25b *= .951
						cvrmul (tmp34b, _951, tmp34b);			// tmp34b *= .951
						cviadd (b, tmp25a, tmp25b);			// b = tmp25a + i*tmp25b
						cvisub (e, tmp25a, tmp25b);			// e = tmp25a - i*tmp25b
						cvisub (c, tmp34a, tmp34b);			// c = tmp34a + i*tmp34b
						cviadd (d, tmp34a, tmp34b);			// d = tmp34a - i*tmp34b
						dest[outstride1*(group + group_member)] = a;
						twidmul     (dest[outstride1*(group + group_member + stride)], b, twiddle1);
						twidmul     (dest[outstride1*(group + group_member + 2*stride)], c, twiddle2);
						twidconjmul (dest[outstride1*(group + group_member + 3*stride)], d, twiddle2);
						twidconjmul (dest[outstride1*(group + group_member + 4*stride)], e, twiddle1);
#else
						VDT sin1_951 = mulpd (_951, twiddle1.sin);
						VDT sin2_951 = mulpd (_951, twiddle2.sin);
						cvrmul   (tmp25a, twiddle1.sin, tmp25a);
						cvrmul   (tmp34a, twiddle2.sin, tmp34a);
						cviaddfm (b, tmp25a, sin1_951, tmp25b);		// b = tmp25a + i*.951*tmp25b
						cvisubfm (e, tmp25a, sin1_951, tmp25b);		// e = tmp25a - i*.951*tmp25b
						cvisubfm (c, tmp34a, sin2_951, tmp34b);		// c = tmp34a + i*.951*tmp34b
						cviaddfm (d, tmp34a, sin2_951, tmp34b);		// d = tmp34a - i*.951*tmp34b
						dest[outstride1*(group + group_member)] = a;
						twidmuldelay     (dest[outstride1*(group + group_member + stride)], b, twiddle1);
						twidmuldelay     (dest[outstride1*(group + group_member + 2*stride)], c, twiddle2);
						twidconjmuldelay (dest[outstride1*(group + group_member + 3*stride)], d, twiddle2);
						twidconjmuldelay (dest[outstride1*(group + group_member + 4*stride)], e, twiddle1);
#endif
					}
				}
			    }
			}

			// Radix-4 FFT
			for ( ; size != 8 && size % 4 == 0; size /= 4, twiddle45_stride *= 4, src = dest, instride1 = 1) {
			    stride = size / 4;
			    if (pass == 1 && stride == 1) dest = (srcarg == 1 ? invec1 : invec2) + blki, outstride1 = pass2_size;
			    for (unsigned int group = 0; group < pass_size; group += stride * 4) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle45_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];
					CVDT d = src[instride1*(group + group_member + 3*stride)];

					CVDT new_c; cvsub (new_c, a, c);			// a - c
					CVDT new_d; cvsub (new_d, b, d);			// b - d
					CVDT new_a; cvadd (new_a, a, c);			// a + c
					CVDT new_b; cvadd (new_b, b, d);			// b + d

					if (twiddle_idx == 0) {	// First group member's twiddles are 1+0i, skip twiddle multiply
						cviadd (c, new_c, new_d);			// c + id
						cvisub (d, new_c, new_d);			// c - id
						cvadd  (a, new_a, new_b);			// a + b
						cvsub  (b, new_a, new_b);			// a - b
						dest[outstride1*(group + group_member)] = a;
						dest[outstride1*(group + group_member + 1*stride)] = c;		// Do not bit-reverse outputs
						dest[outstride1*(group + group_member + 2*stride)] = b;
						dest[outstride1*(group + group_member + 3*stride)] = d;
					}
					else {			// Apply twiddles
						TVDT twiddle1; twidload (twiddle1, &pmdata->twiddles2[4*twiddle_idx]);
						TVDT twiddle2; twidload (twiddle2, &pmdata->twiddles2[4*twiddle_idx+2]);
#if !VFMA
						cviadd (c, new_c, new_d);			// c + id
						cvisub (d, new_c, new_d);			// c - id
						cvadd  (a, new_a, new_b);			// a + b
						cvsub  (b, new_a, new_b);			// a - b
						dest[outstride1*(group + group_member)] = a;
						twidmul     (dest[outstride1*(group + group_member + 1*stride)], c, twiddle1);
						twidmul     (dest[outstride1*(group + group_member + 2*stride)], b, twiddle2);
						twidconjmul (dest[outstride1*(group + group_member + 3*stride)], d, twiddle1);
#else
						cvrmul   (new_c, twiddle1.sin, new_c);
						cvadd  (a, new_a, new_b);			// a + b
						cvsub  (b, new_a, new_b);			// a - b
						cviaddfm (c, new_c, twiddle1.sin, new_d);	// c + id
						cvisubfm (d, new_c, twiddle1.sin, new_d);	// c - id
						dest[outstride1*(group + group_member)] = a;
						twidmuldelay     (dest[outstride1*(group + group_member + 1*stride)], c, twiddle1);
						twidmul          (dest[outstride1*(group + group_member + 2*stride)], b, twiddle2);
						twidconjmuldelay (dest[outstride1*(group + group_member + 3*stride)], d, twiddle1);
#endif
					}
				}
			    }
			}

			// Radix-8 FFT (if needed it is the last radix block done in pass 2)
			if (size == 8) {
			    VDT _707 = broadcastsd (0.70710678118654752440084436210485);
			    stride = 1;
			    for (unsigned int group = 0; group < pass_size; group += stride * 8) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					CVDT a = src[group + group_member];
					CVDT b = src[group + group_member + stride];
					CVDT c = src[group + group_member + 2*stride];
					CVDT d = src[group + group_member + 3*stride];
					CVDT e = src[group + group_member + 4*stride];
					CVDT f = src[group + group_member + 5*stride];
					CVDT g = src[group + group_member + 6*stride];
					CVDT h = src[group + group_member + 7*stride];

					CVDT new_a, new_b, new_c, new_d, new_e, new_f, new_g, new_h, twid_f, twid_h;
					cvadd (new_a, a, e);					// a + e
					cvsub (new_e, a, e);					// a - e
					cvadd (new_b, b, f);					// b + f
					cvsub (new_f, b, f);					// b - f
					cvadd (new_c, c, g);					// c + g
					cvsub (new_g, c, g);					// c - g
					cvadd (new_d, d, h);					// d + h
					cvsub (new_h, d, h);					// d - h

					cvadd  (a, new_a, new_c);				// a + c
					cvsub  (c, new_a, new_c);				// a - c
					cvadd  (b, new_b, new_d);				// b + d
					cvsub  (d, new_b, new_d);				// b - d
					cviadd (e, new_e, new_g);				// e + ig
					cvisub (g, new_e, new_g);				// e - ig
					cviadd (f, new_f, new_h);				// f + ih
					cvisub (h, new_f, new_h);				// f - ih

					cvadd  (new_a, a, b);					// a + b
					cvsub  (new_b, a, b);					// a - b
					cviadd (new_c, c, d);					// c + id
					cvisub (new_d, c, d);					// c - id
					twid_f.real = subpd (f.real, f.imag);			// twiddled f / .707
					twid_f.imag = addpd (f.real, f.imag);
					twid_h.real = addpd (h.real, h.imag);			// -twiddled h / .707
					twid_h.imag = subpd (h.imag, h.real);
#if !VFMA
					cvrmul (twid_f, _707, twid_f);				// twiddled f
					cvrmul (twid_h, _707, twid_h);				// -twiddled h
					cvadd (new_e, e, twid_f);				// e + (sqrthalf+sqrthalf*i)f
					cvsub (new_f, e, twid_f);				// e - (sqrthalf+sqrthalf*i)f
					cvsub (new_g, g, twid_h);				// g + (-sqrthalf+sqrthalf*i)h
					cvadd (new_h, g, twid_h);				// g - (-sqrthalf+sqrthalf*i)h
#else
					cvaddfm (new_e, e, _707, twid_f);			// e + (sqrthalf+sqrthalf*i)f
					cvsubfm (new_f, e, _707, twid_f);			// e - (sqrthalf+sqrthalf*i)f
					cvsubfm (new_g, g, _707, twid_h);			// g + (-sqrthalf+sqrthalf*i)h
					cvaddfm (new_h, g, _707, twid_h);			// g - (-sqrthalf+sqrthalf*i)h
#endif

					dest[group + group_member] = new_a;			// Do not bit reverse outputs
					dest[group + group_member + 4*stride] = new_b;
					dest[group + group_member + 2*stride] = new_c;
					dest[group + group_member + 6*stride] = new_d;
					dest[group + group_member + 1*stride] = new_e;
					dest[group + group_member + 5*stride] = new_f;
					dest[group + group_member + 3*stride] = new_g;
					dest[group + group_member + 7*stride] = new_h;
				}
			    }
			}
		    }
		}

// If pass 1 then there is no pointmul and inverse FFT to do.  If only FFTing the input(s) also skip pointmul and inverse FFT.

		if (pass == 1 || (options & FORWARD_FFT_ONLY)) continue;

// Pointwise multiply

		if (pass == 2) {
		    for (unsigned int j = 0; j < pass_size; j++) {
			CVDT x = invec1[blke * pass2_size + j];
			CVDT y = invec2[blke * pass2_size + j];
			cvmul (outvec[blke * pass2_size + j], x, y);
		    }
		    inv_fft_applied = FALSE;
		    if (powers_of_two == 0) {
			    for (unsigned int j = 0; j < pass_size; j++) cvrmul (outvec[blke * pass2_size + j], inv_fft_size, outvec[blke * pass2_size + j]);
			    inv_fft_applied = TRUE;
		    }
		}

// Inverse FFT
		
		// Pass 2 reads/writes directly from/to output vector.
		// Pass 3 reads from output vector into a scratch area to reduce large strides.  At end, result written back to output vector.
		if (pass == 2) {
			src = dest = outvec + blke * pass2_size;
			instride1 = outstride1 = 1;
		} else {
			src = outvec + blki;
			dest = scratch;
			instride1 = pass2_size;
			outstride1 = 1;
		}

		// Inverse FFT starts with stride = 1 and moves to progressive smaller group sizes with larger strides.
		stride = 1;

		// In pass 2, use the starting twiddle strides.  In pass 3, restore the saved twiddle strides from the end of pass 2.
		if (pass == 2) {
			twiddle3_stride = inv_twiddle3_stride;
			twiddle45_stride = inv_twiddle45_stride;
		} else {
			twiddle3_stride = saved_twiddle3_stride;
			twiddle45_stride = saved_twiddle45_stride;
		}

		// Radix-8 FFT (if needed it is the last radix block done in pass 2)
		if (pass == 2 && powers_of_two & 1) {		// If odd number of powers of two do a radix-8 step
			twiddle45_stride /= 8;		// Next twiddle stride amount
			for (unsigned int group = 0; group < pass_size; group += stride * 8) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					CVDT a = src[group + group_member];
					CVDT b = src[group + group_member + stride];
					CVDT c = src[group + group_member + 2*stride];
					CVDT d = src[group + group_member + 3*stride];
					CVDT e = src[group + group_member + 4*stride];
					CVDT f = src[group + group_member + 5*stride];
					CVDT g = src[group + group_member + 6*stride];
					CVDT h = src[group + group_member + 7*stride];

					CVDT new_a, new_b, new_c, new_d, new_e, new_f, new_g, new_h, twid_f, twid_h;

					// Do the distance = 4 round, mul by inv_fft_size as efficiently as possible
					cvrmul (a, inv_fft_size, a);
#if !VFMA
					cvrmul (e, inv_fft_size, e);
					cvadd (new_a, a, e);					// a + e
					cvsub (new_e, a, e);					// a - e
#else
					cvaddfm (new_a, a, inv_fft_size, e);			// a + e
					cvsubfm (new_e, a, inv_fft_size, e);			// a - e
#endif
					cvadd (new_b, b, f);					// b + f
					cvsub (new_f, b, f);					// b - f
					cvadd (new_c, c, g);					// c + g
					cvsub (new_g, c, g);					// c - g
					cvadd (new_d, d, h);					// d + h
					cvsub (new_h, d, h);					// d - h

					// Do the distance = 2 round, mul by inv_fft_size as efficiently as possible
#if !VFMA
					cvrmul (new_c, inv_fft_size, new_c);
					cvrmul (new_g, inv_fft_size, new_g);
					cvadd  (a, new_a, new_c);				// a + c
					cvsub  (c, new_a, new_c);				// a - c
					cvisub (e, new_e, new_g);				// e - ig
					cviadd (g, new_e, new_g);				// e + ig
#else
					cvaddfm  (a, new_a, inv_fft_size, new_c);		// a + c
					cvsubfm  (c, new_a, inv_fft_size, new_c);		// a - c
					cvisubfm (e, new_e, inv_fft_size, new_g);		// e - ig
					cviaddfm (g, new_e, inv_fft_size, new_g);		// e + ig
#endif
					cvadd  (b, new_b, new_d);				// b + d
					cvsub  (d, new_b, new_d);				// b - d
					cvisub (f, new_f, new_h);				// f - ih
					cviadd (h, new_f, new_h);				// f + ih

					// Do the distance = 1 round, mul by inv_fft_size as efficiently as possible
#if !VFMA
					cvrmul (b, inv_fft_size, b);
					cvrmul (d, inv_fft_size, d);
					cvadd  (new_a, a, b);					// a + b
					cvsub  (new_b, a, b);					// a - b
					cvisub (new_c, c, d);					// c - id
					cviadd (new_d, c, d);					// c + id
#else
					cvaddfm  (new_a, a, inv_fft_size, b);			// a + b
					cvsubfm  (new_b, a, inv_fft_size, b);			// a - b
					cvisubfm (new_c, c, inv_fft_size, d);			// c - id
					cviaddfm (new_d, c, inv_fft_size, d);			// c + id
#endif
					twid_f.real = addpd (f.real, f.imag);			// twiddled f / .707
					twid_f.imag = subpd (f.imag, f.real);
					twid_h.real = subpd (h.real, h.imag);			// -twiddled h / .707
					twid_h.imag = addpd (h.real, h.imag);
#if !VFMA
					cvrmul (twid_f, inv_fft_size_707, twid_f);		// twiddled f
					cvrmul (twid_h, inv_fft_size_707, twid_h);		// -twiddled h
					cvadd (new_e, e, twid_f);				// e + (sqrthalf-sqrthalf*i)f
					cvsub (new_f, e, twid_f);				// e - (sqrthalf-sqrthalf*i)f
					cvsub (new_g, g, twid_h);				// g + (-sqrthalf-sqrthalf*i)h
					cvadd (new_h, g, twid_h);				// g - (-sqrthalf-sqrthalf*i)h
#else
					cvaddfm (new_e, e, inv_fft_size_707, twid_f);		// e + (sqrthalf-sqrthalf*i)f
					cvsubfm (new_f, e, inv_fft_size_707, twid_f);		// e - (sqrthalf-sqrthalf*i)f
					cvsubfm (new_g, g, inv_fft_size_707, twid_h);		// g + (-sqrthalf-sqrthalf*i)h
					cvaddfm (new_h, g, inv_fft_size_707, twid_h);		// g - (-sqrthalf-sqrthalf*i)h
#endif

					dest[group + group_member] = new_a;			// Do not bit reverse outputs
					dest[group + group_member + 4*stride] = new_b;
					dest[group + group_member + 2*stride] = new_c;
					dest[group + group_member + 6*stride] = new_d;
					dest[group + group_member + 1*stride] = new_e;
					dest[group + group_member + 5*stride] = new_f;
					dest[group + group_member + 3*stride] = new_g;
					dest[group + group_member + 7*stride] = new_h;
				}
			}
			stride *= 8;			// Next stride amount
			inv_fft_applied = TRUE;
		}

		// Radix-4 FFT, one fewer complex multiply than radix-2, does more work while data is in registers
		for ( ; (size = pass_size / stride) % 4 == 0; stride *= 4, src = dest, instride1 = 1) {
			twiddle45_stride /= 4;
			if (pass == 3 && stride * 4 == pass_size) dest = outvec + blki, outstride1 = pass2_size;
			for (unsigned int group = 0; group < pass_size; group += stride * 4) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle45_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];
					CVDT d = src[instride1*(group + group_member + 3*stride)];

					// Special case first group member, no twiddles necessary
					CVDT	new_a, new_b, new_c, new_d;
					if (twiddle_idx == 0) {
						if (inv_fft_applied) {
							cvadd (new_a, a, c);		// a + c
							cvsub (new_c, a, c);		// a - c
							cvadd (new_b, b, d);		// b + d
							cvsub (new_d, b, d);		// b - d

							cvadd  (a, new_a, new_b);	// a + b
							cvsub  (b, new_a, new_b);	// a - b
							cvisub (c, new_c, new_d);	// c - id
							cviadd (d, new_c, new_d);	// c + id
						} else {
							cvrmul (a, inv_fft_size, a);
#if !VFMA
							cvrmul (b, inv_fft_size, b);
							cvrmul (c, inv_fft_size, c);
							cvrmul (d, inv_fft_size, d);

							cvadd (new_a, a, c);		// a + c
							cvsub (new_c, a, c);		// a - c
							cvadd (new_b, b, d);		// b + d
							cvsub (new_d, b, d);		// b - d

							cvadd  (a, new_a, new_b);	// a + b
							cvsub  (b, new_a, new_b);	// a - b
							cvisub (c, new_c, new_d);	// c - id
							cviadd (d, new_c, new_d);	// c + id
#else
							cvaddfm (new_a, a, inv_fft_size, c);		// a + c
							cvsubfm (new_c, a, inv_fft_size, c);		// a - c
							cvadd   (new_b, b, d);				// b + d
							cvsub   (new_d, b, d);				// b - d

							cvaddfm  (a, new_a, inv_fft_size, new_b);	// a + b
							cvsubfm  (b, new_a, inv_fft_size, new_b);	// a - b
							cvisubfm (c, new_c, inv_fft_size, new_d);	// c - id
							cviaddfm (d, new_c, inv_fft_size, new_d);	// c + id
#endif
						}
					}
					// Apply the CONJUGATE of the twiddle
					else {
						CVDT	twiddled_b, twiddled_c, twiddled_d;
						TVDT twiddle1; twidload (twiddle1, &pmdata->twiddles2[4*twiddle_idx]);
						TVDT twiddle2; twidload (twiddle2, &pmdata->twiddles2[4*twiddle_idx+2]);
#if !VFMA
						twidconjmul (twiddled_b, b, twiddle1);
						twidconjmul (twiddled_c, c, twiddle2);
						twidmul     (twiddled_d, d, twiddle1);
						cvadd (new_a, a, twiddled_c);				// a + c
						cvsub (new_c, a, twiddled_c);				// a - c
						cvadd (new_b, twiddled_b, twiddled_d);			// b + d
						cvsub (new_d, twiddled_b, twiddled_d);			// b - d

						cvadd  (a, new_a, new_b);				// a + b
						cvsub  (b, new_a, new_b);				// a - b
						cvisub (c, new_c, new_d);				// c - id
						cviadd (d, new_c, new_d);				// c + id
#else
						twidconjmuldelay (twiddled_b, b, twiddle1);
						twidconjmuldelay (twiddled_c, c, twiddle2);
						twidmuldelay     (twiddled_d, d, twiddle1);
						cvaddfm (new_a, a, twiddle2.sin, twiddled_c);		// a + c
						cvsubfm (new_c, a, twiddle2.sin, twiddled_c);		// a - c
						cvadd (new_b, twiddled_b, twiddled_d);			// b + d
						cvsub (new_d, twiddled_b, twiddled_d);			// b - d

						cvaddfm  (a, new_a, twiddle1.sin, new_b);		// a + b
						cvsubfm  (b, new_a, twiddle1.sin, new_b);		// a - b
						cvisubfm (c, new_c, twiddle1.sin, new_d);		// c - id
						cviaddfm (d, new_c, twiddle1.sin, new_d);		// c + id
#endif
					}

					// Save results
					dest[outstride1*(group + group_member)]            = a;		// Do not bit reverse outputs
					dest[outstride1*(group + group_member + 2*stride)] = b;
					dest[outstride1*(group + group_member + 1*stride)] = c;
					dest[outstride1*(group + group_member + 3*stride)] = d;
				}
			}
			inv_fft_applied = TRUE;
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
		for ( ; (size = pass_size / stride) % 5 == 0; stride *= 5, src = dest, instride1 = 1) {
			twiddle45_stride /= 5;
			if (pass == 3 && stride * 5 == pass_size) dest = outvec + blki, outstride1 = pass2_size;
			VDT _309 = broadcastsd (0.30901699437494742410229341718282);
			VDT _809 = broadcastsd (0.80901699437494742410229341718282);
			VDT _951 = broadcastsd (0.95105651629515357211643933337938);
			VDT _588_951 = broadcastsd (0.61803398874989484820458683436564);
			for (unsigned int group = 0; group < pass_size; group += stride * 5) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle45_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];
					CVDT d = src[instride1*(group + group_member + 3*stride)];
					CVDT e = src[instride1*(group + group_member + 4*stride)];

					// Special case first group member, no twiddles necessary
					CVDT	b_plus_e, b_minus_e, c_plus_d, c_minus_d;
					if (twiddle_idx == 0) {
						cvadd (b_plus_e, b, e);				// b + e
						cvsub (b_minus_e, b, e);			// b - e
						cvadd (c_plus_d, c, d);				// c + d
						cvsub (c_minus_d, c, d);			// c - d
					}
					// Apply the CONJUGATE of the twiddle
					else {
						CVDT twiddled_b, twiddled_c, twiddled_d, twiddled_e;
						TVDT twiddle1; twidload (twiddle1, &pmdata->twiddles2[4*twiddle_idx]);
						TVDT twiddle2; twidload (twiddle2, &pmdata->twiddles2[4*twiddle_idx+2]);
						twidconjmul (twiddled_b, b, twiddle1);
						twidconjmul (twiddled_c, c, twiddle2);
#if !VFMA
						twidmul (twiddled_d, d, twiddle2);
						twidmul (twiddled_e, e, twiddle1);
						cvadd (b_plus_e, twiddled_b, twiddled_e);			// b + e
						cvsub (b_minus_e, twiddled_b, twiddled_e);			// b - e
						cvadd (c_plus_d, twiddled_c, twiddled_d);			// c + d
						cvsub (c_minus_d, twiddled_c, twiddled_d);			// c - d
#else
						twidmuldelay (twiddled_d, d, twiddle2);
						twidmuldelay (twiddled_e, e, twiddle1);
						cvaddfm (b_plus_e, twiddled_b, twiddle1.sin, twiddled_e);	// b + e
						cvsubfm (b_minus_e, twiddled_b, twiddle1.sin, twiddled_e);	// b - e
						cvaddfm (c_plus_d, twiddled_c, twiddle2.sin, twiddled_d);	// c + d
						cvsubfm (c_minus_d, twiddled_c, twiddle2.sin, twiddled_d);	// c - d
#endif
					}

					CVDT tmp25a, tmp25b, tmp34a, tmp34b;
					cvaddfm (tmp25a, a, _309, b_plus_e); cvsubfm (tmp25a, tmp25a, _809, c_plus_d);	// tmp25a = a +.309*(b+e) -.809*(c+d)
					cvsubfm (tmp34a, a, _809, b_plus_e); cvaddfm (tmp34a, tmp34a, _309, c_plus_d);	// tmp34a = a -.809*(b+e) +.309*(c+d)
					cvaddfm (tmp25b, b_minus_e, _588_951, c_minus_d);				// tmp25b = delayed .951*(b-e) + 0.588*(c-d)
					cvsubfm (tmp34b, c_minus_d, _588_951, b_minus_e);				// tmp34b = delayed -.588*(b-e) + 0.951*(c-d)
#if !VFMA
					cvrmul (tmp25b, _951, tmp25b);				// tmp25b *= .951
					cvrmul (tmp34b, _951, tmp34b);				// tmp34b *= .951
					cvisub (b, tmp25a, tmp25b);				// b = tmp25a - i*tmp25b
					cviadd (e, tmp25a, tmp25b);				// e = tmp25a + i*tmp25b
					cviadd (c, tmp34a, tmp34b);				// c = tmp34a - i*tmp34b
					cvisub (d, tmp34a, tmp34b);				// d = tmp34a + i*tmp34b
#else
					cvisubfm (b, tmp25a, _951, tmp25b);			// b = tmp25a - i*.951*tmp25b
					cviaddfm (e, tmp25a, _951, tmp25b);			// e = tmp25a + i*.951*tmp25b
					cviaddfm (c, tmp34a, _951, tmp34b);			// c = tmp34a - i*.951*tmp34b
					cvisubfm (d, tmp34a, _951, tmp34b);			// d = tmp34a + i*.951*tmp34b
#endif
					cvadd (a, a, b_plus_e);					// a + b + c + d + e
					cvadd (a, a, c_plus_d);

					// Save results
					dest[outstride1*(group + group_member)]            = a;
					dest[outstride1*(group + group_member + stride)]   = b;
					dest[outstride1*(group + group_member + 2*stride)] = c;
					dest[outstride1*(group + group_member + 3*stride)] = d;
					dest[outstride1*(group + group_member + 4*stride)] = e;
				}
			}
		}

		// Radix-3 FFT
		// A 3-complex inverse FFT is:
		// Res1: (R1+R2+R3) + (I1+I2+I3)i
		// Res2:  (R1-.5R2-.5R3+.866I2-.866I3) + (I1-.5I2-.5I3-.866R2+.866R3)i
		// Res3:  (R1-.5R2-.5R3-.866I2+.866I3) + (I1-.5I2-.5I3+.866R2-.866R3)i
		for ( ; (size = pass_size / stride) % 3 == 0; stride *= 3, src = dest, instride1 = 1) {
			twiddle3_stride /= 3;
			if (pass == 3 && stride * 3 == pass_size) dest = outvec + blki, outstride1 = pass2_size;
			VDT half = broadcastsd (0.5);
			VDT _866 = broadcastsd (0.86602540378443864676372317075294);
			for (unsigned int group = 0; group < pass_size; group += stride * 3) {
				for (unsigned int group_member = 0; group_member < stride; group_member++) {
					unsigned int twiddle_idx = ((group_member * pass_stride) + blki) * twiddle3_stride;
					CVDT a = src[instride1*(group + group_member)];
					CVDT b = src[instride1*(group + group_member + stride)];
					CVDT c = src[instride1*(group + group_member + 2*stride)];

					// Special case first group member, no twiddles necessary
					CVDT	b_plus_c, b_minus_c;
					if (twiddle_idx == 0) {
						cvadd (b_plus_c, b, c);			// b + c
						cvsub (b_minus_c, b, c);		// b - c
					}
					// Apply the CONJUGATE of the twiddle
					else {
						CVDT twiddled_b, twiddled_c;
						TVDT twiddle1; twidload (twiddle1, &pmdata->twiddles1[2*twiddle_idx]);
						twidconjmul  (twiddled_b, b, twiddle1);
#if !VFMA
						twidmul (twiddled_c, c, twiddle1);
						cvadd (b_plus_c, twiddled_b, twiddled_c);			// b + c
						cvsub (b_minus_c, twiddled_b, twiddled_c);			// b - c
#else
						twidmuldelay (twiddled_c, c, twiddle1);
						cvaddfm (b_plus_c, twiddled_b, twiddle1.sin, twiddled_c);	// b + c
						cvsubfm (b_minus_c, twiddled_b, twiddle1.sin, twiddled_c);	// b - c
#endif
					}

					CVDT tmp23; cvsubfm (tmp23, a, half, b_plus_c);			// tmp23 = a - 0.5 * (b + c)
#if !VFMA
					cvrmul (b_minus_c, _866, b_minus_c);
					cvisub (b, tmp23, b_minus_c);					// b = tmp23 - i * .866 * (b - c)
					cviadd (c, tmp23, b_minus_c);					// c = tmp23 + i * .866 * (b - c)
#else
					cvisubfm (b, tmp23, _866, b_minus_c);				// b = tmp23 - i * .866 * (b - c)
					cviaddfm (c, tmp23, _866, b_minus_c);				// c = tmp23 + i * .866 * (b - c)
#endif
					cvadd (a, a, b_plus_c);						// a + b + c

					dest[outstride1*(group + group_member)] = a;
					dest[outstride1*(group + group_member + stride)] = b;
					dest[outstride1*(group + group_member + 2*stride)] = c;
				}
			}
		}
	    }
	    }
	}
}

// outvec += addin vector shifted the proper amount to create a monic polymult result
__inline void monic_line_add_hi (CVDT *fft1, CVDT *addinvec, int addinvec_size, CVDT *outvec, int outvec_size, int circular_size, bool leading_one, bool trailing_one)
{
	// Add (or multiply and add) words starting from the highest word
	int outvec_index = (outvec_size-1) % circular_size;
	for (addinvec += addinvec_size-1; addinvec_size; addinvec--, addinvec_size--, outvec_index--) {
		if (outvec_index < 0) outvec_index = circular_size - 1;
		if (leading_one) { leading_one = FALSE; continue; }
		if (trailing_one && addinvec_size == 1) continue;
		if (fft1 == NULL) {
			cvadd (outvec[outvec_index], outvec[outvec_index], addinvec[0]);
		} else {
			CVDT	tmp;
			cvmul (tmp, *fft1, addinvec[0]);
			cvadd (outvec[outvec_index], outvec[outvec_index], tmp);
		}
	}
}

// outvec += addin vector shifted to the low part of outvec
__inline void monic_line_add_lo (CVDT *fft1, CVDT *addinvec, int addinvec_size, CVDT *outvec, int outvec_size, int circular_size, bool leading_one, bool trailing_one)
{
	// Add (or multiply and add) words starting from the lowest word
	for (int outvec_index = 0; addinvec_size; addinvec++, addinvec_size--, outvec_index++) {
		if (outvec_index == circular_size) outvec_index = 0;
		if (trailing_one) { trailing_one = FALSE; continue; }
		if (leading_one && addinvec_size == 1) continue;
		if (fft1 == NULL) {
			cvadd (outvec[outvec_index], outvec[outvec_index], addinvec[0]);
		} else {
			CVDT	tmp;
			cvmul (tmp, *fft1, addinvec[0]);
			cvadd (outvec[outvec_index], outvec[outvec_index], tmp);
		}
	}
}

// Make adjustments to outvec for any monic input vectors
void monic_line_adjustments (
	gwhandle *gwdata,
	CVDT	*addinvec,		// Polynomial to add to the output vector
	int	addinvec_size,		// Size of the polynomial to add to the output vector
	int	addinvec_options,	// Options associated with the addin vector
	bool	addinvec_monics_stripped, // TRUE if monic ones have been stripped from addin vector
	CVDT	*outvec,		// Output polynomial
	int	outvec_size,		// Size of output polynomial
	int	circular_size,		// Emulating output modulo (X^circular_size - 1)
//GW	     need LSWs_skipped? gw_outvec_size?
	int	line,			// gwnum "line" number 
	int	options)		// Options for the monic vector that had its ones stripped forcing the need to add in the addin_vector
{
	CVDT	fft1, *fft1p;

	// Set flags that will be useful in adding in the addin vector
	bool leading_one = (addinvec_options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC));
	bool trailing_one = (leading_one && (addinvec_options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)));

	// Adjust addin vector pointer so that it looks like monic ones have not been stripped
	if (leading_one && addinvec_monics_stripped) addinvec_size++;			// Increase length to "fake" a leading one
	if (trailing_one && addinvec_monics_stripped) addinvec--, addinvec_size++;	// Change vectaor start and increase length to "fake" a trailing one
   
	// If k != 1, we need to multiply the monic add-ins by FFT(1).
	if (gwdata->k != 1.0) {
		read_line (gwdata, &gwdata->GW_FFT1, 1, line, &fft1, 0, FALSE);
		fft1p = &fft1;
	} else {
		fft1p = NULL;
	}

	// Adjust for an implied high coefficient of one -- add in the poly to the high words
	monic_line_add_hi (fft1p, addinvec, addinvec_size, outvec, outvec_size, circular_size, leading_one, trailing_one);

	// Adjust for an implied low coefficient of one -- add in the poly to the low words
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))
		monic_line_add_lo (fft1p, addinvec, addinvec_size, outvec, outvec_size, circular_size, leading_one, trailing_one);
}

// Read a line that was preprocessed by polymult_preprocess
#ifndef FMA
void read_preprocess_line (gwhandle *gwdata, gwnum *data, int size, int index, CVDT *vec, int options, bool post_process_monics)
{
	preprocessed_poly_header *hdr;  // Header of the preprocessed poly under construction
	CVDT	*first_element;		// First element in the preprocessed poly array

	// Typecast to access the preprocessed poly header
	hdr = (preprocessed_poly_header *) ((char *) data - sizeof (gwarray_header));
	first_element = (CVDT *) round_up_to_multiple_of ((intptr_t) hdr + sizeof (preprocessed_poly_header), 64);

	// If preprocessed data is FFTed, then return the data "as is" optionally decompressing.  Ignore the passed in size, use the poly FFT size
	// that was computed by polymult_line_preprocess.
	if (hdr->options & POLYMULT_PRE_FFT) {
		if (!(hdr->options & POLYMULT_PRE_COMPRESS)) memcpy (vec, (char *) first_element + (intptr_t) index * (intptr_t) hdr->padded_element_size, hdr->element_size);
		else decompress_vec ((char *) first_element + (intptr_t) index * (intptr_t) hdr->padded_element_size, (char *) vec, hdr->fft_size * sizeof (CVDT));
		return;
	}

	// If not post-processing monics and the input is a monic RLP, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP))) {
		if (gwdata->k == 1.0) {
			vec[0].real = broadcastsd (1.0);
			vec[0].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, index, vec, 0, FALSE);
		}
		vec++;
	}

	// If input is an RLP then read the data into the high half of vec.  We will duplicate the data as the last step.
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) vec += (size-1);

	// If we've delayed splitting the two real values in !ALL_COMPLEX_FFT, then read the first array entry for the first two indexes
	int	index_to_read = index;
	if (!gwdata->ALL_COMPLEX_FFT && !gwdata->ZERO_PADDED_FFT) {
		index_to_read--;
		if (index_to_read < 0) index_to_read = 0;
	}

	// Test for valid index
	ASSERTG (index_to_read < hdr->num_lines);

	// Copy the preprocessed entry
	ASSERTG ((hdr->options & POLYMULT_PRE_COMPRESS) || hdr->element_size == size * sizeof (CVDT));
	if (!(hdr->options & POLYMULT_PRE_COMPRESS)) memcpy (vec, (char *) first_element + (intptr_t) index_to_read * (intptr_t) hdr->padded_element_size, hdr->element_size);
	else decompress_vec ((char *) first_element + (intptr_t) index_to_read * (intptr_t) hdr->padded_element_size, (char *) vec, size * sizeof (CVDT));

	// If we've delayed splitting the two real values in !ALL_COMPLEX_FFT turn the two reals into two complex values here
	if (!gwdata->ALL_COMPLEX_FFT && !gwdata->ZERO_PADDED_FFT && index_to_read == 0) {
		for (int j = 0; j < size; j++) {
			union { VDT a; double b[VLEN / sizeof(double)]; } r, i;
			r.a = vec[j].real;
			i.a = vec[j].imag;
			if (index == 0) {
				memset (&r, 0, sizeof (r));
				r.b[0] = i.b[0];
				memset (&i, 0, sizeof (i));
			} else {
				i.b[0] = 0.0;
			}
			vec[j].real = r.a;
			vec[j].imag = i.a;
		}
	}

	// Duplicate RLP data in reverse order
	if (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)) {
		for (int j = 1; j < size; j++) vec[-j] = vec[j];
	}

	// If not post-processing monics and the input is monic, then output a 1+0i
	if (!post_process_monics && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC))) {
		if (gwdata->k == 1.0) {
			vec[size].real = broadcastsd (1.0);
			vec[size].imag = broadcastsd (0.0);
		} else {
			read_line (gwdata, &gwdata->GW_FFT1, 1, index, &vec[size], 0, FALSE);
		}
	}
}
#endif

/* Preprocess a poly that will be used in multiple polymult calls.  Preprocessing can reduce memory consumption or reduce CPU time. */
/* Returns a massaged poly.  Caller should then free the unmassaged poly.  The massaged poly cannot be used in any gwnum calls, it can only be used */
/* in future polymult calls with poly sizes and options that match those passed to this routine. */
gwarray polymult_line_preprocess (	// Returns a plug-in replacement for the input poly
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*gw_invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	int	invec2_size,		// Size of the second input polynomial
	int	outvec_size,		// Size of the output polynomial
	int	options)
{
	int	adjusted_invec1_size, adjusted_invec2_size, adjusted_outvec_size, fft_size, num_lines, element_size, padded_element_size, max_element_size;
	int	true_invec1_size, true_invec2_size, true_outvec_size;
	bool	must_fft, monic_ones_included;
	CVDT	*invec1, *scratch;
	preprocessed_poly_header *hdr;  // Header of the preprocessed poly under construction
	CVDT	*first_element;		// First element in the preprocessed poly array

	// Adjust sizes due to RLPs.  Calculate sizes assuming no monic ones in data and sizes assuming monic ones are in the data.
	adjusted_invec1_size = (options & POLYMULT_INVEC1_RLP) ? 2 * invec1_size - 1 : invec1_size;
	adjusted_invec2_size = (options & POLYMULT_INVEC2_RLP) ? 2 * invec2_size - 1 : invec2_size;
	adjusted_outvec_size = adjusted_invec1_size + adjusted_invec2_size - 1;
	true_invec1_size = adjusted_invec1_size + ((options & POLYMULT_INVEC1_MONIC) ? ((options & POLYMULT_INVEC1_RLP) ? 2 : 1) : 0);
	true_invec2_size = adjusted_invec2_size + ((options & POLYMULT_INVEC2_MONIC) ? ((options & POLYMULT_INVEC2_RLP) ? 2 : 1) : 0);
	true_outvec_size = true_invec1_size + true_invec2_size - 1;

	// Sometimes we must (or it is beneficial) handle monic input polynomials as part of invec1/invec2 as opposed to a post-processing call to
	// monic_line_adjustments.  Find those cases and adjust vector sizes.

	// If both inputs are monic, one of the inputs must have their ones stripped.  It has to be invec1 since once we FFT the poly invec1 data is gone.
	// Stripping ones from invec1 means post proceessing will add invec2.
	if ((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC2_MONIC)) {
		monic_ones_included = FALSE;
	}
	// Next case is a circular multiply where the poly FFT can do the circular operation for free
	else if (options & POLYMULT_CIRCULAR && outvec_size == polymult_fft_size (outvec_size)) {
		monic_ones_included = TRUE;
		adjusted_invec1_size = true_invec1_size;
		adjusted_invec2_size = true_invec2_size;
		adjusted_outvec_size = outvec_size;
	}
	// Next case is where including the monic ones would result in a larger FFT size
	else if (polymult_fft_size (adjusted_outvec_size) < polymult_fft_size (true_outvec_size)) {
		monic_ones_included = FALSE;
	}
	// Finally, include the monic ones
	else {
		monic_ones_included = TRUE;
		adjusted_invec1_size = true_invec1_size;
		adjusted_invec2_size = true_invec2_size;
		adjusted_outvec_size = true_outvec_size;
	}

	// Allocate memory for polymult preprocessing
	must_fft = FALSE;
//GW: The POLYMULT_FFT option is not a commandment to pre-FFT the poly.  It a request to pre-FFT if polymult would normally choose to do an FFT multiply.
//GW: Do we want a POLYMULT_FORCE_FFT option?  or have caller temporarily set BRUTE_BREAK and KARAT_BREAK?
//	must_fft = options & POLYMULT_FORCE_FFT;
	if ((!must_fft && adjusted_outvec_size < pmdata->KARAT_BREAK)) {
		invec1 = (CVDT *) aligned_malloc (adjusted_invec1_size * sizeof (CVDT), 64);
		scratch = NULL;
		options &= ~POLYMULT_PRE_FFT;
	}
	else if (!must_fft && !(options & POLYMULT_CIRCULAR) && adjusted_outvec_size < pmdata->FFT_BREAK) {
		invec1 = (CVDT *) aligned_malloc (adjusted_invec1_size * sizeof (CVDT), 64);
		scratch = NULL;
		options &= ~POLYMULT_PRE_FFT;
	}
	else {
		unsigned int pass1_size, pass2_size;
		fft_size = polymult_fft_size (adjusted_outvec_size);
		pick_pass_sizes (pmdata, fft_size, &pass1_size, &pass2_size);
		invec1 = (CVDT *) aligned_malloc (fft_size * sizeof (CVDT), 64);
		if (pass1_size > 1) scratch = (CVDT *) aligned_malloc (pass1_size * sizeof (CVDT), 64);
		else scratch = NULL;
	}

	// If not-FFTing then clear monic_ones_included so that read_line won't add 1+0i to invec1.  We'll add 1+0i later in read_preprocess_line.
	if (!(options & POLYMULT_PRE_FFT)) monic_ones_included = FALSE;

// Allocate space for the preprocessed poly under construction

	// Compute number of lines the same way polymult does
	num_lines = gwfftlen (pmdata->gwdata);
	if (!pmdata->gwdata->ALL_COMPLEX_FFT) num_lines += 2;
	if (pmdata->gwdata->ZERO_PADDED_FFT) num_lines += 14;
	num_lines = divide_rounding_up (num_lines, sizeof (CVDT) / sizeof (double));

	// Compute size of each line to be written
	if (options & POLYMULT_PRE_FFT) element_size = sizeof (CVDT) * fft_size;
	else element_size = sizeof (CVDT) * invec1_size;
//GW:  shrink element-size if/when we support sparse non-ffted invec1
	padded_element_size = round_up_to_multiple_of (element_size, 64);

	// Allocate and init the preprocessed output
	hdr = (preprocessed_poly_header *) malloc (sizeof (preprocessed_poly_header) + 64 + (intptr_t) num_lines * (intptr_t) padded_element_size);
	if (hdr == NULL) return (NULL);
	memset (hdr, 0, sizeof (preprocessed_poly_header));
	hdr->num_lines = num_lines;
	hdr->element_size = element_size;
	hdr->padded_element_size = padded_element_size;
	hdr->options = options;
	hdr->monic_ones_included = monic_ones_included;
	if (options & POLYMULT_PRE_FFT) hdr->fft_size = fft_size;
	first_element = (CVDT *) round_up_to_multiple_of ((intptr_t) hdr + sizeof (preprocessed_poly_header), 64);

	// When compressing we need to track the size of the largest compressed vector (for resizing the array smaller)
	max_element_size = 0;

	// Grab a "line" of complex values from each gwnum coefficient.  Work and copy that data set.  Then move onto remaining lines in the gwnums.
	char *outptr = (char *) first_element;
	for (int line = 0; line < num_lines; line++) {

		// Read a line.  If we're saving memory by delaying splitting the two reals, then combine the first two lines into one.
		if (line == 0 && !pmdata->gwdata->ALL_COMPLEX_FFT && !pmdata->gwdata->ZERO_PADDED_FFT && !(options & POLYMULT_PRE_FFT)) {
			CVDT *temp = (CVDT *) aligned_malloc (adjusted_invec1_size * sizeof (CVDT), 64);
			read_line (pmdata->gwdata, gw_invec1, invec1_size, 0, temp, options & INVEC1_OPTIONS, !monic_ones_included);
			read_line (pmdata->gwdata, gw_invec1, invec1_size, 1, invec1, options & INVEC1_OPTIONS, !monic_ones_included);
			for (int j = 0; j < adjusted_invec1_size; j++) {
				double	line0_real;
				union { VDT a; double b[VLEN / sizeof(double)]; } x;
				// Extract the real from line 0
				x.a = temp[j].real; line0_real = x.b[0];
				// Put line 0's real into invec1 from line 1
				x.a = invec1[j].imag; x.b[0] = line0_real; invec1[j].imag = x.a;
			}
			line++;
			aligned_free (temp);
		} else
			read_line (pmdata->gwdata, gw_invec1, invec1_size, line, invec1, options & INVEC1_OPTIONS, !monic_ones_included);

		// Optionally forward FFT the line
		if (options & POLYMULT_PRE_FFT) {
			memset (invec1 + adjusted_invec1_size, 0, (fft_size - adjusted_invec1_size) * sizeof (CVDT));
			fft_line (pmdata, invec1, NULL, NULL, scratch, fft_size, FORWARD_FFT_INVEC1 | FORWARD_FFT_ONLY);
		}

		// For non-FFTed RLPs, we only save the upper half of invec1.  Alternatively, we could make read_line not expand RLPs when preprocessing.
		if ((options & POLYMULT_INVEC1_RLP) && !(options & POLYMULT_PRE_FFT)) invec1 += invec1_size-1;

		// If compressing each line, compress the vector of doubles.  Track the maximum compressed vector size.
		if (options & POLYMULT_PRE_COMPRESS) {
			element_size = compress_vec ((char *) invec1, hdr->element_size);
			ASSERTG (element_size <= hdr->element_size);
			if (element_size > max_element_size) max_element_size = element_size;
		}

		// Copy line to the preprocessed poly under construction
		memcpy (outptr, invec1, element_size);
		outptr += padded_element_size;
			    
		// Undo RLP pointer hack
		if ((options & POLYMULT_INVEC1_RLP) && !(options & POLYMULT_PRE_FFT)) invec1 -= invec1_size-1;
	}

	// Free memory
	aligned_free (invec1);
	aligned_free (scratch);

	// Shrink a compressed poly to its minimum possible size using realloc.
	if (options & POLYMULT_PRE_COMPRESS) {
		for (int line = 1; line < num_lines; line++) {
			memmove ((char *) first_element + (intptr_t) line * (intptr_t) max_element_size,
				 (char *) first_element + (intptr_t) line * (intptr_t) padded_element_size, max_element_size);
		}
		hdr = (preprocessed_poly_header *) realloc (hdr, sizeof (preprocessed_poly_header) + 64 + (intptr_t) num_lines * (intptr_t) max_element_size);
		hdr->element_size = max_element_size;
		hdr->padded_element_size = max_element_size;
	}

	// Link into gwnum's list of arrays to free when gwdone is called
	hdr->linkage.next = pmdata->gwdata->array_list;		// Next in linked list of allocated arrays
	hdr->linkage.prev = &pmdata->gwdata->array_list;	// Prev in linked list of allocated arrays
	if (hdr->linkage.next != NULL) {			// Patch new next's prev pointer
		gwarray_header *next_header = (gwarray_header *) ((char *) hdr->linkage.next - sizeof (gwarray_header));
		next_header->prev = &hdr->linkage.next;
	}
	pmdata->gwdata->array_list = (gwarray) &hdr->self_ptr;	// New head of doubly linked list

	// Return the address of the self_ptr that identifies this as a preprocessed poly
	hdr->self_ptr = (gwnum *) &hdr->self_ptr;
	return ((gwnum *) &hdr->self_ptr);
}

// Multiply polynomials by selecting from various polymult algorithms
void polymult_line (
	pmhandle *pmdata)		// Handle for polymult library
{
	int	true_invec1_size, true_invec2_size, true_outvec_size, alloc_scratch_size;
	int	invec1_size, invec2_size, outvec_size, fft_size, fft_options, line, num_lines;
	bool	invec1_preprocessed, invec1_pre_ffted, invec2_preprocessed, invec2_pre_ffted;
	CVDT	*invec1_buffer, *invec2_buffer, *outvec_buffer, *tmpvec_buffer, *scratch_buffer;

	// Grab polymult arguments out of the pmdata structure.  They were copied there so all threads have access to them.
	gwnum	*gw_invec1 = pmdata->invec1;			// First input poly
	int	gw_invec1_size = pmdata->invec1_size;		// Size of the first input polynomial

	// Check for pre-processed input poly
	invec1_preprocessed = is_preprocessed_poly (gw_invec1);
	invec1_pre_ffted = invec1_preprocessed && is_preffted_poly (gw_invec1);

	// Calculate the scratch area size.  This is the one buffer size that cannot be calculated during the planning process.
	alloc_scratch_size = 0;
	for (int i = 0; i < pmdata->num_other_polys; i++) {
		unsigned int pass1_size, pass2_size;
		polymult_plan *plan = &pmdata->plan[i];
		if (plan->impl != POLYMULT_IMPL_FFT) continue;
		pick_pass_sizes (pmdata, plan->fft_size, &pass1_size, &pass2_size);
		if (pass1_size > 1 && (int) pass1_size > alloc_scratch_size) alloc_scratch_size = pass1_size;
	}

	// Allocate memory for polymult.  For the output vector allocate two extra in case there are any monic RLP inputs.
	invec1_buffer = (CVDT *) (pmdata->alloc_invec1_size ? aligned_malloc (pmdata->alloc_invec1_size * sizeof (CVDT), 64) : NULL);
	invec2_buffer = (CVDT *) (pmdata->alloc_invec2_size ? aligned_malloc (pmdata->alloc_invec2_size * sizeof (CVDT), 64) : NULL);
	tmpvec_buffer = (CVDT *) (pmdata->alloc_tmpvec_size ? aligned_malloc (pmdata->alloc_tmpvec_size * sizeof (CVDT), 64) : NULL);
	outvec_buffer = (CVDT *) aligned_malloc (pmdata->alloc_outvec_size * sizeof (CVDT), 64);
	scratch_buffer = (CVDT *) (alloc_scratch_size ? aligned_malloc (alloc_scratch_size * sizeof (CVDT), 64) : NULL);

//GW:   alternative if poly FFT size is large, reread invec1/2 in order to post process addin to save memory! Would this help P-1??? NO.
//GW:		code would be similar to fma processing

	// Grab a "line of complex values" from each gwnum coefficient.  Work that data set.  Then move onto remaining lines in the gwnums.
	num_lines = gwfftlen (pmdata->gwdata);
	if (!pmdata->gwdata->ALL_COMPLEX_FFT) num_lines += 2;
	if (pmdata->gwdata->ZERO_PADDED_FFT) num_lines += 14;
	num_lines = divide_rounding_up (num_lines, sizeof (CVDT) / sizeof (double));
	while ((line = atomic_fetch_incr (pmdata->next_line)) < num_lines) {
	    CVDT *invec1, *invec2, *outvec, *true_outvec, *tmpvec, *scratch;
	    tmpvec = tmpvec_buffer;
	    scratch = scratch_buffer;

	    // Loop over all the other polys to multiply with the first poly
	    int prev_fft_size = 0;
	    for (int i = 0; i < pmdata->num_other_polys; i++) {
		gwnum	*gw_invec2 = pmdata->other_polys[i].invec2;		// Second input poly
		int	gw_invec2_size = pmdata->other_polys[i].invec2_size;	// Size of the second input polynomial
		gwnum	*gw_outvec = pmdata->other_polys[i].outvec;		// Output poly of size invec1_size + invec2_size - 1
		int	gw_outvec_size = pmdata->other_polys[i].outvec_size;	// Size of the output polynomial.  Should be invec1_size + invec2_size (less one if not monic).
		gwnum	*gw_fmavec = pmdata->other_polys[i].fmavec;		// Addin for FMA operations
		int	options = pmdata->options | pmdata->other_polys[i].options;
		polymult_plan *plan = &pmdata->plan[i];

// Reload saved fftsize, post_monics, adjusted_sizes and anything else we need from the plan computed in polymult_several
		bool	strip_monic_from_invec1;	// True if ones are stripped from monic invec1 requiring invec2 to be added in during post processing
		bool	post_monic_addin_invec2;	// True if ones are stripped from monic invec1 requiring invec2 to be added in during post processing
		bool	strip_monic_from_invec2;	// True if ones are stripped from monic invec2 requiring invec1 to be added in during post processing
		bool	post_monic_addin_invec1;	// True if ones are stripped from monic invec2 requiring invec1 to be added in during post processing
		int	circular_size, LSWs_skipped;

		strip_monic_from_invec1 = plan->strip_monic_from_invec1;
		strip_monic_from_invec2 = plan->strip_monic_from_invec2;
		post_monic_addin_invec1 = plan->post_monic_addin_invec1;
		post_monic_addin_invec2 = plan->post_monic_addin_invec2;
		fft_size = plan->fft_size;
		circular_size = plan->circular_size;
		LSWs_skipped = plan->LSWs_skipped;

		// When monic ones are stripped from inputs, the brute/Karatsuba/polyFFT code creates an output that is smaller than the true outvec result size.
		// The unprefixed variables refer to the monic stripped input and output sizes, the true_XXX variables refer to the full input and output sizes.
		invec1_size = plan->adjusted_invec1_size;
		invec2_size = plan->adjusted_invec2_size;
		outvec_size = plan->adjusted_outvec_size;
		true_invec1_size = plan->true_invec1_size;
		true_invec2_size = plan->true_invec2_size;
		true_outvec_size = plan->true_outvec_size;

		// Pad/shift the outvec pointer so that multiplication results land in the correct location of the true outvec.
		true_outvec = outvec_buffer;
		outvec = outvec_buffer + plan->adjusted_shift;

		// Pad output vector on the right with a zero or two to reach the true_outvec_size
		for (int i = 0; i < plan->adjusted_shift; i++) memset (&true_outvec[i], 0, sizeof (CVDT));
		// Pad output vector on the right with a zero or two to reach the true_outvec_size
		for (int i = 1; i <= plan->adjusted_pad; i++) memset (&true_outvec[true_outvec_size-i], 0, sizeof (CVDT));

		// Check for pre-processed other poly and FMA addin vector
		invec2_preprocessed = is_preprocessed_poly (gw_invec2);
		invec2_pre_ffted = invec2_preprocessed && is_preffted_poly (gw_invec2);

		// Read invec1 on processing of first other poly.  For FFT implementation we can usually read invec1 into a buffer shared
		// with tmpvec (the buffer that will hold the FFT of invec1).
		if (i == 0) {
			invec1 = (invec1_buffer == NULL ? tmpvec : invec1_buffer);
			if (!invec1_preprocessed) read_line (pmdata->gwdata, gw_invec1, gw_invec1_size, line, invec1, options & INVEC1_OPTIONS, strip_monic_from_invec1);
			else read_preprocess_line (pmdata->gwdata, gw_invec1, gw_invec1_size, line, invec1, options & INVEC1_OPTIONS, strip_monic_from_invec1);
		}
		// Plan should not change strip_monic_from_invec1 so that we can read invec1 just once (very important since outvecs can output to invec1)
		ASSERTG (i == 0 || strip_monic_from_invec1 == pmdata->plan[i-1].strip_monic_from_invec1);

		// Read invec2.  For FFT implementation we can usually read invec2 into a buffer shared with outvec.
		invec2 = (plan->impl == POLYMULT_IMPL_FFT && !post_monic_addin_invec2) ? outvec : invec2_buffer;
		if (!invec2_preprocessed) read_line (pmdata->gwdata, gw_invec2, gw_invec2_size, line, invec2, options & INVEC2_OPTIONS, strip_monic_from_invec2);
		else read_preprocess_line (pmdata->gwdata, gw_invec2, gw_invec2_size, line, invec2, options & INVEC2_OPTIONS, strip_monic_from_invec2);

		// Now execute a brute force implementation.  Brute_line does implement circular_size, that will not need to be emulated later.
		// Because brute_line implements circular_size, it must get true_outvec buffer.
		// Brute_line does peer into gw_outvec to see if output coefficients need to be computed -- a minor optimization.
		if (plan->impl == POLYMULT_IMPL_BRUTE) {
			brute_line (invec1, invec1_size, invec2, invec2_size, true_outvec, circular_size, gw_outvec, gw_outvec_size, plan->adjusted_shift, LSWs_skipped);
		}

		// Now execute a Karatsuba implementation
		else if (plan->impl == POLYMULT_IMPL_KARATSUBA) {
			karatsuba_line (pmdata, invec1, invec1_size, invec2, invec2_size, tmpvec, outvec);
		}

		// Now execute a poly FFT implementation
		else {
			// FFT invec1 only on first other_poly OR changed fft sizes
			if (invec1_pre_ffted || (i > 0 && fft_size == prev_fft_size)) fft_options = 0;
			else fft_options = FORWARD_FFT_INVEC1;
			// FFT invec2 unless it was done during pre-processing
			if (!invec2_pre_ffted) fft_options |= FORWARD_FFT_INVEC2;
			// If we're FFTing invec1, copy and/or zero pad as necessary
			if (fft_options & FORWARD_FFT_INVEC1) {
				if (invec1 != tmpvec) memcpy (tmpvec, invec1, invec1_size * sizeof (CVDT));
				memset (tmpvec + invec1_size, 0, (fft_size - invec1_size) * sizeof (CVDT));
			}
			// If we're FFTing invec2, copy and/or zero pad as necessary
			if (fft_options & FORWARD_FFT_INVEC2) {
				if (invec2 != outvec) memcpy (outvec, invec2, invec2_size * sizeof (CVDT));
				memset (outvec + invec2_size, 0, (fft_size - invec2_size) * sizeof (CVDT));
			}
			// Now do the necessary forward FFTs, pointwise multiply, and inverse FFT.
			fft_line (pmdata, tmpvec, outvec, outvec, scratch, fft_size, fft_options);
			// Remember poly FFT size for next invec2
			prev_fft_size = fft_size;
		}

		// Emulate circular_size when necessary.  We only do this for coefficients that might actual be returned to the caller.
		if (plan->emulate_circular) {
			int	i, adjusted_i;		// i is index into true_outvec, adjusted_i is i % circular_size
			for (i = circular_size + LSWs_skipped, adjusted_i = LSWs_skipped; i < true_outvec_size; ) {
				if (adjusted_i >= LSWs_skipped + gw_outvec_size) {
					i = divide_rounding_up (i, circular_size) * circular_size + LSWs_skipped;
					adjusted_i = LSWs_skipped;
				} else {
					cvadd (true_outvec[adjusted_i], true_outvec[adjusted_i], true_outvec[i]);
					i++, adjusted_i++;
				}
			}
		}

		// Handle adjustments needed due to monic inputs
//GW: circular monic adjust? mulhi/mid/lo can skip chunks of addin?
		if (post_monic_addin_invec1)
			monic_line_adjustments (pmdata->gwdata, invec1, invec1_size, options & POLYMULT_INVEC1_MONIC_RLP, strip_monic_from_invec1,
						true_outvec, true_outvec_size, circular_size, line, options & POLYMULT_INVEC2_MONIC_RLP);
		if (post_monic_addin_invec2)
			monic_line_adjustments (pmdata->gwdata, invec2, invec2_size, options & POLYMULT_INVEC2_MONIC_RLP, strip_monic_from_invec2,
						true_outvec, true_outvec_size, circular_size, line, options & POLYMULT_INVEC1_MONIC_RLP);

//GWDo we emulate circular before or after monic_line_adjustments???  Do we avoid adjusting for values that are not returned??
//GWwe need to know where the top is when we're not optimizing circular			
//GW	true_outvec_size = 7	circular_size = 7			MSWs_skipped + outvec_size + LSWs_skipped = 7
//GWwe need to know where the top is when we've decided to optimizing circular
//GW	true_outvec_size = 7?   circular_size = 5			MSWs_skipped + outvec_size + LSWs_skipped = 7
//GWwe need to know where the top is when we've decided to ideal circular
//GW	true_outvec_size = 5   circular_size = 5			MSWs_skipped + outvec_size + LSWs_skipped = 5

		// Copy outvec to the output gwnums.  Also apply optional FMA operations.
		write_line (pmdata->gwdata, gw_outvec, gw_outvec_size, line, true_outvec + LSWs_skipped, gw_fmavec, options);
	    }
	}

	// Free memory
	aligned_free (invec1_buffer);
	aligned_free (invec2_buffer);
	aligned_free (outvec_buffer);
	aligned_free (tmpvec_buffer);
	aligned_free (scratch_buffer);
}

#endif


/*--------------------------------------------------------------------------
|			Polymult helper threads
+-------------------------------------------------------------------------*/

#if !defined (SSE2) && !defined (AVX) && !defined (FMA) && !defined (AVX512)

/* Entry point for helper thread to do part of the polymult */
void polymult_thread (void *arg)
{
	pmhandle *pmdata = (pmhandle *) arg;
	gwhandle cloned_gwdata;
	bool	cloned_gwdata_initialized = FALSE;

/* Generate a unique helper thread number (so that callback routine to uniquely identify helper threads) */
/* Call gwnum's optional user provided callback routine so that the caller can set the thread's priority and affinity */

	int thread_num = atomic_incr_fetch (pmdata->next_thread_num);
	if (pmdata->gwdata->thread_callback != NULL)
		(*pmdata->gwdata->thread_callback) (thread_num, 20, pmdata->gwdata->thread_callback_data);

/* Loop doing polymult or user-defined work */

	for ( ; ; ) {
	    int	num_active_helpers;

	    // Wait on the work-to-do event for more work (or termination).  There are two different ways to wait for work.
	    if (pmdata->gwdata->use_spin_wait > thread_num) atomic_spinwait (pmdata->alt_work_to_do, 1);
	    else gwevent_wait (&pmdata->work_to_do, 0);

	    // Terminate helper when polymult_done is called
	    if (pmdata->helpers_must_exit) break;

	    /* IMPORTANT:  Helpers are considered done as soon as all_work_assigned is set and num_active_helpers reaches zero.  When this happens the main */
	    /* thread is allowed to resume.  However, there may be some helpers that have not yet been given a time slice after the work_to_do event was */
	    /* signalled!  We  must catch these workers right here.  If they were to get past this section, they could read some inconsistent state that */
	    /* the main thread has changed in preparation for the next batch of work for the helper threads. */
	    if (pmdata->all_work_assigned) continue;
	    num_active_helpers = atomic_incr_fetch (pmdata->num_active_helpers);
	    if (!pmdata->all_work_assigned) {			// all_work_assigned set should be rare

		// If caller requests using fewer than the maximum number of helper threads, then some of the helpers need to do nothing
		if (num_active_helpers > pmdata->num_threads);

		// Call user-defined routine with a cloned gwdata so the user can do multi-threaded gwnum work
		else if (!pmdata->helpers_doing_polymult) {
			// Clone gwdata for user-defined callbacks to use
			if (!cloned_gwdata_initialized) {
				gwclone (&cloned_gwdata, pmdata->gwdata);
				cloned_gwdata_initialized = TRUE;
			}
			// Call the user-defined helper routine
			(*pmdata->helper_callback) (thread_num, &cloned_gwdata, pmdata->helper_callback_data);
			// After user-defined callback completes sync up all the stats and counters generated by the cloned gwdata
			gwmutex_lock (&pmdata->poly_mutex);
			if (pmdata->stats_gwdata == NULL) pmdata->stats_gwdata = &cloned_gwdata;	// Nothing to add our stats to yet
			else gwclone_merge_stats (pmdata->stats_gwdata, &cloned_gwdata);		// Add our stats to pmdata->stats_gwdata
			gwmutex_unlock (&pmdata->poly_mutex);
		}

		// Polymult optimized routines for AVX512 instructions
#ifdef X86_64
		else if (pmdata->gwdata->cpu_flags & CPU_AVX512F) polymult_line_avx512 (pmdata);
#endif
		// Polymult optimized routines for FMA3 instructions
		else if (pmdata->gwdata->cpu_flags & CPU_FMA3) polymult_line_fma (pmdata);
		// Polymult optimized routines for AVX instructions
		else if (pmdata->gwdata->cpu_flags & CPU_AVX) polymult_line_avx (pmdata);
		// Polymult optimized routines for SSE2 instructions
		else if (pmdata->gwdata->cpu_flags & CPU_SSE2) polymult_line_sse2 (pmdata);
		// Polymult default routines using no vector instructions
#ifndef X86_64
		else polymult_line_dbl (pmdata);
#endif

		// No more work to assign to helper threads
		pmdata->all_work_assigned = TRUE;		// Set flag so any helper threads that have not yet started don't try to do any work
		gwevent_reset (&pmdata->work_to_do);		// Reset event saying there is work for helper threads to do
		if (pmdata->gwdata->use_spin_wait >= 2) atomic_set (pmdata->alt_work_to_do, 0);
	    }

	    // When num active helpers reaches zero, main thread can wake up (either by waiting on all_helpers_done or via a spin wait)
	    // WARNING: It is possible for a straggler thread from the previous work_to_do event to lose its time slice just before the all_helpers_done signal call.
	    // This can cause spurious all_helpers_done signals for the current work_to_do event.  The main thread must be very careful resuming from the
	    // all_helpers_done wait by checking that num_active_helpers is in fact zero.
	    if (atomic_decr_fetch (pmdata->num_active_helpers) == 0) gwevent_signal (&pmdata->all_helpers_done);
	}

/* Call gwnum's optional user provided callback routine so that the caller can do any necessary cleanup */

	if (pmdata->gwdata->thread_callback != NULL)
		(*pmdata->gwdata->thread_callback) (thread_num, 21, pmdata->gwdata->thread_callback_data);

/* Free the cloned gwdata that user-defined callbacks used */

	if (cloned_gwdata_initialized) gwdone (&cloned_gwdata);
}

/* This routine launches the polymult helper threads.  The polymult library uses this routine to do multithreading, users can too! */
void polymult_launch_helpers (
	pmhandle *pmdata)		// Handle for polymult library
{
	// If there are no helper threads we're done
	if (pmdata->num_threads == 1) return;
	// "Neuter" gwnum to single-threaded for user-definied callback
	if (!pmdata->helpers_doing_polymult && !pmdata->helpers_must_exit) {
		pmdata->saved_gwdata_num_threads = gwget_num_threads (pmdata->gwdata);
		gwset_num_threads (pmdata->gwdata, 1);
	}
	// Launch helpers for the first time
	if (pmdata->thread_ids == NULL) {
		atomic_set (pmdata->next_thread_num, 0);
		gwmutex_init (&pmdata->poly_mutex);
		gwevent_init (&pmdata->work_to_do);
		gwevent_init (&pmdata->all_helpers_done);
		atomic_set (pmdata->alt_work_to_do, 0);			// No work for helpers to do yet
		pmdata->thread_ids = (gwthread *) malloc (pmdata->max_num_threads * sizeof (gwthread));
		for (int i = 1; i < pmdata->max_num_threads; i++) gwthread_create_waitable (&pmdata->thread_ids[i], &polymult_thread, (void *) pmdata);
	}
	// Activate helpers.  They're waiting for work to do signal.
	pmdata->stats_gwdata = NULL;
	pmdata->all_work_assigned = FALSE;			// When this is set and num_active_helpers reaches 0, it is safe to signal main thread
	gwevent_signal (&pmdata->work_to_do);			// Start all helper threads
	if (pmdata->gwdata->use_spin_wait >= 2) atomic_set (pmdata->alt_work_to_do, 1);
}

/* This routine waits for the launched polymult helper threads to finish.  This is called by the main thread when it has run out of work to do. */
void polymult_wait_on_helpers (
	pmhandle *pmdata)		// Handle for polymult library
{
	// If there are no helper threads we're done
	if (pmdata->num_threads == 1) return;
	// No more work to assign to helper threads
	pmdata->all_work_assigned = TRUE;
	gwevent_reset (&pmdata->work_to_do);
	if (pmdata->gwdata->use_spin_wait >= 2) atomic_set (pmdata->alt_work_to_do, 0);
	// Wait for helpers to end
	if (pmdata->gwdata->use_spin_wait) atomic_spinwait (pmdata->num_active_helpers, 0);
	else while (atomic_get (pmdata->num_active_helpers)) {
		gwevent_reset (&pmdata->all_helpers_done);
		if (atomic_get (pmdata->num_active_helpers)) gwevent_wait (&pmdata->all_helpers_done, 0);
	}
	// User-defined helpers require syncing up cloned gwdata stats and restoring main thread's gwnum's multithreading capabilities
	if (!pmdata->helpers_doing_polymult) {
		gwclone_merge_stats (pmdata->gwdata, pmdata->stats_gwdata);		// Add our stats to pmdata->stats_gwdata
		gwset_num_threads (pmdata->gwdata, pmdata->saved_gwdata_num_threads);
	}
}

/*--------------------------------------------------------------------------
|			Polymult main routines
+-------------------------------------------------------------------------*/

/* Preprocess a poly that will be used in multiple polymult calls.  Preprocessing can reduce memory consumption or reduce CPU time. */
/* Returns a massaged poly.  Caller should then free the unmassaged poly.  The massaged poly cannot be used in any gwnum calls, it can only be used */
/* in future polymult calls with poly sizes and options that match those passed to this routine. */
gwarray polymult_preprocess (		// Returns a plug-in replacement for the input poly
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*invec1,		// Input poly
	int	invec1_size,		// Size of the input polynomial
	int	invec2_size,		// Size of the other polynomial that will be used in a future polymult call
	int	outvec_size,		// Size of the output polynomial that will be used in a future polymult call
	int	options)		// Future polymult options preprocessing options (FFT, compress -- see polymult.h)
{
	// Cannot preprocess an already preprocessed poly
	ASSERTG (!is_preprocessed_poly(invec1));

	// FFT inputs if needed.  FFT the 7 zero padded FFT values in the header.  Even though we are not allowed to change the input vectors, we
	// can write to the 13 values values in the header reserved for use by the polymult library.
	for (int j = 0; j < invec1_size; j++) {
		if (invec1[j] == NULL) continue;
		gwfft (pmdata->gwdata, invec1[j], invec1[j]);
		if (pmdata->gwdata->ZERO_PADDED_FFT) zpad7_fft (pmdata->gwdata, invec1[j]);
	}

	// Call the instruction set specific version.
#ifdef X86_64
	if (pmdata->gwdata->cpu_flags & CPU_AVX512F) return (polymult_line_preprocess_avx512 (pmdata, invec1, invec1_size, invec2_size, outvec_size, options));
	else
#endif
	if (pmdata->gwdata->cpu_flags & CPU_FMA3) return (polymult_line_preprocess_fma (pmdata, invec1, invec1_size, invec2_size, outvec_size, options));
	else if (pmdata->gwdata->cpu_flags & CPU_AVX) return (polymult_line_preprocess_avx (pmdata, invec1, invec1_size, invec2_size, outvec_size, options));
	else if (pmdata->gwdata->cpu_flags & CPU_SSE2) return (polymult_line_preprocess_sse2 (pmdata, invec1, invec1_size, invec2_size, outvec_size, options));
#ifndef X86_64
	else return (polymult_line_preprocess_dbl (pmdata, invec1, invec1_size, invec2_size, outvec_size, options));
#endif
	return (NULL);			// Can't happen but keeps compiler happy
}

/* Multiply two polynomials producing FFTed gwnums as output.  The output gwnums need to be normalized before use in further gwnum operations. */
/* It is safe to use input gwnums in the output vector.  For a normal polynomial multiply outvec_size must be invec1_size + invec2_size - 1. */
/* For monic polynomial multiply, the leading input coefficients of 1 are omitted as is the leading 1 output coefficient -- thus */
/* outvec_size must be invec1_size + invec2_size. */
void polymult (
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	gwnum	*invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial
	gwnum	*outvec,		// Output poly
	int	outvec_size,		// Size of the output polynomial (or fft_size if POLYMULT_CIRCULAR)
	int	options)
{
	polymult_arg arg;

	// POLYMULT_CIRCULAR and (POLYMULT_MULHI or POLYMUL_MULMID or POLYMULT_MULLO) cannot both be set -- it overloads the meaning of outvec_size.
	ASSERTG (!(options & POLYMULT_CIRCULAR) || !(options & (POLYMULT_MULHI | POLYMULT_MULMID | POLYMULT_MULLO)));
	// FMA options cannot be set
	ASSERTG (!(options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)));

	arg.invec2 = invec2;
	arg.invec2_size = invec2_size;
	arg.outvec = outvec;
	arg.outvec_size = outvec_size;
	arg.fmavec = NULL;
	arg.circular_size = (options & POLYMULT_CIRCULAR) ? outvec_size : 0;
	arg.first_mulmid = 0;
	arg.options = (options & ~INVEC1_OPTIONS);

	polymult_several (pmdata, invec1, invec1_size, &arg, 1, options & INVEC1_OPTIONS);
}

   
/* Multiply two polynomials followed by a fused add or subtract */
void polymult_fma (
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	gwnum	*invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial
	gwnum	*outvec,		// Output poly
	int	outvec_size,		// Size of the output polynomial (or fft_size if POLYMULT_CIRCULAR)
	gwnum	*fmavec,		// FMA poly to add or subtract from poly multiplication result.  Same size as outvec, cannot be preprocessed.
	int	options)
{
	polymult_arg arg;

	// POLYMULT_CIRCULAR and (POLYMULT_MULHI or POLYMUL_MULMID or POLYMULT_MULLO) cannot both be set -- it overloads the meaning of outvec_size.
	ASSERTG (!(options & POLYMULT_CIRCULAR) || !(options & (POLYMULT_MULHI | POLYMULT_MULMID | POLYMULT_MULLO)));
	// FMA options must be set and fmavec cannot be NULL
	ASSERTG ((options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) && fmavec != NULL);

	arg.invec2 = invec2;
	arg.invec2_size = invec2_size;
	arg.outvec = outvec;
	arg.outvec_size = outvec_size;
	arg.fmavec = fmavec;
	arg.circular_size = (options & POLYMULT_CIRCULAR) ? outvec_size : 0;
	arg.first_mulmid = 0;
	arg.options = (options & ~INVEC1_OPTIONS);

	polymult_several (pmdata, invec1, invec1_size, &arg, 1, options & INVEC1_OPTIONS);
}

/* Multiply two polynomials.  Supports every possible polymult option.  The original polymult interface overloaded the outvec_size argument. */
void polymult2 (
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial
	gwnum	*invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial
	gwnum	*outvec,		// Output poly
	int	outvec_size,		// Size of the output polynomial.  If POLYMULT_MULHI or POLYMULT_MULLO is set this is the number of coefficients to return.
	gwnum	*fmavec,		// FMA poly to add or subtract from poly multiplication result.  Same size as outvec, cannot be preprocessed.
	int	circular_size,		// If POLYMULT_CIRCULAR set, compute poly result modulo (X^circular_size - 1)
	int	first_mulmid,		// If POLYMULT_MULMID set, this is the number of least significant coefficients that will not be returned
	int	options)
{
	polymult_arg arg;

	// POLYMULT_CIRCULAR requires circular_size to be set
	ASSERTG (((options & POLYMULT_CIRCULAR) && circular_size > 0) || (!(options & POLYMULT_CIRCULAR) && circular_size == 0));
	// POLYMUL_MULMID requires first_mulmid be set
	ASSERTG (((options & POLYMULT_MULMID) && first_mulmid >= 0) || (!(options & POLYMULT_MULMID) && first_mulmid == 0));
	// FMA options must be set and fmavec cannot be NULL
	ASSERTG (((options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) && fmavec != NULL) ||
		 (!(options & (POLYMULT_FMADD | POLYMULT_FMSUB | POLYMULT_FNMADD)) && fmavec == NULL));

	arg.invec2 = invec2;
	arg.invec2_size = invec2_size;
	arg.outvec = outvec;
	arg.outvec_size = outvec_size;
	arg.fmavec = fmavec;
	arg.circular_size = circular_size;
	arg.first_mulmid = first_mulmid;
	arg.options = (options & ~INVEC1_OPTIONS);

	polymult_several (pmdata, invec1, invec1_size, &arg, 1, options & INVEC1_OPTIONS);
}


/* Multiply one poly with several other polys.  This yields a small performance gain in that the one poly is read and FFTed just once. */
void polymult_several (		// Multiply one poly with several polys	
	pmhandle *pmdata,	// Handle for polymult library
	gwnum	*invec1,	// First input poly
	int	invec1_size,	// Size of the first input polynomial
	polymult_arg *other_polys, // Array of "other poly descriptors" to multiply with first input poly (describes each second input poly and output poly)
	int	num_other_polys,// Number of other polys to multiply with first input poly
	int	options)	// Poly #1 options.  Options not associated with poly #1 are applied to all other polys.
{
	polymult_plan *plan;		// Plan for how to implement each of the invec1 by invec2 multiplications
	int	invec1_options;		// Options that only applied to invec1
	int	global_invec2_options;	// Options that apply to all invec2s

	// Split up options.  For convenience to the programmer, invec2 options that apply to all other_polys can be specified in options argument.
	invec1_options = options & INVEC1_OPTIONS;
	global_invec2_options = options & ~INVEC1_OPTIONS;

	// Allocate array to plan each poly multiplication
	pmdata->plan = plan = (polymult_plan *) malloc (num_other_polys * sizeof (polymult_plan));

	// Check for pre-processed input poly
	bool invec1_preprocessed = is_preprocessed_poly (invec1);
	bool invec1_pre_ffted = invec1_preprocessed && is_preffted_poly (invec1);

	// For better or worse, if "other polys" would select different FFT sizes we try to make them all use the largest FFT size.  That way invec1 will
	// be forward FFTed only once.  We also need to make sure the strip_monic_from_invec1 is consistent (invec1 is only read once and that is where
	// monic ones are stripped).  To make this happen we track the largest FFT size and strip_monic_from_invec1 settings.  We way need to make two
	// passes through the planning process to make this happen.
	int	largest_fft_size = 0;	// FFT size that we can use for all the other poly plans where we are free to choose the FFT size
	bool	largest_strip_monic_from_invec1; // Whether invec1 has its monic one stripped when using largest_fft_size
	bool	must_replan;		// TRUE if replanning will be required
	bool	exists_pre_ffted_invec2 = FALSE; // TRUE if any of the invec2s are pre-FFTed.  If so, we may not be able to strip monic ones from invec1.
	bool	pass2_strip_monic_from_invec1; // In second pass, all other polys must use the same strip_monic_from_invec1 setting.

	// Determine the best plan (brute/Karatsuba/FFT, FFT size, post_process_monic, 1*1 addins, etc.).  Sometimes this takes two passes.
	for (int pass = 1; pass == 1 || must_replan; pass++) {
	    must_replan = FALSE;

	    // Track the largest buffer sizes that must be allocated
	    pmdata->alloc_invec1_size = 0;
	    pmdata->alloc_invec2_size = 0;
	    pmdata->alloc_outvec_size = 0;
	    pmdata->alloc_tmpvec_size = 0;

	    // Loop over all the other polys to determine the plan (brute/Karatsuba/FFT, FFT size, post_process_monic, 1*1 addins, etc.) needed for each
	    for (int i = 0; i < num_other_polys; i++) {
		gwnum	*invec2 = other_polys[i].invec2;		// Second input poly
		int	invec2_size = other_polys[i].invec2_size;	// Size of the second input polynomial
		gwnum	*outvec = other_polys[i].outvec;		// Output poly of size invec1_size + invec2_size - 1
		int	outvec_size = other_polys[i].outvec_size;	// Size of the output polynomial.  Should be invec1_size + invec2_size (less one if not monic).
		gwnum	*fmavec = other_polys[i].fmavec;		// Addin for FMA operations
		int	circular_size = other_polys[i].circular_size;	// Return result modulo (X^circular_size - 1)
		int	options;					// Merged invec1 and invec2 options to make life a little simpler
		bool	must_fft = FALSE;				// TRUE if there are any pre-FFTed poly inputs
		int	fft_size;

		// Apply any global invec2 options to all of the other polys.  Merge options into one variable.
		other_polys[i].options |= global_invec2_options;	// Merge with any non-poly1 global options
		options = invec1_options | other_polys[i].options;

		// Don't allow specifying any invec1 options in other poly array
		ASSERTG (!(other_polys[i].options & INVEC1_OPTIONS));
		// Cannot use combinations of POLYMULT_MULHI, POLYMULT_MULMID, and POLYMULT_MULLO
		ASSERTG (!(options & POLYMULT_MULHI) || !(options & (POLYMULT_MULMID | POLYMULT_MULLO)));
		ASSERTG (!(options & POLYMULT_MULMID) || !(options & (POLYMULT_MULHI | POLYMULT_MULLO)));
		ASSERTG (!(options & POLYMULT_MULLO) || !(options & (POLYMULT_MULHI | POLYMULT_MULMID)));
		// If POLYMULT_CIRCULAR option is consistent with circular_size setting
		ASSERTG (((options & POLYMULT_CIRCULAR) && circular_size > 0) || (!(options & POLYMULT_CIRCULAR) && circular_size == 0));
		// If first_mulmid can only be set when POLYMULT_MULMID is set
		ASSERTG ((options & POLYMULT_MULMID) || other_polys[i].first_mulmid == 0);
		// FMA vector cannot be pre-processed
		ASSERTG (fmavec == NULL || !is_preprocessed_poly (fmavec));

//GW: Make this a subroutine that polymult_preprocess can use???

		// Check for pre-processed other poly
		bool invec2_preprocessed = is_preprocessed_poly (invec2);
		bool invec2_pre_ffted = invec2_preprocessed && is_preffted_poly (invec2);
		exists_pre_ffted_invec2 |= invec2_pre_ffted;

		// Check for and validate pre-FFTed polys
		if (invec1_pre_ffted) {
			must_fft = TRUE;
			fft_size = preprocessed_fft_size (invec1);
		}
		if (invec2_pre_ffted) {
			ASSERTG (!must_fft || fft_size == preprocessed_fft_size (invec2));
			must_fft = TRUE;
			fft_size = preprocessed_fft_size (invec2);
		}

		// Adjust sizes due to RLPs and monics.  Roughly speaking, the adjusted_XXX variables are for monics with their ones stripped thus requiring
		// monic post processing, and the true_XXX variables are for monics handled by the poly FFT code.
		int adjusted_invec1_size = (options & POLYMULT_INVEC1_RLP) ? 2 * invec1_size - 1 : invec1_size;
		int adjusted_invec2_size = (options & POLYMULT_INVEC2_RLP) ? 2 * invec2_size - 1 : invec2_size;
		int adjusted_outvec_size = adjusted_invec1_size + adjusted_invec2_size - 1;
		int true_invec1_size = adjusted_invec1_size + ((options & POLYMULT_INVEC1_MONIC) ? ((options & POLYMULT_INVEC1_RLP) ? 2 : 1) : 0);
		int true_invec2_size = adjusted_invec2_size + ((options & POLYMULT_INVEC2_MONIC) ? ((options & POLYMULT_INVEC2_RLP) ? 2 : 1) : 0);
		int true_outvec_size = true_invec1_size + true_invec2_size - 1;

		// We now know the true input and output vector sizes.  THIS WILL NOT CHANGE.  We may do a brute/Karatsuba/FFT multiply that is smaller
		// than the true sizes and do some post-monic-addins to create the true result, but that will be reflected in other variables.
		// Save the true sizes to the plan.
		plan[i].true_invec1_size = true_invec1_size;
		plan[i].true_invec2_size = true_invec2_size;
		plan[i].true_outvec_size = true_outvec_size;

		// If circular size is larger than the polymult result, then the modulo operation is a no-op.
		if (circular_size >= true_outvec_size) options &= ~POLYMULT_CIRCULAR;
		// Set circular_size even if POLYMULT_CIRCULAR option is not set.  This just makes later code simpler.
		if (!(options & POLYMULT_CIRCULAR)) circular_size = true_outvec_size;

		// Monic * monic is a mine field of edge cases.  We may have several options of including the monic one in the poly FFT data or using
		// monic_line_adjustments after the polynomial multiplication.  Consider a monic-RLP times a monic-RLP, there are 4 possible implentations:
		//
		// Example:  1xxx1 * 1y1, impl1 is the obvious:  1xxx1*1y1
		// Impl2 is:	(xxx*y)<<1 + 1xxx1<<2 + 1xxx1 + 1y1<<4 + 1y1 but don't add the ones in twice
		// Impl3 is:	(1xxx1*y)<<1 + 1xxx1<<2 + 1xxx1
		// Impl4 is:	(xxx*1y1)<<1 + 1y1<<4 + 1y1
		//
		// Furthermore, polymult must not calculate 1*1 because one is not a random number.  Consider a brute force multiplication of two
		// polynomials that have one as a coefficient.  When two ones are multiplied together it adds the tiny signal FFT(1) to the result.
		// This is a worst case scenario for round off error as it is easy for the FFT(1) signal to be "drowned out" by the summing of dozens
		// or hundreds of multiplied signals of random numbers in the same "column".  There can be up to four 1*1 locations in the output
		// vector (monic RLP * monic RLP).  Calculate the locations of these possible 1*1s.  These locations will need to be "patched" later
		// and can affect the algorithm chosen to perform the polymult.
		int	addin1, addin2, addin3, addin4, subout;		// Locations of the up to five 1*1 values may need to be patched later
		addin1 = addin2 = addin3 = addin4 = subout = -1;	// Assume no addins or subouts of 1*1 are needed

		// We also need to track which of the true_outvec_size coefficients will be returned.  For example, a monic * monic will not return the
		// leading one.  A monic-RLP * monic-RLP returns a monic-RLP where the leading one and entire lower half of outvec are not returned.
		// Also, the caller can request a specific number of the most significant or least significant coefficients using MULHI/MULMID/MULLO.
		int	LSWs_skipped;		// Number of least significant coefficients that will not be returned to the caller
		int	MSWs_skipped;		// Number of most significant coefficients that will not be returned to the caller

		// In examples in the comments that follow, we show how two length two invecs can be processed

		// The simplest case is two non-monic inputs giving a non-monic result
		// Compute invec1 * invec2 as		(ab * cd)
		if (!(options & POLYMULT_INVEC1_MONIC) && !(options & POLYMULT_INVEC2_MONIC)) {
			LSWs_skipped = 0;
			MSWs_skipped = 0;
		}

		// Next up are monic * non-monic cases giving a non-monic result
		else if ((options & POLYMULT_INVEC1_MONIC) && !(options & POLYMULT_INVEC2_MONIC)) {
			// monic-RLP * non-monic
			// Compute invec1 * invec2 as		(aba * cd)		OR	(1aba1 * cd)
			// Then post_monic_add			(aba * cd)<<1 + cd<<4 + cd
			if (options & POLYMULT_INVEC1_RLP) {
				LSWs_skipped = 0;
				MSWs_skipped = 0;
			}
			// monic * non-monic
			// Compute invec1 * invec2 as		(ab * cd)		OR	(1ab * cd)
			// Then post_monic_add			(ab * cd) + cd<<2
			else {
				LSWs_skipped = 0;
				MSWs_skipped = 0;
			}
		}

		// and similarly the non-monic * monic cases giving a non-monic result
		else if (!(options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC2_MONIC)) {
			// non-monic * monic-RLP
			// Compute invec1 * invec2 as		(ab * cdc)		OR	(ab * 1cdc1)
			// Then post_monic_add			(ab * cdc)<<1 + ab<<4 + ab
			if (options & POLYMULT_INVEC2_RLP) {
				LSWs_skipped = 0;
				MSWs_skipped = 0;
			}
			// non-monic * monic
			// Compute invec1 * invec2 as		(ab * cd)		OR	(ab * 1cd)
			// Then post_monic_add			(ab * cd) + ab<<2
			else {
				LSWs_skipped = 0;
				MSWs_skipped = 0;
			}
		}

		// Both inputs are monic raising the possibility of forbidden 1*1 calculations
		// The monic * monic case gives us a monic result unless POLYMULT_CIRCULAR is set (leading 1*1 not output unless POLYMULT_CIRCULAR is set)
		// Compute invec1 * invec2 as			(ab * 1cd)		OR	(1ab * cd)		OR	(1ab * 1cd)
		// Then post_monic_add				(ab * 1cd) + cd<<2	OR	(1ab * cd) + ab<<2
		// One can also start with (ab * cd) and do two post monic adds:    (ab * cd) + cd<<2 + ab<<2
		// Then post add the 1*1 if necessary
		else if (!(options & POLYMULT_INVEC1_RLP) && !(options & POLYMULT_INVEC2_RLP)) {
			addin1 = true_outvec_size - 1;				// Location of the leading 1*1
			MSWs_skipped = 1;					// The leading 1*1 is not output
			LSWs_skipped = 0;
		}

		// The monic-RLP * monic-RLP case gives us a monic-RLP result unless POLYMULT_CIRCULAR is set
		// Compute invec1 * invec2 as		(aba * 1cdc1)				OR	(1aba1 * cdc)
		// Then post_monic_add			(aba * 1cdc1)<<1 + cdc<<5 + cdc		OR	(1aba1 * cdc)<<1 + aba<<5 + aba
		// One can also start with (aba * cdc) and do four post monic adds:	(aba * cdc)<<2 + cdc<<5 + cdc + aba<<5 + aba
		// Then post add two middle 1*1s OR post add four 1*1s if circular
		else if ((options & POLYMULT_INVEC1_RLP) && (options & POLYMULT_INVEC2_RLP)) {
			addin1 = true_outvec_size - 1;				// Location of the leading 1*1
			addin2 = true_invec2_size - 1;				// Location of middle 1*1 (leading 1 of invec2 * trailing 1 of invec1)
			addin3 = true_invec1_size - 1;				// Location of middle 1*1 (leading 1 of invec1 * trailing 1 of invec2)
			addin4 = 0;						// Location of the trailing 1*1
			MSWs_skipped = 1;					// The leading 1*1 is not output
			LSWs_skipped = 0;					// This will be changed later
		}

		// The monic-RLP * monic case gives us a monic result unless POLYMULT_CIRCULAR is set
		// Compute invec1 * invec2 as		(aba * 1cd)			OR	(1aba1 * cd)
		// Then post_monic_add			(aba * 1cd)<<1 + cd<<4 + cd		(1aba1 * cd) + aba<<3
		// One can also start with (aba * cd) and do three post monic adds:	(aba * cd)<<1 + cd<<4 + cd + aba<<3
		// Then post add one middle 1*1 OR post add both 1*1s if circular
		else if (!(options & POLYMULT_INVEC2_RLP)) {
			addin1 = true_outvec_size - 1;				// Location of the leading 1*1
			addin2 = true_invec2_size - 1;				// Location of middle 1*1 (leading 1 of invec2 * trailing 1 of invec1)
			MSWs_skipped = 1;					// The leading 1*1 is not output
			LSWs_skipped = 0;
		}

		// Similarly, the monic * monic-RLP case gives us a monic result unless POLYMULT_CIRCULAR is set
		// Compute invec1 * invec2 as		(ab * 1cdc1)			OR	(1ab * cdc)
		// Then post_monic_add			(ab * 1cdc1) + cdc<<3			(1ab * cdc)<<1 + ab<<4 + ab
		// One can also start with (ab * cdc) and do three post monic adds:	(ab * cdc)<<1 + cdc<<3 + ab<<4 + ab
		// Then post add one middle 1*1 OR post add both 1*1s if circular
		else if (!(options & POLYMULT_INVEC2_RLP)) {
			addin1 = true_outvec_size - 1;				// Location of the leading 1*1
			addin3 = true_invec1_size - 1;				// Location of middle 1*1 (leading 1 of invec1 * trailing 1 of invec2)
			MSWs_skipped = 1;					// The leading 1*1 is not output
			LSWs_skipped = 0;
		}

		// The RLP * RLP case gives us a RLP result unless POLYMULT_CIRCULAR is set.  Do not return the half of the coefficients.
		if ((options & POLYMULT_INVEC1_RLP) && (options & POLYMULT_INVEC2_RLP)) {
			LSWs_skipped = true_outvec_size / 2;			// The trailing reciprocal coefficients are not output
		}

		// For a circular polymult every coefficient is computed
		if (options & POLYMULT_CIRCULAR) {
			MSWs_skipped = 0;
			LSWs_skipped = 0;
		}

		// Validate outvec size argument
		ASSERTG ((!(options & (POLYMULT_MULHI | POLYMULT_MULMID | POLYMULT_MULLO)) && outvec_size == circular_size - MSWs_skipped - LSWs_skipped) ||
			 ((options & (POLYMULT_MULHI | POLYMULT_MULLO)) && outvec_size <= circular_size - MSWs_skipped - LSWs_skipped) ||
			 ((options & POLYMULT_MULMID) && outvec_size <= circular_size - MSWs_skipped - LSWs_skipped - other_polys[i].first_mulmid));

		// Adjust the output words skipped based on caller's MULHI/MULMID/MULLO request
		if (options & POLYMULT_MULHI) LSWs_skipped += circular_size - MSWs_skipped - LSWs_skipped - outvec_size;
		if (options & POLYMULT_MULMID) LSWs_skipped += other_polys[i].first_mulmid;
		if (options & POLYMULT_MULLO) MSWs_skipped += circular_size - MSWs_skipped - LSWs_skipped - outvec_size;

		// If emulating POLYMULT_CIRCULAR, apply modulo to every addin location
		if (addin1 >= circular_size) addin1 %= circular_size;
		if (addin2 >= circular_size) addin2 %= circular_size;
		if (addin3 >= circular_size) addin3 %= circular_size;
		if (addin4 >= circular_size) addin4 %= circular_size;

		// Make the 1*1 addins an index into outvec rather than the current index into true_outvec.  That is, adjust for LSWs skipped.
		// We do this because 1*1 addins are added as the very last step before returning outvec.
		addin1 -= LSWs_skipped;
		addin2 -= LSWs_skipped;
		addin3 -= LSWs_skipped;
		addin4 -= LSWs_skipped;

		// If any of the 1*1s are for a coefficient that is not returned to the caller, then clear the addin flag
		if (addin1 >= 0 && (addin1 >= outvec_size || outvec[addin1] == NULL)) addin1 = -1;
		if (addin2 >= 0 && (addin2 >= outvec_size || outvec[addin2] == NULL)) addin2 = -1;
		if (addin3 >= 0 && (addin3 >= outvec_size || outvec[addin3] == NULL)) addin3 = -1;
		if (addin4 >= 0 && (addin4 >= outvec_size || outvec[addin4] == NULL)) addin4 = -1;

		// Save addins and skipped outputs in the plan
		plan[i].addin[0] = addin1;
		plan[i].addin[1] = addin2;
		plan[i].addin[2] = addin3;
		plan[i].addin[3] = addin4;
		plan[i].LSWs_skipped = LSWs_skipped;
		plan[i].MSWs_skipped = MSWs_skipped;

		// If any 1*1s affect returned data then one (or both) of the monic polys must be stripped of their one coefficients and results
		// patched as a final step.  Prefer to strip the one coefficients of the longer input (so we patch by adding in the shorter input).
		// We track whether stripping ones from a monic invec is required.  Plus track if either invec can have its ones stripped without requiring
		// monic post processing (very rare, but an early prime95 P-1 implementation actually did this).
		bool	must_strip_a_monic;		// TRUE if at least one of the two invecs must have its monic ones stripped to avoid a dangerous 1*1
		bool	invec1_can_be_stripped_cheaply;	// TRUE if invec1 can have its ones stripped without requiring monic post processing
		bool	invec2_can_be_stripped_cheaply;	// TRUE if invec1 can have its ones stripped without requiring monic post processing
		bool	strip_monic_from_invec1, strip_monic_from_invec2, post_monic_addin_invec1, post_monic_addin_invec2;

		// If any 1*1 calculations land in a returned coefficient, then one of the monics must have its ones stripped.
		must_strip_a_monic = (addin1 >= 0 || addin2 >= 0 || addin3 >= 0 || addin4 >= 0);

		// Stripping the ones from a monic normally requires post processing that adds the other invec to the result.  However, if all the non-one
		// coefficients in the other invec land in columns that are not returned, then the post processing step can be skipped (at most just the
		// final 1*1 addin is needed to create the correct result).
		invec1_can_be_stripped_cheaply = (options & POLYMULT_INVEC1_MONIC) &&
			true_invec2_size - (((options & POLYMULT_INVEC2_MONIC) && (options & POLYMULT_INVEC2_RLP)) ? 1 : 0) <= MSWs_skipped &&
			(!(options & POLYMULT_INVEC1_RLP) || true_invec2_size - ((options & POLYMULT_INVEC2_MONIC) ? 1 : 0) <= LSWs_skipped);
		invec2_can_be_stripped_cheaply = (options & POLYMULT_INVEC2_MONIC) &&
			true_invec1_size - (((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC1_RLP)) ? 1 : 0) <= MSWs_skipped &&
			(!(options & POLYMULT_INVEC2_RLP) || true_invec1_size - ((options & POLYMULT_INVEC1_MONIC) ? 1 : 0) <= LSWs_skipped);

		// If a pre-FFTed input had its monic ones stripped then we obviously must strip a monic
		must_strip_a_monic |= (invec1_pre_ffted && !preprocessed_monics_included (invec1));
		must_strip_a_monic |= (invec2_pre_ffted && !preprocessed_monics_included (invec2));

		// Under ideal conditions, a circular polymult can be done for free by the poly FFT.  Otherwise it must be emulated.
		// Post-monic processing can handle a circular adjusted polymult result, but handling a shifted circular result would be tricky (invecs would
		// need shifting or don't delete the trailing ones, brute_line would need changing).  Since this tricky combination only occurs in the rare
		// monic * monic-RLP case, we simply force emulating circular size.
		bool	can_use_circular_fft_size = FALSE;		// TRUE if using circular_size smaller than the true_outvec_size as poly FFT size will work
		bool	can_use_circular_if_monic_stripped = FALSE;	// TRUE if using circular_size smaller than the true_outvec_size works if a monic is stripped
		int	circular_fft_size;
		if (options & POLYMULT_CIRCULAR) {
			// Verify that a pre-FFTed input vector is compatible with the circular size
			ASSERTG (!must_fft || (fft_size >= circular_size && (fft_size <= circular_size + LSWs_skipped + 1 || fft_size >= true_outvec_size)));
			// Calculate polymult FFT size we hope we can use for a free POLYMULT_CIRCULAR operation
			circular_fft_size = must_fft ? fft_size : polymult_fft_size (circular_size);
			// We can use circular_fft_size for a free POLYMULT_CIRCULAR when circular_fft_size equals circular size OR when the wrapped coefficients
			// pollute only the LSW_skipped coefficients (i.e. none of wrapped coefficients belonged in the final result).
			// This means true_outvec_size - circular_size <= LSWs_skipped
			// If we strip monics we can let one wrapped coefficient be in the final result: true_outvec_size - circular_size <= LSWs_skipped + 1
			if (circular_size >= true_invec1_size && circular_size >= true_invec2_size &&
			    circular_fft_size < true_outvec_size &&
			    (circular_fft_size == circular_size ||
			     (true_outvec_size - circular_size <= LSWs_skipped && !must_strip_a_monic) ||
			     (true_outvec_size - circular_size <= LSWs_skipped + 1 && (options & (POLYMULT_INVEC1_MONIC | POLYMULT_INVEC2_MONIC)) &&
			      !((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC2_MONIC) && (options & (POLYMULT_INVEC1_RLP | POLYMULT_INVEC2_RLP)))))) {
				can_use_circular_fft_size = TRUE;
				can_use_circular_if_monic_stripped = !(circular_fft_size == circular_size ||
								       (true_outvec_size - circular_size <= LSWs_skipped && !must_strip_a_monic));
			}
		}

		// If some of the coefficients are not returned, we can use a smaller poly FFT size as the long as the wrapped
		// coefficients pollute the coefficients that are not being returned.  For example, consider a true_outvec_size of 18:
		// 4 MSWs skipped, 10 returned, 4 LSWs skipped => A circular mult=14 works (0 MSWs skipped, 10 returned, 4 LSWs skipped)
		// 3 MSWs skipped, 10 returned, 5 LSWs skipped => A circular mult=15 works (0 MSWs skipped, 10 returned, 5 LSWs skipped)
		// 5 MSWs skipped, 10 returned, 3 LSWs skipped => A circular mult=15 works (2 MSWs skipped, 10 returned, 3 LSWs skipped)
		//GW NOTE:  We can do one better by writing monic_addin code that understands MSWs skipped.  That is a monic size9 * non-monic size9 (true_size=18)
		// 5 MSWs skipped, 9 returned, 4 LSWs skipped => Strip the monic9 and use circular mult=13 (0 MSWs skipped, 9 returned, 4 LSWs skipped
		// with only the lower 4 coefficients of the non-monic size9 added back in to the 9 returned coefficeints)
		else if (MSWs_skipped > 0 && LSWs_skipped > 0 && (!must_fft || fft_size < true_outvec_size)) {
			int savings = (LSWs_skipped <= MSWs_skipped) ? LSWs_skipped : MSWs_skipped;
			circular_fft_size = polymult_fft_size (true_outvec_size - savings);
			if (must_fft && fft_size > circular_fft_size) circular_fft_size = fft_size;
			if (circular_fft_size < true_outvec_size) can_use_circular_fft_size = TRUE;
//GW: Pre-fft knows nothing about this optimization -- need a method for pre-fft to output this smaller circular size??
		}

		// Figure out how we want to handle monics.  If implementation is POLYMULT_IMPL_FFT, post processing monics is an O(N) cost whereas
		// having the poly FFT do some or all of the monic work is free --- as long as it does not require a larger poly FFT size.
		// For other implementations, the choice between post-processing monics or handling as part of brute or Karatsuba is pretty much irrelevant.
		int	best_choice = 0;
		int	best_fft_size = 0;
		int	best_poly_addin_size = 0;

		// There are 4 possible monic implementations.  1) All monic ones included in FFT data.  2) Monic ones stripped from invec2 requiring invec1
		// to be added with post processing.  3) Monic ones stripped from invec1 requiring invec2 to be added with post processing.  4) Monic ones
		// stripped from both invec1 and invec2 requiring invec1 and invec2 to be added with post processing.
		// Select the implementation that gives us the smallest poly FFT size and fewest number of post-processing adds.

//GW: best monic impl may well depend on how circular_size is implemented.  worrying too much about circular impl with monic addins is fairly pointless

//GW:  Does code below work well for fft_size < FFT_BREAK?  we'll be selecting brute force impl.  Not sure if post_process is slower or faster.
//GW:  This is especially important for detecting different strip_monic_invec1 settings

		// If handling monics without stripping ones is possible, make that our best poly fft choice thusfar
		if (!must_strip_a_monic && (pass == 1 || !pass2_strip_monic_from_invec1)) {
					int possible_fft_size = polymult_fft_size (true_outvec_size);
			if (can_use_circular_fft_size && !can_use_circular_if_monic_stripped && circular_fft_size < possible_fft_size) possible_fft_size = circular_fft_size;
			if (must_fft && fft_size > possible_fft_size) possible_fft_size = fft_size;
			if (!must_fft && pass == 2 && possible_fft_size < largest_fft_size) possible_fft_size = largest_fft_size;
			if (!must_fft || possible_fft_size <= fft_size) {
				best_choice = 1;
				best_fft_size = possible_fft_size;
				best_poly_addin_size = 0;
			}
		}

		// An experimental case where we allow a single 1*1 calculation that we think will not trigger roundoff errors.  If wrapping just the 1*1
		// in a monic times monic polymult results in a smaller FFT size, then use the smaller FFT size and subtract 1*1 from the least significant
		// coefficient.  We think this will be OK because the poly FFT will be adding the 1*1 to a single subproduct (LSW-of-invec1 * LSW-of-invec2).
		if (!must_strip_a_monic && (pass == 1 || !pass2_strip_monic_from_invec1) &&
		    (options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC2_MONIC) && LSWs_skipped == 0) {
			int possible_fft_size = polymult_fft_size (true_outvec_size - 1);
			if ((must_fft || possible_fft_size >= pmdata->FFT_BREAK) &&
			    (best_choice == 0 || possible_fft_size < best_fft_size) &&
			    (pass == 1 || possible_fft_size == largest_fft_size) &&
			    (!must_fft || (must_fft && fft_size == possible_fft_size))) {
				best_choice = 2;
				best_fft_size = possible_fft_size;
				best_poly_addin_size = 0;
			}
		}

		// If handling monics by removing ones from invec2 and post-adding invec1 is possible, then see if that is a better choice
		if ((options & POLYMULT_INVEC2_MONIC) &&
		    (pass == 1 || !pass2_strip_monic_from_invec1) &&
		    (!invec2_pre_ffted || !preprocessed_monics_included (invec2)) &&
		    (invec2_can_be_stripped_cheaply || !invec1_pre_ffted)) {
			int possible_fft_size = polymult_fft_size (true_outvec_size - (true_invec2_size - adjusted_invec2_size));
			int possible_poly_addin_size = (invec2_can_be_stripped_cheaply ? 0 : adjusted_invec1_size);
			if (can_use_circular_fft_size && circular_fft_size < possible_fft_size) possible_fft_size = circular_fft_size;
			if (!must_fft && pass == 2 && possible_fft_size < largest_fft_size) possible_fft_size = largest_fft_size;
			if (must_fft && fft_size > possible_fft_size) possible_fft_size = fft_size;
			if ((!must_fft || possible_fft_size <= fft_size) &&
			    (best_choice == 0 || possible_fft_size < best_fft_size ||
			     (possible_fft_size == best_fft_size && possible_poly_addin_size < best_poly_addin_size))) {
				best_choice = 3;
				best_fft_size = possible_fft_size;
				best_poly_addin_size = possible_poly_addin_size;
			}
		}

		// If handling monics by removing ones from invec1 and post-adding invec2 is possible, then see if that is a better choice
		if ((options & POLYMULT_INVEC1_MONIC) &&
		    (pass == 1 || pass2_strip_monic_from_invec1) &&
		    (!invec1_pre_ffted || !preprocessed_monics_included (invec1)) &&
		    (invec1_can_be_stripped_cheaply || !invec2_pre_ffted)) {
			int possible_fft_size = polymult_fft_size (true_outvec_size - (true_invec1_size - adjusted_invec1_size));
			int possible_poly_addin_size = (invec1_can_be_stripped_cheaply ? 0 : adjusted_invec2_size);
			if (can_use_circular_fft_size && circular_fft_size < possible_fft_size) possible_fft_size = circular_fft_size;
			if (!must_fft && pass == 2 && possible_fft_size < largest_fft_size) possible_fft_size = largest_fft_size;
			if (must_fft && fft_size > possible_fft_size) possible_fft_size = fft_size;
			if ((!must_fft || possible_fft_size <= fft_size) &&
			    (best_choice == 0 || possible_fft_size < best_fft_size ||
			     (possible_fft_size == best_fft_size && possible_poly_addin_size < best_poly_addin_size))) {
				best_choice = 4;
				best_fft_size = possible_fft_size;
				best_poly_addin_size = possible_poly_addin_size;
			}
		}

		// If handling both monics by removing ones is possible, then see if that is a better choice
		if ((options & POLYMULT_INVEC1_MONIC) && (options & POLYMULT_INVEC2_MONIC) &&
		    (pass == 1 || pass2_strip_monic_from_invec1) &&
		    (!invec1_pre_ffted || !preprocessed_monics_included (invec1)) &&
		    (!invec2_pre_ffted || !preprocessed_monics_included (invec2)) &&
		    (invec1_can_be_stripped_cheaply || !invec2_pre_ffted) &&
		    (invec2_can_be_stripped_cheaply || !invec1_pre_ffted)) {
		    	int possible_fft_size = polymult_fft_size (adjusted_outvec_size);
			int possible_poly_addin_size = (invec1_can_be_stripped_cheaply ? 0 : adjusted_invec2_size) +
						       (invec2_can_be_stripped_cheaply ? 0 : adjusted_invec1_size);
			if (can_use_circular_fft_size && circular_fft_size < possible_fft_size) possible_fft_size = circular_fft_size;
			if (!must_fft && pass == 2 && possible_fft_size < largest_fft_size) possible_fft_size = largest_fft_size;
			if (must_fft && fft_size > possible_fft_size) possible_fft_size = fft_size;
			if ((!must_fft || possible_fft_size <= fft_size) &&
			    (best_choice == 0 || possible_fft_size < best_fft_size ||
			     (possible_fft_size == best_fft_size && possible_poly_addin_size < best_poly_addin_size))) {
				best_choice = 5;
				best_fft_size = possible_fft_size;
				best_poly_addin_size = possible_poly_addin_size;
			}
		}

		// Make sure one of the 5 implementations was possible
		ASSERTG (best_choice != 0);

		// Set flags and sizes based on the best monic implementation determined above
		if (best_choice == 1) {		// Monics handled without post processing
			strip_monic_from_invec1 = post_monic_addin_invec2 = FALSE;
			strip_monic_from_invec2 = post_monic_addin_invec1 = FALSE;
			adjusted_invec1_size = true_invec1_size;
			adjusted_invec2_size = true_invec2_size;
		}
		if (best_choice == 2) {		// Monics handled without post processing by wrapping a 1*1 and later subtracting it back out
			strip_monic_from_invec1 = post_monic_addin_invec2 = FALSE;
			strip_monic_from_invec2 = post_monic_addin_invec1 = FALSE;
			adjusted_invec1_size = true_invec1_size;
			adjusted_invec2_size = true_invec2_size;
			circular_size = best_fft_size;
//GW:  Do we need to set or adjust anything else?
			subout = 0;		// Subtract one from least significant coefficient
		}
		if (best_choice == 3) {		// Monics handled by strippng ones from invec2 and usually adding invec1 during post processing
			strip_monic_from_invec1 = post_monic_addin_invec2 = FALSE;
			strip_monic_from_invec2 = TRUE; post_monic_addin_invec1 = !invec2_can_be_stripped_cheaply;
			adjusted_invec1_size = true_invec1_size;
		}
		if (best_choice == 4) {		// Monics handled by strippng ones from invec1 and usually adding invec2 during post processing
			strip_monic_from_invec1 = TRUE; post_monic_addin_invec2 = !invec1_can_be_stripped_cheaply;
			strip_monic_from_invec2 = post_monic_addin_invec1 = FALSE;
			adjusted_invec2_size = true_invec2_size;
		}
		if (best_choice == 5) {		// Monics handled by stripping ones from both invecs
			strip_monic_from_invec1 = TRUE; post_monic_addin_invec2 = !invec1_can_be_stripped_cheaply;
			strip_monic_from_invec2 = TRUE; post_monic_addin_invec1 = !invec2_can_be_stripped_cheaply;
		}
		adjusted_outvec_size = adjusted_invec1_size + adjusted_invec2_size - 1;

		// Make sure every "other poly" uses the same FFT size if possible.  In pass 1 we find the largest FFT size.
		// This is not optimal but not worth optimizing, after all how often will someone want to have many other polys of different size!
		// Optimal would likely order the other polys by FFT size and then minimize the number of FFT sizes used.
		if (pass == 1) {
			if (best_fft_size != largest_fft_size) {
				// We need a second pass if we encounter multiple FFT sizes in the first pass
				if (largest_fft_size != 0) must_replan = TRUE;
				// Remember the largest FFT size encountered
				if (best_fft_size > largest_fft_size) {
					largest_fft_size = best_fft_size;
					largest_strip_monic_from_invec1 = strip_monic_from_invec1;
				}
			} else {
				// Remember if any usage of the largest FFT size requires stripping monics from invec1
				if (strip_monic_from_invec1 != largest_strip_monic_from_invec1) {
					largest_strip_monic_from_invec1 = TRUE;
					must_replan = TRUE;
				}
			}
		}
		// In pass 2, we should be using the largest FFT size whenever possible
		ASSERTG (pass == 1 || must_fft || best_fft_size == largest_fft_size);

		// Choose a brute force implementation if we are not returning many coefficients or one of the input polys is very small
		// NOTE: No studies have been done regarding a close-to-optimal definition of "very small".
		if (!must_fft && (outvec_size < pmdata->KARAT_BREAK || adjusted_invec1_size < pmdata->KARAT_BREAK || adjusted_invec2_size < pmdata->KARAT_BREAK)) {
			plan[i].impl = POLYMULT_IMPL_BRUTE;
		}
		// Choose Karatsuba
		else if (!must_fft && !can_use_circular_fft_size && adjusted_outvec_size < pmdata->FFT_BREAK) {
			plan[i].impl = POLYMULT_IMPL_KARATSUBA;
		}
		// FFT implementation.  There may be several possible FFT sizes to choose from.  We prefer to choose the smallest.
		else {
//GW: Does this code deal with the must_fft flag???  
			plan[i].impl = POLYMULT_IMPL_FFT;
			plan[i].fft_size = fft_size = best_fft_size;
		}						

		// If we selected the circular fft size, then save some necessary plan variables
		if (can_use_circular_fft_size && plan[i].impl == POLYMULT_IMPL_FFT && fft_size == circular_fft_size) {
			circular_size = circular_fft_size;
//GW:	MSWs_skipped = MSWs_skipped - LSWs_skipped; --- no, see monic_line_adjustments
// GW:  LSWs_skipped was already written to plan

//GW: Pre-fft knows nothing about this optimization -- need a method for pre-fft to output this smaller circular size
//GW: Rotation not allowed if both inputs pre-ffted
			// Verify that a pre-FFTed input vector is compatible with the circular size
//GW			ASSERTG (!must_fft || fft_size >= circular_size);
		}


		// Calculate how much to pad the adjusted polymult result.  Padding may be required on the left and right to reach true_outvec_size.
//		if (true_outvec_size != circular_size) {
//GW: BUG????
//do this before the brute/karat/fft above?
//plan[i].adjusted_pad = (int) post_monic_addin_invec1 + (int) post_monic_addin_invec2;		// Padding on the left
//plan[i].adjusted_shift = true_outvec_size - plan[i].adjusted_pad - adjusted_outvec_size;	// Padding on the right
//GW: this is wrong if poly FFT is implementing circular????? s.b. 0,0???
//		} else {
			plan[i].adjusted_pad = (int) post_monic_addin_invec1 + (int) post_monic_addin_invec2;		// Padding on the left
			plan[i].adjusted_shift = true_outvec_size - plan[i].adjusted_pad - adjusted_outvec_size;	// Padding on the right
//		}


		// Save adjusted sizes and monic addins in the plan
		plan[i].adjusted_invec1_size = adjusted_invec1_size;
		plan[i].adjusted_invec2_size = adjusted_invec2_size;
		plan[i].adjusted_outvec_size = adjusted_outvec_size;
		plan[i].post_monic_addin_invec1 = post_monic_addin_invec1;
		plan[i].post_monic_addin_invec2 = post_monic_addin_invec2;
		plan[i].strip_monic_from_invec1 = strip_monic_from_invec1;
		plan[i].strip_monic_from_invec2 = strip_monic_from_invec2;
		plan[i].subout = subout;

		// Set plan's emulate_circular flag.  Brute force and ideal circular FFT size both implement circular_size for free.
		plan[i].emulate_circular = ((options & POLYMULT_CIRCULAR) && circular_size != true_outvec_size &&
					    ((plan[i].impl == POLYMULT_IMPL_FFT && circular_size != fft_size) || plan[i].impl == POLYMULT_IMPL_KARATSUBA));
		plan[i].circular_size = circular_size;

//GW:			when given a choice pre-process should include the ones unless it forces a larger poly FFT size.
//GW:			unless poly 2 will be pre-ffted, then we MUST include the ones  (need a preprocess option for that? POLYMULT_OTHER_POLY_PRE_FFTED?)

//GW  ASSERT if pre-ffted invec1 or invec2 is the wrong FFT size (or if we have a choice of sizes and pre-ffted size is one of them then use it)

//GW: POST PROCESS to find cases where different FFT sizes are used and and switch all to larger one.  Could happen in non-power-of-two tree, but likely very rare:
//gw: tree == 769 breaks down to 385 and 384 length.  It's not clear that using fftlen=385(rounds up to 512) for both would be faster    

		// Track the largest buffer sizes that must be allocated
		if (plan[i].impl == POLYMULT_IMPL_BRUTE) {
			if (adjusted_invec1_size > pmdata->alloc_invec1_size) pmdata->alloc_invec1_size = adjusted_invec1_size;
			if (adjusted_invec2_size > pmdata->alloc_invec2_size) pmdata->alloc_invec2_size = adjusted_invec2_size;
			if (true_outvec_size > pmdata->alloc_outvec_size) pmdata->alloc_outvec_size = true_outvec_size;
		}
		else if (plan[i].impl == POLYMULT_IMPL_KARATSUBA) {
			if (adjusted_invec1_size > pmdata->alloc_invec1_size) pmdata->alloc_invec1_size = adjusted_invec1_size;
			if (adjusted_invec2_size > pmdata->alloc_invec2_size) pmdata->alloc_invec2_size = adjusted_invec2_size;
			if (true_outvec_size > pmdata->alloc_outvec_size) pmdata->alloc_outvec_size = true_outvec_size;
			if (adjusted_outvec_size + 32 > pmdata->alloc_tmpvec_size) pmdata->alloc_tmpvec_size = adjusted_outvec_size + 32; // May need extra for each recursion
		}
		else {
			// Invec1 will be read into tmpvec unless post_monic_addin_invec1 is set or multiple FFT sizes/stripping will be used
			if ((post_monic_addin_invec1 || fft_size != largest_fft_size || strip_monic_from_invec1 != largest_strip_monic_from_invec1) &&
			    true_invec1_size > pmdata->alloc_invec1_size)
				pmdata->alloc_invec1_size = true_invec1_size;
			// The FFT of invec1 will be stored in tmpvec
			if (fft_size > pmdata->alloc_tmpvec_size) pmdata->alloc_tmpvec_size = fft_size;
			// Invec2 will be read into outvec unless post_monic_addin_invec2 is set
			if (post_monic_addin_invec2 && true_invec2_size > pmdata->alloc_invec2_size) pmdata->alloc_invec2_size = true_invec2_size;
			// The FFT of invec2 as well as the polymult result will be stored in outvec
			if (fft_size + plan[i].adjusted_pad + plan[i].adjusted_shift > pmdata->alloc_outvec_size)
				pmdata->alloc_outvec_size = fft_size + plan[i].adjusted_pad + plan[i].adjusted_shift;
		}
	    }

	    // For the optional second pass, if any of the times we selected the largest FFT size required stripping monic one from invec1, then always strip.
	    pass2_strip_monic_from_invec1 = largest_strip_monic_from_invec1;
	    // For the optional second pass, if any of the other polys were pre-FFTed then assume we cannot strip monic one from invec1.
	    // If using the largest FFT size required stripping monics then we will have to use the next larger FFT size.
	    if (exists_pre_ffted_invec2) {
		    if (largest_strip_monic_from_invec1) largest_fft_size = polymult_fft_size (largest_fft_size + 1);
		    pass2_strip_monic_from_invec1 = FALSE;
	    }
	}
//GW		Change poly_preprocess to either fft_size of polymult_fft_size(insz1 + insz2 - 1) OR polymult_fft_size (outsz) if polymult circular set.
//GW		This offloads some of this monic decision crap to the caller.  And lets him manage the exact fft size in case we would have decided wrong fft size.

				
	// Process each of the other polys
	for (int i = 0; i < num_other_polys; i++) {
		gwnum	*invec2 = other_polys[i].invec2;
		int	invec2_size = other_polys[i].invec2_size;
		gwnum	*outvec = other_polys[i].outvec;
		int	outvec_size = other_polys[i].outvec_size;
		gwnum	*fmavec = other_polys[i].fmavec;
		int	circular_size = other_polys[i].circular_size;
		int	options = other_polys[i].options + invec1_options;

		// FFT invec2 if needed.  FFT the 7 zero padded FFT values in the header.  Even though we are not allowed to change the input vectors, we
		// can write to the 13 values values in the header reserved for use by the polymult library.
		if (! is_preprocessed_poly (invec2)) {
			for (int j = 0; j < invec2_size; j++) {
				if (invec2[j] == NULL) continue;
				gwfft (pmdata->gwdata, invec2[j], invec2[j]);
				if (pmdata->gwdata->ZERO_PADDED_FFT) zpad7_fft (pmdata->gwdata, invec2[j]);
			}
		}
		// Also FFT the optional FMA vector
		if (fmavec != NULL) {
			for (int j = 0; j < other_polys[i].outvec_size; j++) {
				if (fmavec[j] == NULL) continue;
				gwfft (pmdata->gwdata, fmavec[j], fmavec[j]);
				if (pmdata->gwdata->ZERO_PADDED_FFT) zpad7_fft (pmdata->gwdata, fmavec[j]);
			}
		}
	}

	// FFT invec1 if needed.  FFT the 7 zero padded FFT values in the header.  Even though we are not allowed to change the input vectors, we
	// can write to the 13 values values in the header reserved for use by the polymult library.
	if (! is_preprocessed_poly (invec1)) {
		for (int j = 0; j < invec1_size; j++) {
			if (invec1[j] == NULL) continue;
			gwfft (pmdata->gwdata, invec1[j], invec1[j]);
			if (pmdata->gwdata->ZERO_PADDED_FFT) zpad7_fft (pmdata->gwdata, invec1[j]);
		}
	}

	// Prepare for polymult in parallel
	pmdata->invec1 = invec1;
	pmdata->invec1_size = invec1_size;
	pmdata->other_polys = other_polys;
	pmdata->num_other_polys = num_other_polys;
	pmdata->options = invec1_options;
	pmdata->next_line = 0;

	// Fire up each helper thread
	pmdata->helpers_doing_polymult = TRUE;			// The helpers are doing polymult work, not user defined work
	polymult_launch_helpers (pmdata);

	// This thread works on the polymult too.  Select which optimized routines to use.
#ifdef X86_64
	if (pmdata->gwdata->cpu_flags & CPU_AVX512F) polymult_line_avx512 (pmdata); else
#endif
	if (pmdata->gwdata->cpu_flags & CPU_FMA3) polymult_line_fma (pmdata);
	else if (pmdata->gwdata->cpu_flags & CPU_AVX) polymult_line_avx (pmdata);
	else if (pmdata->gwdata->cpu_flags & CPU_SSE2) polymult_line_sse2 (pmdata);
#ifndef X86_64
	else polymult_line_dbl (pmdata);
#endif

	// Wait for all helper threads to complete
	polymult_wait_on_helpers (pmdata);
	pmdata->helpers_doing_polymult = FALSE;			// The helpers are not doing polymult work, user is free to have the threads do user work

	// Process all of the output polys.  UnFFT results and add in the needed 1*1 values.
	for (int i = 0; i < num_other_polys; i++) {
		gwnum	*outvec = other_polys[i].outvec;
		int	outvec_size = other_polys[i].outvec_size;
		int	options = other_polys[i].options + invec1_options;
		int	addin_value = 1;

		// If POLYMULT_FNMADD is set we will be subtracting the 1*1 value
		if (options & POLYMULT_FNMADD) addin_value = -1;

		// Pre-planning figured out which output entries (if any) need to get 1.0 added to the result because a monic poly is multiplied by a monic poly.
		// We used to add FFT(1) in monic_adjustments but this led to excessive round off error as FFT(1) is NOT a random signal.  It only
		// affects the least significant bits and is apt to get drowned out by the addition of many multiplied large signals.

		// Unfft all output coefficients
		for (int j = 0; j < outvec_size; j++) {
			if (outvec[j] == NULL) continue;
			if (pmdata->gwdata->ZERO_PADDED_FFT) zpad7_unfft (pmdata->gwdata, outvec[j]);
			FFT_state (outvec[j]) = (pmdata->gwdata->k == 1.0 ? FULLY_FFTed : FFTed_FOR_FMA);
			bool j_affected_by_addin = (j == plan[i].addin[0] || j == plan[i].addin[1] ||
						    j == plan[i].addin[2] || j == plan[i].addin[3] || j == plan[i].subout);
			if ((options & POLYMULT_NO_UNFFT) && !j_affected_by_addin) continue;
//GW: use gwsetaddin rather than gwaddsmall?
			gwunfft2 (pmdata->gwdata, outvec[j], outvec[j],
				  ((options & (POLYMULT_STARTNEXTFFT | POLYMULT_NEXTFFT)) && !j_affected_by_addin) ? GWMUL_STARTNEXTFFT : 0);
			if (j == plan[i].addin[0]) gwaddsmall (pmdata->gwdata, outvec[j], addin_value);
			if (j == plan[i].addin[1]) gwaddsmall (pmdata->gwdata, outvec[j], addin_value);
			if (j == plan[i].addin[2]) gwaddsmall (pmdata->gwdata, outvec[j], addin_value);
			if (j == plan[i].addin[3]) gwaddsmall (pmdata->gwdata, outvec[j], addin_value);
			if (j == plan[i].subout) gwaddsmall (pmdata->gwdata, outvec[j], -addin_value);
			if (options & POLYMULT_NEXTFFT) gwfft (pmdata->gwdata, outvec[j], outvec[j]);
		}
	}

	// Free array that planned each poly multiplication
	free (plan);
}




/*------------------------------------------------------------------------------------------------
|	Easy multithreading using polymult threads.   Example code!!
+------------------------------------------------------------------------------------------------*/

const int HELPER_OPCODE_COPY = 1;
const int HELPER_OPCODE_FFT = 2;
const int HELPER_OPCODE_UNFFT = 3;
const int HELPER_OPCODE_UNFFT_FFT = 4;

// Define a structure to pass all needed data to polymult_helper_example
struct poly_helper_example_data {
	pmhandle *pmdata;		// Polymult handle -- mainly to access pmdata->poly_mutex if needed
	int	opcode;			// Work type for poly_helper_example
	gwatomic counter;		// Counter that lets each thread work on independent data
	union {
		struct {		// Arguments to a poly_copy call
			gwnum	*invec;
			gwnum	*outvec;
			int	poly_size;
		} copy_args;
		struct {		// Arguments to a poly_fft_coefficients, poly_unfft_coefficients, poly_unfft_fft_coefficients call
			gwnum	*vec;
			int	poly_size;
		} fft_args;
	};
};
   
// Copy invec to outvec
void poly_copy (pmhandle *pmdata, gwnum *invec, gwnum *outvec, int poly_size) {
	struct poly_helper_example_data helper_data;
	// Copy arguments to helper_data so that every thread can access them
	helper_data.opcode = HELPER_OPCODE_COPY;			// Work type for polymult_helper_example to execute
	atomic_set (helper_data.counter, 0);				// Reset counter to no-work-done-yet state
	helper_data.pmdata = pmdata;					// Polymult structure (not needed for this simple example)
	helper_data.copy_args.invec = invec;				// Poly to copy
	helper_data.copy_args.outvec = outvec;				// Destination
	helper_data.copy_args.poly_size = poly_size;			// Size of the poly
	// Multithreaded poly copy
	pmdata->helper_callback = &poly_helper_example;
	pmdata->helper_callback_data = &helper_data;
	polymult_launch_helpers (pmdata);				// Launch helper threads
	poly_helper_example (0, pmdata->gwdata, &helper_data);		// Have this thread help too
	polymult_wait_on_helpers (pmdata);				// Wait for helper threads to finish
}

// UNFFT and/or FFT all the coefficients in a poly.  Why might one want to do this?  In P-1/ECM Stage 2, we start with a large number of size 1 polys.  We multiply
// pairs to create size 2 polys.  We multiply pairs again to create size 4 polys, and so on.  Say you are working with small numbers where the gwnum library
// cannot multithread (or poorly multithreads) gwunfft or gwfft calls.  Say your machine supports 16 threads and polymult is asked to multiply two size 2 polys.
// If polymult, fires up helper threads to gwfft inputs or gwunfft outputs, there are only 4 coefficients to work on leaving 12 threads sitting idle.  The solution,
// is to combine the large number size 2 polys into one gigantic poly and use the routines below to gwunfft and/or gwfft all the coefficients in one batch,
// keeping all 16 threads busy at once.  This is done in conjunction with the POLYMULT_NO_UNFFT polymult option.  Note that using poly_unfft_fft_coefficients
// in the P-1/ECM scenario is more efficient because the gwunfft followed immediately by gwfft is likely to find the gwnum still in the CPU caches.
void poly_fft_coefficients (pmhandle *pmdata, gwnum *vec, int poly_size) {
	struct poly_helper_example_data helper_data;
	// Copy arguments to helper_data so that every thread can access them
	helper_data.opcode = HELPER_OPCODE_FFT;				// Work type for polymult_helper_example to execute
	atomic_set (helper_data.counter, 0);				// Reset counter to no-work-done-yet state
	helper_data.pmdata = pmdata;					// Polymult structure (not needed for this simple example)
	helper_data.fft_args.vec = vec;					// Poly to fft coefficients
	helper_data.fft_args.poly_size = poly_size;			// Size of the poly
	// Multithreaded poly fft coefficients
	pmdata->helper_callback = &poly_helper_example;
	pmdata->helper_callback_data = &helper_data;
	polymult_launch_helpers (pmdata);				// Launch helper threads
	poly_helper_example (0, pmdata->gwdata, &helper_data);		// Have this thread help too
	polymult_wait_on_helpers (pmdata);				// Wait for helper threads to finish
}
void poly_unfft_coefficients (pmhandle *pmdata, gwnum *vec, int poly_size) {
	struct poly_helper_example_data helper_data;
	// Copy arguments to helper_data so that every thread can access them
	helper_data.opcode = HELPER_OPCODE_UNFFT;			// Work type for polymult_helper_example to execute
	atomic_set (helper_data.counter, 0);				// Reset counter to no-work-done-yet state
	helper_data.pmdata = pmdata;					// Polymult structure (not needed for this simple example)
	helper_data.fft_args.vec = vec;					// Poly to unfft coefficients
	helper_data.fft_args.poly_size = poly_size;			// Size of the poly
	// Multithreaded poly unfft coefficients
	pmdata->helper_callback = &poly_helper_example;
	pmdata->helper_callback_data = &helper_data;
	polymult_launch_helpers (pmdata);				// Launch helper threads
	poly_helper_example (0, pmdata->gwdata, &helper_data);		// Have this thread help too
	polymult_wait_on_helpers (pmdata);				// Wait for helper threads to finish
}
void poly_unfft_fft_coefficients (pmhandle *pmdata, gwnum *vec, int poly_size) {
	struct poly_helper_example_data helper_data;
	// Copy arguments to helper_data so that every thread can access them
	helper_data.opcode = HELPER_OPCODE_UNFFT_FFT;			// Work type for polymult_helper_example to execute
	atomic_set (helper_data.counter, 0);				// Reset counter to no-work-done-yet state
	helper_data.pmdata = pmdata;					// Polymult structure (not needed for this simple example)
	helper_data.fft_args.vec = vec;					// Poly to unfft/fft coefficients
	helper_data.fft_args.poly_size = poly_size;			// Size of the poly
	// Multithreaded poly unfft/fft coefficients
	pmdata->helper_callback = &poly_helper_example;
	pmdata->helper_callback_data = &helper_data;
	polymult_launch_helpers (pmdata);				// Launch helper threads
	poly_helper_example (0, pmdata->gwdata, &helper_data);		// Have this thread help too
	polymult_wait_on_helpers (pmdata);				// Wait for helper threads to finish
}


// Example helper routine
void poly_helper_example (
	int	helper_num,		// 0 = main thread, 1+ = helper thread num
	gwhandle *gwdata,		// Single-threaded and *MOSTLY* thread-safe gwdata (probably cloned) to use.  Each helper thread gets its own gwdata.
	void	*info)			// Can point to anything you choose.  In this example it points to struct polymult_helper_example_data
{
	struct poly_helper_example_data *helper_data = (struct poly_helper_example_data *) info;
	pmhandle *pmdata = helper_data->pmdata;

/* Copy a poly.  This example demonstrates using an atomic counter to get next available coefficient to work on. */

	if (helper_data->opcode == HELPER_OPCODE_COPY) {
		// Loop copying coefficients from invec to outvec
		gwnum *invec = helper_data->copy_args.invec;
		gwnum *outvec = helper_data->copy_args.outvec;
		for ( ; ; ) {
			// Get poly coefficient number we are to copy
			int coeff = atomic_fetch_incr (helper_data->counter);

			// If the other helper threads have completed copying invec to outvec, then we're done
			if (coeff >= helper_data->copy_args.poly_size) break;

			// Copy the gwnum
			gwcopy (gwdata, invec[coeff], outvec[coeff]);
		}
	}

/* FFT a poly's coefficients */

	else if (helper_data->opcode == HELPER_OPCODE_FFT) {
		// Loop FFTing coefficients in vec
		gwnum *vec = helper_data->fft_args.vec;
		for ( ; ; ) {
			// Get poly coefficient number to FFT
			int coeff = atomic_fetch_incr (helper_data->counter);

			// If the other helper threads have completed FFTing vec, then we're done
			if (coeff >= helper_data->fft_args.poly_size) break;

			// FFT the gwnum
			gwfft (gwdata, vec[coeff], vec[coeff]);
		}
	}

/* UNFFT a poly's coefficients */

	else if (helper_data->opcode == HELPER_OPCODE_UNFFT) {
		// Loop UNFFTing coefficients in vec
		gwnum *vec = helper_data->fft_args.vec;
		for ( ; ; ) {
			// Get poly coefficient number to UNFFT
			int coeff = atomic_fetch_incr (helper_data->counter);

			// If the other helper threads have completed UNFFTing vec, then we're done
			if (coeff >= helper_data->fft_args.poly_size) break;

			// UNFFT the gwnum
			gwunfft (gwdata, vec[coeff], vec[coeff]);
		}
	}

/* UNFFT then FFT a poly's coefficients */

	else if (helper_data->opcode == HELPER_OPCODE_UNFFT_FFT) {
		// Loop UNFFTing coefficients in vec
		gwnum *vec = helper_data->fft_args.vec;
		for ( ; ; ) {
			// Get poly coefficient number to UNFFT/FFT
			int coeff = atomic_fetch_incr (helper_data->counter);

			// If the other helper threads have completed UNFFT/FFTing vec, then we're done
			if (coeff >= helper_data->fft_args.poly_size) break;

			// UNFFT/FFT the gwnum
			gwunfft2 (gwdata, vec[coeff], vec[coeff], GWMUL_STARTNEXTFFT);
			gwfft (gwdata, vec[coeff], vec[coeff]);
		}
	}
}

#endif

