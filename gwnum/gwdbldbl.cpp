/**************************************************************
 *
 *  gwdbldbl.cpp
 *
 *  This file contains all the gwnum initialization routines that require
 *  extended precision floating point.  We want to initialize our sin/cos
 *  and FFT weight arrays with doubles that are as accurate as possible.
 *
 *  This is the only C++ routine in the gwnum library.  Since gwnum is
 *  a C based library, we declare all routines here as extern "C".
 * 
 *  Copyright 2005 Just For Fun Software, Inc.
 *  All Rights Reserved.
 *
 **************************************************************/

/* Include files */

#include "gwdbldbl.h"

/* Pick which doubledouble package we will use. */

#define QD
//#define KEITH_BRIGGS

/* Turn on #define that will disable extended precision floating point */
/* registers, as they wreak havoc with double-double library routines. */

#ifndef X86_64
#define x86
#endif

/* Use Hida, Li & Bailey's QD doubledouble C++ package. */

#ifdef QD
#include "dd.cc"
#endif

/* Use Keith Briggs' doubledouble C++ package.  I find the QD package */
/* a better choice. */

#ifdef KEITH_BRIGGS
#define DD_INLINE
#include "doubledouble.cc"
#include "math.cc"
#define	dd_real doubledouble
#define	_2pi TwoPi
#define	_log2 Log2
#endif


/* Now write all the routines that use the dd_real package. */

// Utility routine to compute a sin/cos premultiplier or a set of 3
// sine-cosine values.
// This is used during setup, formerly written in assembly language to take
// advantage of the extra precision in the FPU's 80-bit registers. 
// NOTE: When computing cosine / sine, divide by the 64-bit sine
// not the 80-bit sine since macros will multiply by the 64-bit sine.

extern "C"
void gwsincos (
	unsigned long x,
	unsigned long N,
	double	*results)
{
	dd_real arg, sine, cosine;

	x86_FIX
	arg = dd_real::_2pi * (double) x / (double) N;
	sincos (arg, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[0] = sine;
	results[1] = cosine / results[0];
	END_x86_FIX
}

extern "C"
void gwsincos3 (
	unsigned long x,
	unsigned long N,
	double	*results)
{
	dd_real arg1, arg2, arg3, sine, cosine;

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[0] = sine;
	results[1] = cosine / results[0];
	arg2 = arg1 + arg1;
	sincos (arg2, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[2] = sine;
	results[3] = cosine / results[2];
	arg3 = arg2 + arg1;
	sincos (arg3, sine, cosine);
	sine += 1E-200;			/* Protect against divide by zero */
	results[4] = sine;
	results[5] = cosine / results[4];
	END_x86_FIX
}


//
// Utility routines to compute fft weights
//
// The FFT weight for the j-th FFT word doing a 2^q+c weighted transform is
//	2 ^ (ceil (j*q/FFTLEN) - j*q/FFTLEN)   *    abs(c) ^ j/FFTLEN
//

static	dd_real gw__bits_per_word;
static	int gw__c_is_one;
static	dd_real gw__log2_abs_c_div_fftlen;
static	dd_real gw__fftlen_inverse;
static	dd_real gw__over_fftlen;

extern "C"
void gwfft_weight_setup (
	int	zero_pad,
	double	k,
	unsigned long n,
	signed long c,
	unsigned long fftlen)
{
	x86_FIX
	gw__fftlen_inverse = dd_real (1.0) / dd_real ((double) fftlen);
	if (zero_pad) {
		gw__bits_per_word =
			dd_real ((double) (n + n)) * gw__fftlen_inverse;
		gw__c_is_one = 1;
		gw__over_fftlen = dd_real (2.0) * gw__fftlen_inverse;
	} else {
		gw__bits_per_word =
			(dd_real ((double) n) +
			 log (dd_real (k)) / dd_real::_log2) *
			gw__fftlen_inverse;
		gw__c_is_one = (abs ((int) c) == 1);
		gw__log2_abs_c_div_fftlen =
			log (dd_real (abs ((int) c))) / dd_real::_log2 *
			gw__fftlen_inverse;
		gw__over_fftlen = dd_real (k * 2.0) * gw__fftlen_inverse;
	}
	END_x86_FIX
}

extern "C"
double gwfft_weight (
	unsigned long j)
{
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_real::_log2 * twopow);
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.

extern "C"
double gwfft_weight_sloppy (
	unsigned long j)
{
	dd_real temp, twopow;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	END_x86_FIX
	return (pow (2.0, double (twopow)));
}

// Compute the inverse of the fft weight

extern "C"
double gwfft_weight_inverse (
	unsigned long j)
{
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_real::_log2 * -twopow);
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.

extern "C"
double gwfft_weight_inverse_sloppy (
	unsigned long j)
{
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	END_x86_FIX
	return (pow (2.0, - double (twopow)));
}

// This computes the inverse FFT weight multiplied by the appropriate constant
// to produce an integer during an FFT multiply's normalize stage.  This
// constant is 2/FFTLEN for a zero-padded FFT and k*2/FFTLEN for a
// non-zero-padded FFT.

extern "C"
double gwfft_weight_inverse_over_fftlen (
	unsigned long j)
{
	dd_real temp, twopow, result;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	result = exp (dd_real::_log2 * -twopow) * gw__over_fftlen;
	END_x86_FIX
	return (double (result));
}

// This computes the three FFT weights in one call.  It is faster than
// calling the above individually.

extern "C"
void gwfft_weights3 (
	unsigned long j,
	double	*fft_weight,
	double	*fft_weight_inverse,
	double	*fft_weight_inverse_over_fftlen)
{
	dd_real temp, twopow, weight;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	if (! gw__c_is_one)
		twopow += gw__log2_abs_c_div_fftlen * dd_real ((double) j);
	weight = exp (dd_real::_log2 * twopow);
	*fft_weight = double (weight);
	weight = dd_real (1.0) / weight;
	if (fft_weight_inverse != NULL)
		*fft_weight_inverse = double (weight);
	if (fft_weight_inverse_over_fftlen != NULL)
		*fft_weight_inverse_over_fftlen = double (weight * gw__over_fftlen);
	END_x86_FIX
}

extern "C"
double gwfft_weight_exponent (
	unsigned long j)
{
	dd_real temp, twopow;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = ceil (temp) - temp;
	END_x86_FIX
	return (double (twopow));
}

//
// Utility routine to compute fft base for j-th fft word
//
// The FFT base for the j-th FFT word doing a 2^q+c weighted transform is
//	ceil (j*q/FFTLEN)
// This routine returns ceil (j*q/FFTLEN) taking great care to return a
// value accurate to 53 bits.  This is important when j*q is really close to
// a multiple of FFTLEN (admittedly quite rare).  It would be really bad if
// rounding differences caused this routine to compute ceil (j*q/FFTLEN)
// differently than the weighting functions.
//

extern "C"
unsigned long gwfft_base (
	unsigned long j)
{
	dd_real temp;
	unsigned long twopow;

	x86_FIX
	temp = dd_real ((double) j) * gw__bits_per_word;
	twopow = (int) ceil (temp);
	END_x86_FIX
	return (twopow);
}
