/*----------------------------------------------------------------------
| polymult.h
|
| This file contains the headers and definitions that are used in the
| polynomial multiplication built upon the gwnum library.
| 
|  Copyright 2021 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

#ifndef _POLYMULT_H
#define _POLYMULT_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

// These routines can be used before polymult_init.  The information is used to select a gwnum FFT size with sufficient safety_margin
// to do polymults of the desired size.

// Get the needed safety_margin required for an invec1_size by invec2_size polymult
double polymult_safety_margin (int invec1_size, int invec2_size);

// Get the FFT size that will be used for an n = invec1_size + invec2_size polymult
int polymult_fft_size (int n);


/* Callers of the polymult routines must first allocate a pmhandle (on the heap or stack) and pass it to all polymult routines. */
typedef struct pmhandle_struct pmhandle;

// Initialize a polymult handle
void polymult_init (
	pmhandle *pmdata,		// Handle for polymult library that needs initializing
	gwhandle *gwdata);		// Handle for gwnum FFT library

// Control how many threads to use during polymult.  This must be set prior to first call to polymult.  If gwnum library has a thread_callback
// routine set, the polymult library will use the same callback routine with action code 20 and 21. */
#define polymult_set_num_threads(h,n)	(h)->num_threads = (n)

// Terminate use of a polymult handle.  Free up memory.
void polymult_done (
	pmhandle *pmdata);		// Handle for polymult library


/* Multiply two polynomials.  It is safe to use input gwnums in the output vector.  For a normal polynomial multiply outvec_size must be */
/* invec1_size + invec2_size - 1.  For monic polynomial multiply, the leading input coefficients of 1 are omitted as is the leading 1 output */
/* coefficient -- thus outvec_size must be invec1_size + invec2_size. */
/* Entries in outvec can be NULL if caller is not interested in those output coefficients. */
/* The POLYMULT_FORWARD_FFT option allows the FFT of invec1 to be used over and over again, for outvec the caller must know the FFT size that will */
/* be used (see polymult_fft_size).  An FFTed input vector cannot be multiplied with a monic polynomial.  The downside to reusing input vector FFTs is */
/* they consumes more memory. */
void polymult (
	pmhandle *pmdata,		// Handle for polymult library
	gwnum	*invec1,		// First input poly
	int	invec1_size,		// Size of the first input polynomial (or fft size if POLYMULT_INVEC1_FFTED)
	gwnum	*invec2,		// Second input poly
	int	invec2_size,		// Size of the second input polynomial (or fft size if POLYMULT_INVEC2_FFTED)
	gwnum	*outvec,		// Output poly
	int	outvec_size,		// Size of the output polynomial (or fft_size if POLYMULT_FORWARD_FFT)
	int	options);

#define POLYMULT_INVEC1_MONIC	0x1	// Invec1 is a monic polynomial.  Leading coefficient of one is implied.
#define POLYMULT_INVEC2_MONIC	0x2	// Invec2 is a monic polynomial.  Leading coefficient of one is implied.
#define POLYMULT_INVEC1_RLP	0x4	// Invec1 is an RLP (Reciprocal Laurent Polynomial).  Needs only half the storage.
#define POLYMULT_INVEC2_RLP	0x8	// Invec2 is an RLP (Reciprocal Laurent Polynomial).  Needs only half the storage.
#define POLYMULT_INVEC1_NEGATE	0x10	// Invec1 coefficients should be negated.
#define POLYMULT_INVEC2_NEGATE	0x20	// Invec2 coefficients should be negated.
#define	POLYMULT_INVEC1_FFTED	0x40	// First polynomial has been forward FFTed
#define	POLYMULT_INVEC2_FFTED	0x80	// Second polynomial has been forward FFTed
#define	POLYMULT_FORWARD_FFT	0x100	// Forward FFT the first polynomial, the second polynomial is ignored
#define POLYMULT_STARTNEXTFFT	0x200	// Similar to GWMUL_STARTNEXTFFT.  Applied to all output coefficients.
//#define GWMUL_ADDINCONST	0x0400		/* Addin the optional gwsetaddin value to the multiplication result */
//#define GWMUL_MULBYCONST	0x0800		/* Multiply the final result by the small mulbyconst */
#define POLYMULT_INVEC1_MONIC_RLP	(POLYMULT_INVEC1_MONIC | POLYMULT_INVEC1_RLP)
#define POLYMULT_INVEC2_MONIC_RLP	(POLYMULT_INVEC2_MONIC | POLYMULT_INVEC2_RLP)
//GW: option to accept compressed poly inputs?
//GW: invec1_ffted and invec2_ffted not needed if it is part of the FFT_STATE.

/* The pmhandle structure containing all of polymult's "global" data. */
struct pmhandle_struct {
	gwhandle *gwdata;		// Handle for gwnum FFT library
	int	num_threads;		// Number of threads to use computing polymults
	gwevent work_to_do;		// Event to signal polymult helper threads there is work to do
	gwevent helpers_done;		// Event to signal polymult helper threads are done
	gwmutex	poly_mutex;		// Mutex to make polymult thread safe when multi-threading
	int	next_thread_num;	// Lets us generate a unique id for each helper thread
	int	next_line;		// Next line for a thread to process
	int	helpers_active;		// Count of helper threads that are still active
	bool	termination_in_progress; // Flag for helper threads to exit
	gwthread *thread_ids;		// Thread ids for the spawned threads
//	bool	will_use_polymult_forward_fft;   Affects GW_HEADER_SIZE (thus a gwinit option rather than a polymult option?)
	int	twiddles_initialized;	// size of the twiddle tables
	double	*twiddle_real;		// Radix-2 sin/cos table
	double	*twiddle_imag;
	double	*twiddle_real4;		// Radix-4 sin/cos table
	double	*twiddle_imag4;
	// Brute force from outvec_size [1..KARAT_BREAK), Karatsuba from [KARAT_BREAK..FFT_BREAK). FFT from [FFT_BREAK..infinity)
	int	KARAT_BREAK;		// Output vector size where we switch from brute force to Karatsuba
	int	FFT_BREAK;		// Output vector size where we switch from Karatsuba to FFTs
	// Arguments to the current polymult call.  Copied here so that helper threads can access the arguments.
	gwnum	*invec1;		// First input poly
	int	invec1_size;		// Size of the first input polynomial (or fft size if POLYMULT_INVEC1_FFTED)
	gwnum	*invec2;		// Second input poly
	int	invec2_size;		// Size of the second input polynomial (or fft size if POLYMULT_INVEC2_FFTED)
	gwnum	*outvec;		// Output poly
	int	outvec_size;		// Size of the output polynomial (or fft_size if POLYMULT_FORWARD_FFT)
	int	options;
};


#ifdef __cplusplus
}
#endif

#endif
