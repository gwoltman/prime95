/*----------------------------------------------------------------------
| Copyright 2020-2021 Mersenne Research, Inc.  All rights reserved
|
| Auxiliary routines to exponentiate gwnums
+---------------------------------------------------------------------*/

#ifndef _EXPONENTIATE_H
#define _EXPONENTIATE_H

/* This is used by C and C++ code.  If used in a C++ program, don't let the C++ compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

// Raise a gwnum to a power
void simple_exponentiate (gwhandle *gwdata, gwnum x, uint64_t power);

// Raise a gwnum to a power
void exponentiate (gwhandle *gwdata, gwnum x, uint64_t power);

// Raise a gwnum to a mpz power
void exponentiate_mpz (gwhandle *gwdata, gwnum x, mpz_t power);

#ifdef __cplusplus
}
#endif

#endif
