/*----------------------------------------------------------------------
| gwdbldbl.h
|
| This file contains the headers for the gwnum helper routines that use
| extended-precision floats.
| 
|  Copyright 2005 Just For Fun Software, Inc.
|  All Rights Reserved.
+---------------------------------------------------------------------*/

#ifndef _GWDBLDBL_H
#define _GWDBLDBL_H

/* Include common definitions */

#include "common.h"

/* Extended precision helper routines */

EXTERNC void gwfft_weight_setup (int, double, unsigned long, signed long, unsigned long);
EXTERNC double gwfft_weight (unsigned long);
EXTERNC double gwfft_weight_inverse (unsigned long);
EXTERNC double gwfft_weight_inverse_over_fftlen (unsigned long);
EXTERNC void gwfft_weights3 (unsigned long, double *, double *, double *);
EXTERNC double gwfft_weight_exponent (unsigned long);
EXTERNC unsigned long gwfft_base (unsigned long);
EXTERNC void gwsincos (unsigned long, unsigned long, double *);
EXTERNC void gwsincos3 (unsigned long, unsigned long, double *);

#endif
