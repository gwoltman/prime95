/*----------------------------------------------------------------------
| common.h
|
| This file contains handy #defines that I use in all my projects
| 
|  Copyright 2005 Just For Fun Software, Inc.
|  All Rights Reserved.
+---------------------------------------------------------------------*/

#ifndef _COMMON_H
#define _COMMON_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#define TRUE	1
#define FALSE	0

/* Define the ASSERT macro I use while debugging */

#ifdef GDEBUG
#include <assert.h>
#define ASSERTG		assert
#else
#define ASSERTG(a)
#endif

/* If you compile the gwnum library without STPCHK defined, then the */
/* code that lets prime95 ESC out of a long GCD operation will not */
/* be compiled in. */

#ifdef STPCHK
EXTERNC int stopCheck();
#else
#define stopCheck()	FALSE
#endif

/* When debugging gwnum and giants, I sometimes write code that "cheats" */
/* by calling a routine that is part of prime95 rather than the gwnum and */
/* giants library.  Make that forward declaration here. */

EXTERNC void OutputBoth (char *);

#endif
