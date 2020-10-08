/*----------------------------------------------------------------------
| radix.h
|
| This file contains the C routine definitions for radix conversion
| when required by gianttogw or gwtogiant.
| 
|  Copyright 2020 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

void nonbase2_gianttogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	g,		/* Input giant */
	gwnum	x);		/* Output gwnum */


