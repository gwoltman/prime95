/*----------------------------------------------------------------------
| This file contains various utility routines that may be used by gwnum
| routines, prime95, or PRP.
+---------------------------------------------------------------------*/

#ifndef _GWUTIL_H
#define _GWUTIL_H

/* Aligned malloc routines.  MSVC 8 supports these in the C runtime library. */
/* Emulate these routines for other ports. */

void * aligned_offset_malloc (size_t size, size_t alignment, size_t mod);
void * aligned_malloc (size_t size, size_t alignment);
void  aligned_free (void *ptr);

#endif
