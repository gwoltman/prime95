; Copyright 1995-2004 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements creates some unused data items to get
; segments aligned on 32-byte boundaries in Linux.
;

	TITLE   setup

	.386

_TEXT32 SEGMENT PAGE USE32 PUBLIC 'DATA'
_TEXT32 ENDS

END
