; Copyright 1995-2005 Just For Fun Software, Inc., all rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements creates some unused data items to get
; segments aligned on 32-byte boundaries in Linux.
;

	TITLE   setup

	.386
	.MODEL	FLAT

_TEXT SEGMENT
dummy1a		DD	0
dummy2a		DD	0
dummy3a		DD	0
dummy4a		DD	0
dummy5a		DD	0
dummy6a		DD	0
dummy7a		DD	0
_TEXT ENDS

END
