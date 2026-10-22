# FreeBSD's make reads this file before Makefile.  Makefile needs GNU make,
# so this runs gmake with the same targets, variables and -j.  MAKEFLAGS is
# cleared because it holds FreeBSD make's own flags, which gmake rejects.

GMAKE ?= gmake
_GMAKE = env -u MAKEFLAGS -u MFLAGS ${GMAKE} ${.MAKE.JOBS:D-j${.MAKE.JOBS}} ${.MAKEOVERRIDES:O:u:@v@${v}=${${v}:Q}@}

all: .PHONY
	@${_GMAKE}

.DEFAULT: .PHONY
	@${_GMAKE} ${.TARGET}
