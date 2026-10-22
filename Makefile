# Builds mprime on Linux, FreeBSD and macOS: the same as "make" in linux/,
# which builds gwnum first.  See BUILDING.md, and linux/Makefile for the
# options (BITS=32, RELEASE=1, BOOST=...).
#
# This is a GNU make makefile.  On FreeBSD, BSDmakefile makes plain "make"
# run gmake (pkg install gmake).

all mprime:
	$(MAKE) -C linux

gwnum:
	$(MAKE) -C gwnum

clean distclean:
	$(MAKE) -C linux $@

.PHONY: all mprime gwnum clean distclean
