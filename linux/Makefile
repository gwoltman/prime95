# Makefile for mprime on Linux (64-bit and 32-bit), FreeBSD (64-bit) and
# macOS (Intel)
#
#	make		builds mprime for this machine, building gwnum first
#	make BITS=32	builds 32-bit Linux mprime on 64-bit Linux
#	make RELEASE=1	links the way official releases do (see below)
#	make clean	deletes what make built here and in gwnum, except the
#			assembled code
#	make distclean	deletes that too
#
# Needs GMP, hwloc, libcurl and the Boost headers (see BUILDING.md).  BOOST
# names the Boost directory if it is not where the compiler looks already:
#	make BOOST=~/boost_1_76_0
# In a git clone, gwnum's and factor's assembly code is assembled with uasm
# and objconv; set UASM and OBJCONV to use them from somewhere other than the
# PATH.
#
# By default mprime links with the shared hwloc and curl libraries.  Official
# releases (RELEASE=1) link those statically, so that mprime runs on more
# systems, and print a link map.  On Linux that needs static builds of them,
# for example hwloc configured with
#	--enable-static --disable-shared --disable-io
# and libcurl configured with
#	--without-ssl --without-libssh2 --without-libidn --disable-ldap --disable-ldaps --without-gssapi --enable-ipv6 --without-krb4
# Official Linux releases are built on CentOS 7: a program linked against a
# newer glibc will not run on older Linux versions.  On macOS, RELEASE=1 links
# /usr/local/lib/libhwloc.a and libgmp.a.
#
# This is a GNU make makefile.  On FreeBSD, BSDmakefile makes plain "make"
# run gmake (pkg install gmake).

include ../gwnum/platform.mk

BOOST = ../../boost_1_76_0
GWLIBS = ../gwnum/gwnum.a ../gwnum/polymult.a

ifeq ($(PLATFORM),linux64)
CC = gcc
CFLAGS = -I.. -I../gwnum -I/usr/local/include -std=gnu99 -DX86_64 -Wno-unused-result -O2 -DPORT=8
CXX = g++
CXXFLAGS = -I.. -I../gwnum -I$(BOOST) -DX86_64 -std=c++11 -O2 -DPORT=8
LDFLAGS = -Wl,-L/usr/local/lib
ifdef RELEASE
LDFLAGS = -Wl,-M -Wl,-L/usr/local/lib
LIBS = -lm -lpthread -Wl,-Bstatic -lhwloc -Wl,-Bstatic -lcurl -Wl,-Bdynamic -lrt -lstdc++ -ldl -lgmp
else
LIBS = -lm -lpthread -lhwloc -lcurl -lrt -lstdc++ -ldl -lgmp
endif
FACTOR = factor64.o
endif

ifeq ($(PLATFORM),linux32)
CC = gcc
CFLAGS = -I.. -I../gwnum -I/usr/local/include -std=gnu99 -Wno-unused-result -march=i486 -malign-double -O2 -DPORT=2 $(M32)
CXX = g++
CXXFLAGS = -I.. -I../gwnum -I$(BOOST) -Wno-unused-result -march=i486 -malign-double -std=c++11 -O2 -DPORT=2 $(M32)
LDFLAGS = -Wl,-L/usr/local/lib $(M32)
ifdef RELEASE
LDFLAGS = -Wl,-M -Wl,-L/usr/local/lib $(M32)
LIBS = -lm -lpthread -Wl,-Bstatic -lhwloc -Wl,-Bstatic -lcurl -Wl,-Bdynamic -lrt -lstdc++ -ldl -lgmp -latomic
else
LIBS = -lm -lpthread -lhwloc -lcurl -lrt -lstdc++ -ldl -lgmp -latomic
endif
FACTOR = factor32.o
endif

ifeq ($(PLATFORM),freebsd64)
# gcc, not FreeBSD's default clang: a clang build of mprime crashes at once.
# The C++ library is linked statically because it is gcc's, not the clang
# one FreeBSD installs.
CC = gcc
CFLAGS = -O2 -pipe -I.. -I../gwnum -I$(LOCALBASE)/include -std=gnu99 -DX86_64 -Wno-unused-result -O2 -DPORT=12
CXX = g++
CXXFLAGS = -O2 -pipe -I.. -I../gwnum -I$(LOCALBASE)/include -I$(BOOST) -DX86_64 -Wno-unused-result -std=c++11 -DPORT=12
LDFLAGS = -fpie -L$(LOCALBASE)/lib
ifdef RELEASE
LDFLAGS = -fpie -L$(LOCALBASE)/lib -Wl,-M
endif
LIBS = -lm -lpthread -lhwloc -lcurl -Wl,-Bstatic -lstdc++ -Wl,-Bdynamic -lcompat -lgmp
FACTOR = factor64.o
endif

ifeq ($(PLATFORM),macos64)
CC = gcc
CFLAGS = -I.. -I../gwnum -I../linux -std=gnu99 -O2 -DX86_64 -DCOMMAND_LINE_MPRIME -m64 -DPORT=10
CXX = g++
CXXFLAGS = -I.. -I../gwnum -I../linux -I$(BOOST) -std=c++11 -O2 -DX86_64 -DCOMMAND_LINE_MPRIME -m64 -DPORT=10
LDFLAGS = -m64 -Wl,-no_pie -L/usr/local/lib
ifdef RELEASE
LIBS = -lm -lpthread /usr/local/lib/libhwloc.a /usr/local/lib/libgmp.a -lcurl -framework IOKit -framework CoreFoundation -lc++
else
LIBS = -lm -lpthread -lhwloc -lgmp -lcurl -framework IOKit -framework CoreFoundation -lc++
endif
FACTOR = ../prime95/macosx64/factor64.o
endif

OBJS = $(addprefix $(OBJDIR)/,prime.o menu.o cJSON.o ecm.o exponentiate.o pair.o pm1prob.o)
vpath %.c .. ../third_party/cjson
vpath %.cpp ..

# Builds without the private security files (all builds but official ones)
# use empty ones and identify themselves as "Untrusted Prime95"
SECFILES = ../security.h ../security.c ../secure5.c

#########################################################################

mprime: $(OBJS) $(FACTOR) $(GWLIBS) $(PLATSTAMP)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(FACTOR) $(GWLIBS) $(LIBS)

$(OBJDIR)/%.o: %.c | $(SECFILES)
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.cpp | $(SECFILES)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(SECFILES):
	touch $@

-include $(wildcard $(OBJDIR)/*.d)

# gwnum: its own makefile decides whether anything needs rebuilding, and
# mprime is relinked only if gwnum.a or polymult.a actually changed
../gwnum/gwnum.a: FORCE
	$(MAKE) -C ../gwnum
../gwnum/polymult.a: ../gwnum/gwnum.a

#########################################################################
# factor, the trial factoring assembly code in prime95/.  As with gwnum,
# uasm assembles it into a COFF object and objconv converts that to ELF or
# Mach-O.

../prime95/linux64/factor64.obj: ../prime95/factor64.asm ../prime95/factor64.mac ../gwnum/unravel.mac
	@command -v $(UASM) > /dev/null || { echo "Assembling factor64 needs uasm, which is not on PATH; see BUILDING.md (or set UASM=/path/to/uasm)" >&2; exit 1; }
	mkdir -p ../prime95/linux64
	cd ../prime95 && $(UASM) -c -nologo -DX86_64 -DLINUX64 -DARCH=FMA3 -I../gwnum -win64 -safeseh -Folinux64/factor64.obj factor64.asm
../prime95/linux64/factor64.o: ../prime95/linux64/factor64.obj
	$(OBJCONV) -felf64 $< $@
../prime95/macosx64/factor64.o: ../prime95/linux64/factor64.obj
	mkdir -p ../prime95/macosx64
	$(OBJCONV) -fmacho64 -nu+ $< $@

../prime95/linux/factor32.obj: ../prime95/factor32.asm ../prime95/factor32.mac ../gwnum/unravel.mac
	@command -v $(UASM) > /dev/null || { echo "Assembling factor32 needs uasm, which is not on PATH; see BUILDING.md (or set UASM=/path/to/uasm)" >&2; exit 1; }
	mkdir -p ../prime95/linux
	cd ../prime95 && $(UASM) -c -nologo -coff -safeseh -I../gwnum -Folinux/factor32.obj factor32.asm
../prime95/linux/factor32.o: ../prime95/linux/factor32.obj
	$(OBJCONV) -felf32 $< $@

# Linux and FreeBSD link a copy here, as factor64.o or factor32.o: the
# linker records that name in mprime
factor64.o: ../prime95/linux64/factor64.o
	cp $< $@
factor32.o: ../prime95/linux/factor32.o
	cp $< $@

#########################################################################

clean:
	rm -rf build mprime factor64.o factor32.o
	$(MAKE) -C ../gwnum clean

distclean:
	rm -rf build mprime factor64.o factor32.o ../prime95/linux64 ../prime95/linux ../prime95/macosx64
	$(MAKE) -C ../gwnum distclean

FORCE:

.PHONY: clean distclean FORCE
.SECONDARY:
