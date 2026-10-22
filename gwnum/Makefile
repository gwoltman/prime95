# Makefile for the gwnum library on Linux (64-bit and 32-bit), FreeBSD
# (64-bit) and macOS (Intel)
#
#	make		builds gwnum.a and polymult.a for this machine
#	make BITS=32	builds them for 32-bit Linux on 64-bit Linux
#	make asm	assembles just the assembly code
#	make clean	deletes what make built, except the assembled code
#	make distclean	deletes that too (only when the assembly source is here)
#
# Programs that use gwnum link with gwnum.a and polymult.a.
#
# gwnum.a is a copy of the assembled code, linux64/gwnum.a (64-bit Linux and
# FreeBSD), linux/gwnum.a (32-bit Linux) or macosx64/gwnum.a, with the C code
# added.  In a git clone, make assembles that code first, with uasm and
# objconv (see BUILDING.md); set UASM and OBJCONV to use them from somewhere
# other than the PATH.  The gwnum package comes with the assembled code and
# no assembly source, so there make needs neither.
#
# This is a GNU make makefile.  On FreeBSD, BSDmakefile makes plain "make"
# run gmake (pkg install gmake).

include platform.mk

AR = ar

ifeq ($(PLATFORM),linux64)
CC = gcc
CFLAGS = -I.. -I/usr/local/include -std=gnu99 -DX86_64 -Wno-unused-result -O2
CXX = g++
CXXFLAGS = -I.. -std=c++11 -DX86_64 -O2
ASMLIB = linux64/gwnum.a
endif
ifeq ($(PLATFORM),linux32)
CC = gcc
CFLAGS = -I.. -I/usr/local/include -std=gnu99 -Wno-unused-result -march=i486 -malign-double -O2 $(M32)
CXX = g++
CXXFLAGS = -I.. -std=c++11 -O2 -march=i486 -malign-double $(M32)
ASMLIB = linux/gwnum.a
endif
ifeq ($(PLATFORM),freebsd64)
# gcc, not FreeBSD's default clang: a clang build of mprime crashes
CC = gcc
CFLAGS = -O2 -pipe -I.. -I$(LOCALBASE)/include -std=gnu99 -DX86_64
CXX = g++
CXXFLAGS = -O2 -pipe -I.. -I$(LOCALBASE)/include -DX86_64 -std=c++11
ASMLIB = linux64/gwnum.a
endif
ifeq ($(PLATFORM),macos64)
CC = gcc
CFLAGS = -I.. -O2 -std=gnu99 -DX86_64 -m64 -fno-exceptions
CXX = g++
CXXFLAGS = -I.. -O2 -DX86_64 -std=c++11 -m64 -fno-exceptions -fno-rtti -Wno-stdlibcxx-not-found
ASMLIB = macosx64/gwnum.a
endif

COBJS = cpuid gwnum gwtables gwthread gwini gwbench gwutil gwdbldbl giants radix
# polymult.c compiled once plainly and once for each instruction set
POLYOBJS = polymult polymult_sse2 polymult_avx polymult_fma
ifneq ($(PLATFORM),linux32)
POLYOBJS += polymult_avx512
endif
COBJFILES = $(addprefix $(OBJDIR)/,$(addsuffix .o,$(COBJS)))
POLYOBJFILES = $(addprefix $(OBJDIR)/,$(addsuffix .o,$(POLYOBJS)))

#########################################################################

all: gwnum.a polymult.a

gwnum.a: $(ASMLIB) $(COBJFILES) $(PLATSTAMP)
	cp $(ASMLIB) $@
	$(AR) -rs $@ $(COBJFILES)

polymult.a: $(POLYOBJFILES) $(PLATSTAMP)
	rm -f $@
	$(AR) -rs $@ $(POLYOBJFILES)

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/polymult_sse2.o: polymult.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -DSSE2 -msse2 -c polymult.c -o $@
$(OBJDIR)/polymult_avx.o: polymult.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -DAVX -mavx -c polymult.c -o $@
$(OBJDIR)/polymult_fma.o: polymult.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -DFMA -mfma -c polymult.c -o $@
$(OBJDIR)/polymult_avx512.o: polymult.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -DAVX512 -mavx512f -c polymult.c -o $@

-include $(wildcard $(OBJDIR)/*.d)

#########################################################################
# The assembly code.  Each file is assembled by uasm into a COFF object (the
# same object the Windows build uses) and converted by objconv into ELF or
# Mach-O.

ifneq ($(wildcard mult.asm),)

asm: $(ASMLIB)

# Any change to a macro file reassembles everything
MACS = $(wildcard *.mac)

# --- 64-bit: linux64/ (ELF, for Linux and FreeBSD) and macosx64/ (Mach-O)

UASMFLAGS64 = -c -nologo -DX86_64 -DLINUX64 -win64 -safeseh
# Objects named NAME_arch are assembled from NAME.asm with ARCH set to match
SUFFIXOBJS64 = xmult2a_core xmult2a_k8 xmult3a_core xmult3a_k8 \
	       hg_blend hg_core hg_p4 hg_k8 hg_k10 \
	       r4_core r4_p4 r4_k8 r4_k10 \
	       r4delay_core r4delay_p4 r4delay_k8 r4delay_k10 \
	       r4delay3_core r4delay3_p4 r4delay3_k8 r4delay3_k10 \
	       r4delay5_core r4delay5_p4 r4delay5_k8 r4delay5_k10 \
	       r4dwpn_core r4dwpn_p4 r4dwpn_k8 r4dwpn_k10 \
	       r4dwpn3_core r4dwpn3_p4 r4dwpn3_k8 r4dwpn3_k10 \
	       r4dwpn5_core r4dwpn5_p4 r4dwpn5_k8 r4dwpn5_k10 \
	       ymult1ax_core ymult1ax_fma3 yr4_core yr4_fma3 ymult3a_core ymult3a_fma3 \
	       yr4dwpn_core yr4dwpn_fma3 yr4dwpn3_core yr4dwpn3_fma3 \
	       yr4dwpn4_core yr4dwpn4_fma3 yr4dwpn5_core yr4dwpn5_fma3 \
	       zmult1ax_skx zr4_skx zmult3a_skx zr4dwpn_skx zr4dwpn2_skx
# Objects assembled from the same-named .asm with ARCH=BLEND
BLENDOBJS64 = xmult1ax xmult2 xmult2ax xmult3 xmult3ax ymult3 ymult3ax zmult3 zmult3ax
# Objects that take no ARCH
PLAINOBJS64 = cpuidhlp hyperhlp mult gianthlp
OBJS64 = $(SUFFIXOBJS64) $(BLENDOBJS64) $(PLAINOBJS64) timeit
ELFOBJS64 = $(addprefix linux64/,$(addsuffix .o,$(OBJS64)))
# The macOS library leaves out the K8 and K10 objects, and uses a mult.obj
# assembled with APPLE defined
MACHOOBJS = $(addprefix macosx64/,$(addsuffix .o,$(filter-out %_k8 %_k10,$(OBJS64))))

linux64/gwnum.a: $(ELFOBJS64)
	rm -f $@
	ar rcs $@ $(ELFOBJS64)

linux64/%.o: linux64/%.obj
	$(OBJCONV) -felf64 $< $@

linux64/%_blend.obj: %.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -DARCH=BLEND -Fo$@ $<
linux64/%_core.obj: %.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -DARCH=CORE -Fo$@ $<
linux64/%_p4.obj: %.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -DARCH=P4 -Fo$@ $<
linux64/%_k8.obj: %.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -DARCH=K8 -Fo$@ $<
linux64/%_k10.obj: %.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -DARCH=K10 -Fo$@ $<
linux64/%_fma3.obj: %.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -DARCH=FMA3 -Fo$@ $<
linux64/%_skx.obj: %.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -DARCH=SKX -Fo$@ $<
$(addprefix linux64/,$(addsuffix .obj,$(BLENDOBJS64))): linux64/%.obj: %.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -DARCH=BLEND -Fo$@ $<
$(addprefix linux64/,$(addsuffix .obj,$(PLAINOBJS64))): linux64/%.obj: %.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -Fo$@ $<
linux64/timeit.obj: timeit.asm $(MACS) | linux64
	$(UASM) $(UASMFLAGS64) -DARCH=FMA3 -Fo$@ $<

macosx64/gwnum.a: $(MACHOOBJS)
	rm -f $@
	ar rcs $@ $(MACHOOBJS)

# -nu+ adds the leading underscore macOS symbol names need
macosx64/mult.o: macosx64/mult.obj
	$(OBJCONV) -fmacho64 -wd1050 -nu+ $< $@
macosx64/%.o: linux64/%.obj | macosx64
	$(OBJCONV) -fmacho64 -wd1050 -nu+ $< $@
macosx64/mult.obj: mult.asm $(MACS) | macosx64
	$(UASM) $(UASMFLAGS64) -DAPPLE -Fo$@ $<

# --- 32-bit: linux/ (ELF)

UASMFLAGS32 = -c -nologo -coff -safeseh
SUFFIXOBJS32 = xmult2a_core xmult2a_k8 xmult3a_core xmult3a_k8 \
	       hg_blend hg_core hg_p4 hg_p4tp hg_k8 hg_k10 \
	       r4_core r4_p4 r4_p4tp r4_k8 r4_k10 \
	       r4delay_p4tp r4delay5_p4tp \
	       r4dwpn_core r4dwpn_p4 r4dwpn_p4tp r4dwpn_k8 r4dwpn_k10 \
	       r4dwpn3_core r4dwpn3_p4 r4dwpn3_p4tp r4dwpn3_k8 r4dwpn3_k10 \
	       r4dwpn5_core r4dwpn5_p4 r4dwpn5_p4tp r4dwpn5_k8 r4dwpn5_k10 \
	       ymult1ax_core yr4_core ymult3a_core \
	       yr4dwpn_core yr4dwpn3_core yr4dwpn4_core yr4dwpn5_core
BLENDOBJS32 = xmult1ax xmult2 xmult2ax xmult3 xmult3ax ymult3 ymult3ax
# Objects assembled from the same-named .asm with PPRO defined
PPROOBJS32 = mult1 mult2 mult3 mult3a mult4
# NAMEp objects are assembled from NAME.asm with PPRO and PFETCH defined
PFETCHOBJS32 = mult2p mult3p mult3ap mult4p
# Objects that take no defines
PLAINOBJS32 = mult mult1aux mult2a mult2aux gianthlp cpuidhlp hyperhlp
OBJS32 = $(PLAINOBJS32) $(PPROOBJS32) $(PFETCHOBJS32) $(BLENDOBJS32) $(SUFFIXOBJS32) timeit
ELFOBJS32 = $(addprefix linux/,$(addsuffix .o,$(OBJS32)))

linux/gwnum.a: $(ELFOBJS32)
	rm -f $@
	ar rcs $@ $(ELFOBJS32)

linux/%.o: linux/%.obj
	$(OBJCONV) -felf32 $< $@

linux/%_blend.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DARCH=BLEND -Fo$@ $<
linux/%_core.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DARCH=CORE -Fo$@ $<
linux/%_p4.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DARCH=P4 -Fo$@ $<
linux/%_p4tp.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DARCH=P4TP -Fo$@ $<
linux/%_k8.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DARCH=K8 -Fo$@ $<
linux/%_k10.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DARCH=K10 -Fo$@ $<
$(addprefix linux/,$(addsuffix .obj,$(BLENDOBJS32))): linux/%.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DARCH=BLEND -Fo$@ $<
$(addprefix linux/,$(addsuffix .obj,$(PPROOBJS32))): linux/%.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DPPRO -Fo$@ $<
$(addprefix linux/,$(addsuffix .obj,$(PFETCHOBJS32))): linux/%p.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DPPRO -DPFETCH -Fo$@ $<
$(addprefix linux/,$(addsuffix .obj,$(PLAINOBJS32))): linux/%.obj: %.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -Fo$@ $<
linux/timeit.obj: timeit.asm $(MACS) | linux
	$(UASM) $(UASMFLAGS32) -DARCH=CORE -Fo$@ $<

# Creating an output directory is the first step of assembling in a fresh
# clone, so check for the tools here
linux64 linux macosx64:
	@command -v $(UASM) > /dev/null || { echo "Assembling gwnum needs uasm, which is not on PATH; see BUILDING.md (or set UASM=/path/to/uasm)" >&2; exit 1; }
	@command -v $(OBJCONV) > /dev/null || { echo "Assembling gwnum needs objconv, which is not on PATH; see BUILDING.md (or set OBJCONV=/path/to/objconv)" >&2; exit 1; }
	mkdir -p $@

asmclean:
	rm -rf linux64 linux macosx64

else

asm:
	@echo "This copy of gwnum comes with its assembly code already assembled (in $(ASMLIB))"

$(ASMLIB):
	@echo "$(ASMLIB) is missing, and this copy of gwnum has no assembly source to build it from" >&2; exit 1

asmclean:

endif

#########################################################################

clean:
	rm -rf build gwnum.a polymult.a

distclean: clean asmclean

.PHONY: all asm clean asmclean distclean
.SECONDARY:
