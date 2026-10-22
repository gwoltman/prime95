# Included by gwnum/Makefile and linux/Makefile: works out which platform
# this is and sets PLATFORM to linux64, linux32, freebsd64 or macos64.
#
# On 64-bit Linux, "make BITS=32" builds for 32-bit Linux instead (gcc then
# needs its 32-bit support, for example Debian's gcc-multilib package).

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S),Linux)
  ifneq ($(filter i386 i486 i586 i686,$(UNAME_M)),)
    NATIVE_BITS := 32
  else ifeq ($(UNAME_M),x86_64)
    NATIVE_BITS := 64
  else
    $(error gwnum needs an x86 processor; this machine is $(UNAME_M))
  endif
  BITS ?= $(NATIVE_BITS)
  ifeq ($(BITS),64)
    ifneq ($(NATIVE_BITS),64)
      $(error BITS=64 needs a 64-bit Linux system)
    endif
    PLATFORM := linux64
  else ifeq ($(BITS),32)
    PLATFORM := linux32
    # -m32 only when building 32-bit on a 64-bit system
    ifeq ($(NATIVE_BITS),64)
      M32 := -m32
    endif
  else
    $(error BITS must be 32 or 64)
  endif
else ifeq ($(UNAME_S),FreeBSD)
  ifneq ($(UNAME_M),amd64)
    $(error only 64-bit FreeBSD (amd64) is supported; this machine is $(UNAME_M))
  endif
  PLATFORM := freebsd64
  # Where packages are installed.  The ports system sets this.
  LOCALBASE ?= /usr/local
else ifeq ($(UNAME_S),Darwin)
  ifneq ($(UNAME_M),x86_64)
    $(error only Intel Macs are supported; this machine is $(UNAME_M))
  endif
  PLATFORM := macos64
  # Targeting releases older than 10.9 (Mavericks) fails: they have no C++ atomics
  export MACOSX_DEPLOYMENT_TARGET = 10.9
else
  $(error no build is defined for $(UNAME_S); see BUILDING.md)
endif

# Object files go in build/$(PLATFORM).  build/platform names the platform
# last built, so that switching platforms (BITS=32, say) rebuilds the
# libraries and programs, not just the objects.
OBJDIR := build/$(PLATFORM)
PLATSTAMP := build/platform
$(shell mkdir -p $(OBJDIR) && { [ "`cat $(PLATSTAMP) 2>/dev/null`" = "$(PLATFORM)" ] || echo $(PLATFORM) > $(PLATSTAMP); })

# Compiler options that also write a .d file listing the headers each object
# depends on, so that changing a header rebuilds what uses it
DEPFLAGS = -MMD -MP

# uasm reads extra command line arguments from an environment variable named
# UASM, so never pass this variable to the commands make runs
UASM = uasm
OBJCONV = objconv
unexport UASM
