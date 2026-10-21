# Building Prime95 and gwnum

This covers Windows with Microsoft Visual C++, and 64-bit and 32-bit Linux.
FreeBSD and macOS instructions will be added as those builds are brought up
to date.

## Windows

### Prerequisites

| What | Default location | Override |
|---|---|---|
| Visual Studio 2019 or later (any edition) with the *Desktop development with C++* workload and the *C++ MFC* component for x86 and x64 | found automatically | |
| UASM 2.58 with the fix from [gwoltman/UASM-2.58-gw1](https://github.com/gwoltman/UASM-2.58-gw1) (see below), as `uasm64.exe` | `c:\uasm_258gw1` | `UASM_DIR` |
| [objconv](https://www.agner.org/optimize/#objconv) | `c:\objconv\objconv.exe` | `OBJCONV` (full path, without `.exe`) |
| GMP, specially built as `libgmp-gw1` (`gmp.h`, `libgmp-gw1.lib`) | `\gmp64`, `c:\gmp32` | `P95_GMP64_DIR`, `P95_GMP32_DIR` |
| [hwloc](https://www.open-mpi.org/projects/hwloc/) 2.8.0 Windows builds | `\hwloc-win64-build-2.8.0`, `c:\hwloc-win32-build-2.8.0` | `P95_HWLOC64_DIR`, `P95_HWLOC32_DIR` |
| [cURL](https://curl.se/windows/) 8.5.0 mingw builds | `\curl-8.5.0-win64-mingw`, `\curl-8.5.0-win32-mingw` | `P95_CURL64_DIR`, `P95_CURL32_DIR` |
| [Boost](https://www.boost.org/) 1.76.0 headers | `\boost_1_76_0` | `P95_BOOST_DIR` |

A location starting with `\` is at the root of the drive holding the source.

SQLite 3.42.0, cJSON, QD and the other small libraries prime95 uses are
included in this repository.

### Building UASM

Unmodified UASM 2.57 and 2.58 fail on gwnum's `zr4.asm` and `timeit.asm`, and
2.52 and older mis-assemble part of `factor64.asm`.  Build the fixed version
from a Developer Command Prompt and copy it into place:

```
git clone -b fix-stale-evex-decorators https://github.com/gwoltman/UASM-2.58-gw1.git
cd UASM-2.58-gw1
msbuild UASM.vcxproj -p:Configuration=Release -p:Platform=x64
mkdir c:\uasm_258gw1
copy x64\Release\UASM.exe c:\uasm_258gw1\uasm64.exe
```

The same 64-bit `uasm64.exe` assembles both the 64-bit and 32-bit code.

### Using other locations

The `UASM_DIR` and `OBJCONV` environment variables are read by the batch files
below.  The `P95_*` locations are set in `prime95\deps.props`; to change them,
either set environment variables with the same names before starting Visual
Studio or MSBuild, or create `prime95\deps.user.props` (it is not committed):

```xml
<Project xmlns="http://schemas.microsoft.com/developer/msbuild/2003">
  <PropertyGroup>
    <P95_GMP64_DIR>D:\libs\gmp64</P95_GMP64_DIR>
  </PropertyGroup>
</Project>
```

### Build steps

Run these from a normal Command Prompt.  Each batch file sets up the Visual
Studio environment itself (or uses the current one, if you are already in a
Developer Command Prompt for the same architecture).

1. Build gwnum.  `c64.bat` builds the 64-bit libraries, `c.bat` the 32-bit ones.
   Each builds both release and debug libraries, plus the gwnum archives the
   Linux and macOS builds link with.

   ```
   cd gwnum
   c64.bat
   c.bat
   ```

2. Build prime95's factoring code, again `c64.bat` for 64-bit and `c.bat` for
   32-bit.

   ```
   cd ..\prime95
   c64.bat
   c.bat
   ```

3. Build `prime95\prime95.sln` in Visual Studio, or from a Developer Command
   Prompt:

   ```
   msbuild prime95.sln -p:Configuration=Release -p:Platform=x64
   ```

   The executable is written to `prime95\x64\Release\prime95.exe` (64-bit) or
   `prime95\Release\prime95.exe` (32-bit).  To run it, put the GMP, hwloc and
   cURL DLLs (for 64-bit: `libgmp-gw1.dll`, `libhwloc-15.dll`,
   `libcurl-x64.dll`) next to it.

### Official builds

Builds made from this repository identify themselves to PrimeNet as
"Untrusted Prime95".  Official builds also need the private `security.h`,
`security.c` and `secure5.c`, copied into the top-level directory.  They are
listed in `.gitignore` and must never be committed.  Without them, the build
creates empty placeholders.  The same applies to the Linux build below.

## Linux (64-bit)

Official Linux builds are made on CentOS 7, so that mprime runs on older Linux
distributions too (an executable linked against a newer glibc will not).

### Prerequisites

- gcc and g++ new enough for AVX-512 (`-mavx512f`).  On CentOS 7, use
  devtoolset-8: `source /opt/rh/devtoolset-8/enable`.
- GNU make.
- GMP, hwloc and libcurl, installed in `/usr/local`.  hwloc and libcurl are
  linked statically; `linux/make64` has the configure options used.
- Boost 1.76.0 headers.  By default `linux/make64` looks for them in
  `../../boost_1_76_0`, next to this repository; pass `BOOST=/path/to/boost`
  to use another location.
- UASM, built from <https://github.com/gwoltman/UASM-2.58-gw1> (branch
  `fix-stale-evex-decorators`): UASM 2.58 with a fix that gwnum's AVX-512 code
  needs.  Unmodified UASM 2.57 and 2.58 fail on `zr4.asm` and `timeit.asm`.

  ```
  git clone -b fix-stale-evex-decorators https://github.com/gwoltman/UASM-2.58-gw1.git
  cd UASM-2.58-gw1
  make -f Makefile-Linux-GCC-64.mak
  ```

  This builds `GccUnixR/uasm`.  With an older gcc that defaults to C89,
  such as gcc 4.8 on 32-bit CentOS 7, add C99 mode:
  `make -f Makefile-Linux-GCC-64.mak extra_c_flags="-DNDEBUG -O2 -funsigned-char -Wwrite-strings -std=gnu99"`.
- [objconv](https://www.agner.org/optimize/#objconv), built from its
  `source.zip`: `g++ -O2 -o objconv *.cpp`.

The assembly code is assembled exactly as on Windows: uasm writes a COFF
object, and objconv converts it to ELF.  uasm's own ELF output is not used.

### Build steps

Pass `UASM` and `OBJCONV` if the tools are not on your PATH.

```
cd gwnum
make -f makeasm64 UASM=~/tools/uasm OBJCONV=~/tools/objconv
make -f make64
cd ../prime95
make -f makeasm64 UASM=~/tools/uasm OBJCONV=~/tools/objconv
cd ../linux
make -f make64
```

`gwnum/makeasm64` builds `gwnum/linux64/gwnum.a` (the assembly code) and
`gwnum/make64` adds the C code to make `gwnum/gwnum.a` and
`gwnum/polymult.a`.  `prime95/makeasm64` builds the factoring code,
`prime95/linux64/factor64.o`.  The result is `linux/mprime`.

## Linux (32-bit)

This build is tested on 32-bit CentOS 7, whose gcc 4.8 is new enough: the
32-bit code has no AVX-512.  The prerequisites are otherwise the
same as for 64-bit Linux, and so is UASM (built as a 32-bit program on a
32-bit system).  The steps use the 32-bit makefiles:

```
cd gwnum
make -f makeasm32 UASM=~/tools/uasm OBJCONV=~/tools/objconv
make -f makefile
cd ../prime95
make -f makeasm32 UASM=~/tools/uasm OBJCONV=~/tools/objconv
cd ../linux
make -f makefile
```

`gwnum/makeasm32` builds `gwnum/linux/gwnum.a` from the same 32-bit objects
the Windows build uses, and `prime95/makeasm32` builds
`prime95/linux/factor32.o`.  The result is `linux/mprime`.
