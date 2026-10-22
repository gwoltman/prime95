# Building Prime95 and gwnum

This covers Windows with Microsoft Visual C++, 64-bit and 32-bit Linux,
64-bit FreeBSD, and mprime on Intel macOS.

## Windows

### Prerequisites

| What | Default location | Override |
|---|---|---|
| Visual Studio 2019 or later (any edition) with the *Desktop development with C++* workload and the *C++ MFC* component for x86 and x64 | found automatically | |
| UASM 2.58 with the fix from [gwoltman/UASM-2.58-gw1](https://github.com/gwoltman/UASM-2.58-gw1) (see below), as `uasm64.exe` | `c:\uasm_258gw1` | `UASM_DIR` |
| [objconv](https://www.agner.org/optimize/#objconv) | `c:\objconv\objconv.exe` | `OBJCONV` (full path, without `.exe`) |
| GMP 6.2.0 built as `libgmp-gw1` (`gmp.h`, `libgmp-gw1.dll`, `libgmp-gw1.lib`), from the [prime95-windows-deps release](https://github.com/gwoltman/prime95-windows-deps/releases/tag/v1) | `\gmp64`, `c:\gmp32` | `P95_GMP64_DIR`, `P95_GMP32_DIR` |
| [hwloc](https://www.open-mpi.org/projects/hwloc/) 2.8.0 Windows builds | `\hwloc-win64-build-2.8.0`, `c:\hwloc-win32-build-2.8.0` | `P95_HWLOC64_DIR`, `P95_HWLOC32_DIR` |
| cURL 8.5.0 mingw builds (curl-for-win), from the [prime95-windows-deps release](https://github.com/gwoltman/prime95-windows-deps/releases/tag/v1) | `\curl-8.5.0-win64-mingw`, `\curl-8.5.0-win32-mingw` | `P95_CURL64_DIR`, `P95_CURL32_DIR` |
| [Boost](https://www.boost.org/) 1.76.0 headers | `\boost_1_76_0` | `P95_BOOST_DIR` |

A location starting with `\` is at the root of the drive holding the source.

SQLite 3.42.0, cJSON, QD and the other small libraries prime95 uses are
included in this repository, in `third_party/`.

### Getting UASM

Unmodified UASM 2.57 and 2.58 fail on gwnum's `zr4.asm` and `timeit.asm`, and
2.52 and older mis-assemble part of `factor64.asm`.  Use the fixed version:
download `uasm-2.58-gw1-windows-x64.zip` from the
[v2.58-gw1 release](https://github.com/gwoltman/UASM-2.58-gw1/releases/tag/v2.58-gw1)
and copy `uasm64.exe` to `c:\uasm_258gw1`.  Release builds for Linux and
FreeBSD are there too.

Or build it from a Developer Command Prompt and copy it into place:

```
git clone -b v2.58-gw1 https://github.com/gwoltman/UASM-2.58-gw1.git
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

## Linux, FreeBSD and macOS

mprime builds on 64-bit and 32-bit x86 Linux, 64-bit FreeBSD and Intel Macs
(command-line mprime only, not the Prime95 GUI).  Once the prerequisites
below are installed, the build is the same everywhere, in the top-level
directory:

```
make
```

This assembles gwnum's and factor's assembly code, builds gwnum
(`gwnum/gwnum.a` and `gwnum/polymult.a`), and then builds `linux/mprime`.
`make` works out the platform itself.  Running it again rebuilds only what
changed, including after a header change.  `make` in `gwnum/` builds just
gwnum, and `make` in `linux/` does the same as at the top level.

Options, given on the `make` command line:

| Option | Effect |
|---|---|
| `BOOST=/path/to/boost` | Where the Boost headers are, if not in `../../boost_1_76_0` (next to this repository) or where the compiler looks already |
| `UASM=/path/to/uasm`, `OBJCONV=/path/to/objconv` | The assembler tools, if they are not on `PATH` |
| `BITS=32` | Build 32-bit Linux mprime on 64-bit Linux (gcc then needs its 32-bit support, for example the `gcc-multilib` package) |
| `RELEASE=1` | Link the way official releases do (see below) |
| `-j4` | Run 4 compiles at a time |

`make clean` deletes what was built, except the assembled code;
`make distclean` deletes that too.

The assembly code is assembled exactly as on Windows: uasm writes a COFF
object, and objconv converts it to ELF or Mach-O.  uasm's own ELF and Mach-O
output is not used.  FreeBSD uses the same assembled code as 64-bit Linux.

By default mprime links with the shared hwloc and curl libraries, so it
builds with the libraries a distribution provides.  `RELEASE=1` links them
statically, as official releases do, so that mprime runs on more systems; that
needs static builds of them (the notes at the top of `linux/Makefile` give
the configure options used).  It also prints a link map.

### Prerequisites on every platform

- UASM 2.58 with a fix that gwnum's AVX-512 code needs (unmodified UASM 2.57
  and 2.58 fail on `zr4.asm` and `timeit.asm`).  Download
  `uasm-2.58-gw1-linux-x86_64.tar.gz`, `-linux-i686` or
  `-freebsd-x86_64.tar.gz` from the
  [v2.58-gw1 release](https://github.com/gwoltman/UASM-2.58-gw1/releases/tag/v2.58-gw1),
  or build it from source:

  ```
  git clone -b v2.58-gw1 https://github.com/gwoltman/UASM-2.58-gw1.git
  cd UASM-2.58-gw1
  make -f Makefile-Linux-GCC-64.mak
  ```

  This builds `GccUnixR/uasm`.  With an older gcc that defaults to C89,
  such as gcc 4.8 on 32-bit CentOS 7, add C99 mode:
  `make -f Makefile-Linux-GCC-64.mak extra_c_flags="-DNDEBUG -O2 -funsigned-char -Wwrite-strings -std=gnu99"`.
  On FreeBSD use `gmake`; for macOS see below.
- [objconv](https://www.agner.org/optimize/#objconv), built from its
  `source.zip`: `g++ -O2 -o objconv *.cpp` (`clang++` on macOS).
- GMP, hwloc and libcurl, with their development headers.
- The Boost headers (1.76.0 is what official builds use).  A distribution's
  Boost package (such as Debian's `libboost-dev`) works too.
- GNU make.

The gwnum package (see [The gwnum package](#the-gwnum-package)) comes with
the assembly code already assembled, so building gwnum from it needs neither
UASM nor objconv.

### Linux

- gcc and g++ new enough for AVX-512 (`-mavx512f`), which gcc 8 is.  On CentOS
  7, use devtoolset-8: `source /opt/rh/devtoolset-8/enable`.  32-bit builds
  have no AVX-512 code, and 32-bit CentOS 7's gcc 4.8 is enough.
- On Debian or Ubuntu, the libraries are
  `apt install libgmp-dev libhwloc-dev libcurl4-openssl-dev libboost-dev`.

Official Linux builds are made on CentOS 7, so that mprime runs on older Linux
distributions too (an executable linked against a newer glibc will not), with
GMP, hwloc and libcurl built from source in `/usr/local`.

### FreeBSD

- gcc and g++ from ports or packages (`lang/gcc`): mprime is built with gcc
  rather than the system clang, because a clang build of mprime crashes.
- GNU make (`pkg install gmake`).  The makefiles need GNU make; a small
  `BSDmakefile` beside each makes FreeBSD's own `make` run `gmake`, so plain
  `make` still works.
- The GMP, curl and hwloc packages, in `/usr/local` (set `LOCALBASE` for
  another location).

Official FreeBSD builds are made on FreeBSD 12.

### macOS (Intel)

- The Xcode Command Line Tools (`xcode-select --install`): clang, which
  answers to `gcc`/`g++`, and GNU make 3.81.
- GMP, hwloc and libcurl in `/usr/local` (Homebrew's work too).  Official
  builds link GMP and hwloc statically from `/usr/local/lib/libgmp.a` and
  `libhwloc.a`.
- UASM built from source.  Its macOS makefile compiles in C89 mode, which
  current UASM no longer supports, so pass C99 instead:

  ```
  git clone -b v2.58-gw1 https://github.com/gwoltman/UASM-2.58-gw1.git
  cd UASM-2.58-gw1
  make -f Makefile-OSX-Clang-64.mak extra_c_flags="-DNDEBUG -O2 -std=gnu99 -funsigned-char -fwritable-strings"
  ```

mprime runs on macOS 10.9 and later.  This build is tested on macOS 12.1
(Intel) with Apple clang 11 from the Command Line Tools.

## Official Windows releases

`release\windows.ps1` builds the four official Windows zips
(`p95v<version>.win64.zip`, `.win32.zip`, `.win64.service.zip` and
`.win32.service.zip`) from a fresh clone of a commit or tag, with the private
security files from the `prime95-security` repository (by default a
`prime95-security` checkout next to this one):

```
powershell -ExecutionPolicy Bypass -File release\windows.ps1 -Ref <tag or commit> -OutDir <where the zips go>
```

It refuses to package an executable that is still an "Untrusted Prime95"
build, writes a `SHA256SUMS`, and deletes the security files from its work
directory when done.

## Official Linux, FreeBSD and macOS releases

`release/unix.sh` builds `p95v<version>.linux64.tar.gz`, `.linux32.tar.gz`,
`.FreeBSD64.tar.gz` or `.MacOSX.tar.gz`, depending on the machine it runs
on.  It builds the tree it is in, which must be fresh, adds the security
files (`-s`, by default a `prime95-security` directory next to the tree),
checks that mprime is a trusted build of the version in `commonc.h`, and
packages mprime with the text files and, except on macOS where GMP is linked
statically, the GMP library it was linked against (`-g`, default
`/usr/local/lib`).  uasm and objconv (and gmake on FreeBSD) must be on
`PATH`; `-b` gives the Boost directory (default `~/boost_1_76_0`).

From Windows, `release\vm.ps1` does all of this on a build machine over SSH:
it exports the commit or tag, copies it and the security files there, runs
`unix.sh`, fetches the tarball into the same output directory that
`release\windows.ps1` uses, adds it to `SHA256SUMS`, and deletes the security
files from the build machine:

```
powershell -ExecutionPolicy Bypass -File release\vm.ps1 -Target george@127.0.0.1 -Port 2222 -Ref <tag or commit>
powershell -ExecutionPolicy Bypass -File release\vm.ps1 -Target george@127.0.0.1 -Port 2223 -Ref <tag or commit>
powershell -ExecutionPolicy Bypass -File release\vm.ps1 -Target george@127.0.0.1 -Port 2224 -Ref <tag or commit>
powershell -ExecutionPolicy Bypass -File release\vm.ps1 -Target george@192.168.0.39 -Ref <tag or commit>
```

(the 64-bit and 32-bit CentOS VMs, the FreeBSD VM and the Mac).  Official
Linux builds are made on CentOS 7 (glibc 2.17), so they run on old
distributions too, with GMP 6.2 built from source in `/usr/local/lib`;
FreeBSD builds on FreeBSD 12, and Mac builds on macOS 12 for macOS 10.9 and
later.

## Making a release

Set `VERSION` and `BUILD_NUM` in `commonc.h`, commit, and tag the commit
`v<version>b<build>`, for example `v31.7b2` for version 31.7 build 2.  The
release scripts refuse a tag of that form that does not match `commonc.h`.
Then build everything at once:

```
powershell -ExecutionPolicy Bypass -File release\all.ps1 -Ref v31.7b2
```

`release\all.ps1` runs `release\windows.ps1` and a `release\vm.ps1` for each
Unix build machine at the same time, starting the VirtualBox VMs that are
powered off (and shutting them down afterwards, unless `-KeepVMs`).  It
leaves the eight release files and one `SHA256SUMS` in
`prime95-release-<version>` next to this repository, with the build logs in
`%TEMP%\prime95-release-logs`.  The build machines are listed at the top of
the script; `-Only windows,linux64` builds just some of them.

## The gwnum package

`gwnum\package.ps1` packages gwnum for programs that build on it, so they
need neither UASM nor objconv: gwnum's C sources, headers, Unix makefiles and
`tutorial.txt`, the complete Windows libraries, the pre-assembled
`linux64/gwnum.a`, `linux/gwnum.a` and `macosx64/gwnum.a`, and the QD and
SQLite sources gwnum's C code needs (`third_party/qd` and
`third_party/sqlite`).  On Linux, FreeBSD and macOS, users run `make` in its
gwnum directory, which compiles gwnum's C code onto the pre-assembled code;
on Windows they link the libraries directly.  Run it after `gwnum\c64.bat` and `gwnum\c.bat`:

```
powershell -ExecutionPolicy Bypass -File gwnum\package.ps1
```

It writes `p95v<version>b<build>.gwnum.zip`, named from `VERSION` and
`BUILD_NUM` in `commonc.h`.  CI builds it on every run too, as a downloadable
artifact.

## Continuous integration

`.github/workflows/build.yml` builds on GitHub on every push and pull request:
everything on 64-bit Linux, followed by a PRP test of M86243, and everything
on Windows, including `prime95.exe` in all four configurations, the two
Windows services and the gwnum package.  It downloads
UASM from the v2.58-gw1 release, objconv from agner.org, GMP (`libgmp-gw1`)
and cURL 8.5.0 from the prime95-windows-deps release, and hwloc and Boost from
their own sites, checking each against a known SHA-256 hash.

agner.org serves only the latest objconv.  Version 2.57 converts every object
exactly as 2.51 did, except that in 64-bit Mach-O output it writes a different
value into the header's unused `reserved` field.
