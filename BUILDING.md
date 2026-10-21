# Building Prime95 and gwnum

This covers Windows with Microsoft Visual C++.  Linux, FreeBSD and macOS
instructions will be added as those builds are brought up to date.

## Windows

### Prerequisites

| What | Default location | Override |
|---|---|---|
| Visual Studio 2019 or later (any edition) with the *Desktop development with C++* workload and the *C++ MFC* component for x86 and x64 | found automatically | |
| [UASM](http://www.terraspace.co.uk/uasm.html) 2.48 for gwnum, 2.52 and 2.49 for prime95's 64-bit and 32-bit factoring code | `c:\uasm_248`, `c:\uasm_252`, `c:\uasm_249` | `UASM_DIR` |
| [objconv](https://www.agner.org/optimize/#objconv) | `c:\objconv\objconv.exe` | `OBJCONV` (full path, without `.exe`) |
| GMP, specially built as `libgmp-gw1` (`gmp.h`, `libgmp-gw1.lib`) | `\gmp64`, `c:\gmp32` | `P95_GMP64_DIR`, `P95_GMP32_DIR` |
| [hwloc](https://www.open-mpi.org/projects/hwloc/) 2.8.0 Windows builds | `\hwloc-win64-build-2.8.0`, `c:\hwloc-win32-build-2.8.0` | `P95_HWLOC64_DIR`, `P95_HWLOC32_DIR` |
| [cURL](https://curl.se/windows/) 8.5.0 mingw builds | `\curl-8.5.0-win64-mingw`, `\curl-8.5.0-win32-mingw` | `P95_CURL64_DIR`, `P95_CURL32_DIR` |
| [Boost](https://www.boost.org/) 1.76.0 headers | `\boost_1_76_0` | `P95_BOOST_DIR` |

A location starting with `\` is at the root of the drive holding the source.

SQLite 3.42.0, cJSON, QD and the other small libraries prime95 uses are
included in this repository.

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
creates empty placeholders.
