# Package gwnum for programs that build on it.
#
# The package holds gwnum's C sources, headers, Unix makefiles and tutorial,
# with the assembly code already assembled, so its users need neither uasm nor
# objconv:
#   - the Windows libraries (gwnum64.lib, gwnum.lib, polymult64.lib, ... and
#     their debug versions), complete;
#   - linux64/gwnum.a, linux/gwnum.a and macosx64/gwnum.a, the assembled code
#     that make64, makefile, makebsd64 and makemac start from;
#   - qd and the SQLite amalgamation, which gwnum's C code includes.
#
# Run it after gwnum\c64.bat and gwnum\c.bat have built everything.  It writes
# p95v<version>b<build>.gwnum.zip (for example p95v3107b02.gwnum.zip) to the
# top-level directory, or to -OutDir.

param([string]$OutDir)
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.IO.Compression, System.IO.Compression.FileSystem

$gw = $PSScriptRoot
$root = Split-Path -Parent $gw
if (-not $OutDir) { $OutDir = $root }

# Name the package from VERSION and BUILD_NUM in commonc.h
$h = Get-Content (Join-Path $root 'commonc.h') -Raw
$v = [regex]::Match($h, '#define\s+VERSION\s+"(\d+)\.(\d+)"')
$b = [regex]::Match($h, '#define\s+BUILD_NUM\s+"(\d+)"')
if (-not $v.Success -or -not $b.Success) { throw 'VERSION or BUILD_NUM not found in commonc.h' }
$ver = '{0}{1:D2}b{2:D2}' -f $v.Groups[1].Value, [int]$v.Groups[2].Value, [int]$b.Groups[1].Value
$zipPath = Join-Path $OutDir "p95v$ver.gwnum.zip"

# Committed source files, so no build leftovers get in
$tracked = git -C $root ls-files gwnum qd sqlite-amalgamation-3420000
if ($LASTEXITCODE -ne 0) { throw 'git ls-files failed' }
$wanted = '\.(c|cpp|h|cc)$|^gwnum/(make64|makefile|makebsd64|makemac|makehaiku|makemsys|tutorial\.txt)$'
$files = @($tracked | Where-Object { $_ -match $wanted })

# The assembled code, from the Windows build
$built = 'gwnum/gwnum64.lib', 'gwnum/gwnum64d.lib', 'gwnum/polymult64.lib', 'gwnum/polymult64d.lib',
         'gwnum/gwnum.lib', 'gwnum/gwnumd.lib', 'gwnum/polymult.lib', 'gwnum/polymultd.lib',
         'gwnum/linux64/gwnum.a', 'gwnum/linux/gwnum.a', 'gwnum/macosx64/gwnum.a'
$missing = @($built | Where-Object { -not (Test-Path (Join-Path $root $_)) })
if ($missing) { throw "Not built yet (run gwnum\c64.bat and gwnum\c.bat): $missing" }

$readme = @"
gwnum $($v.Groups[1].Value).$($v.Groups[2].Value) build $($b.Groups[1].Value)
==================

The gwnum library from prime95 (https://github.com/gwoltman/prime95), with
its assembly code already assembled.  See gwnum/tutorial.txt for how to use
gwnum, and license.txt for the terms.

Windows (Microsoft Visual C++)
  Link with gwnum64.lib and polymult64.lib (64-bit) or gwnum.lib and
  polymult.lib (32-bit), from the gwnum directory.  They are built with the
  static C runtime (/MT); gwnum64d.lib, polymult64d.lib, gwnumd.lib and
  polymultd.lib are debug versions built with /MTd.

Linux, FreeBSD and macOS
  Compile gwnum's C code, which adds it to the assembled code, by running
  one of these in the gwnum directory:
      make -f make64        64-bit Linux
      make -f makefile      32-bit Linux
      make -f makebsd64     64-bit FreeBSD
      make -f makemac       macOS (Intel)
  This produces gwnum.a and polymult.a; link with both.  On macOS, also link
  with -Wl,-no_pie: the assembly code uses absolute addresses.

Building gwnum from scratch, including the assembly code, needs the full
prime95 source and the tools described in its BUILDING.md.
"@

New-Item -ItemType Directory -Force $OutDir | Out-Null
if (Test-Path $zipPath) { Remove-Item $zipPath }
$zip = [System.IO.Compression.ZipFile]::Open($zipPath, 'Create')
try {
    foreach ($f in $files + $built) {
        [void][System.IO.Compression.ZipFileExtensions]::CreateEntryFromFile($zip, (Join-Path $root $f), $f, 'Optimal')
    }
    [void][System.IO.Compression.ZipFileExtensions]::CreateEntryFromFile($zip, (Join-Path $root 'prime95\license.txt'), 'license.txt', 'Optimal')
    $e = $zip.CreateEntry('README.txt')
    $w = New-Object System.IO.StreamWriter($e.Open())
    $w.Write(($readme -replace "`r?`n", "`r`n"))
    $w.Close()
} finally {
    $zip.Dispose()
}
Write-Host ("{0}: {1} files, {2:N1} MB" -f $zipPath, ($files.Count + $built.Count + 2), ((Get-Item $zipPath).Length / 1MB))
