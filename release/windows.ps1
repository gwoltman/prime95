# Build the official Windows release of prime95.
#
# Produces, in -OutDir (default: prime95-release-<version> next to this
# repository):
#   p95v<version>.win64.zip           prime95.exe (64-bit)
#   p95v<version>.win32.zip           prime95.exe (32-bit)
#   p95v<version>.win64.service.zip   ntprime64.exe
#   p95v<version>.win32.service.zip   ntprime.exe
#   SHA256SUMS
# where <version> comes from commonc.h, for example 3107b02.
#
# It clones -Ref (a tag, branch or commit; default HEAD) of this repository
# into a fresh work directory, so only committed source is built, adds the
# private security files from -SecurityDir, builds everything, checks that
# every executable is a trusted build, and packages each one with its DLLs and
# text files.  The security files are deleted from the work directory at the
# end.
#
# Needs everything BUILDING.md lists for Windows.  Library locations come from
# the same P95_* environment variables as prime95\deps.props, with the same
# defaults.
#
#   powershell -ExecutionPolicy Bypass -File release\windows.ps1 -Ref v31.7b2

param(
    [string]$Ref = 'HEAD',
    [string]$SecurityDir,
    [string]$OutDir,
    [string]$WorkDir
)
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.IO.Compression, System.IO.Compression.FileSystem

$repo = Split-Path -Parent $PSScriptRoot
if (-not $SecurityDir) { $SecurityDir = Join-Path (Split-Path -Parent $repo) 'prime95-security' }
if (-not $WorkDir) { $WorkDir = Join-Path $env:TEMP 'prime95-release' }
# Some locked-down environments stop cmd from running batch files in the
# current directory; the build's batch files call each other that way.
$env:NoDefaultCurrentDirectoryInExePath = $null

function Step($text) { Write-Host "== $text" }

# Library locations, as in prime95\deps.props
function Dep($name, $default) { $v = [Environment]::GetEnvironmentVariable($name); if ($v) { $v } else { $default } }
$dirs = @{
    gmp64 = Dep P95_GMP64_DIR '\gmp64';            gmp32 = Dep P95_GMP32_DIR 'c:\gmp32'
    hwloc64 = Dep P95_HWLOC64_DIR '\hwloc-win64-build-2.8.0'; hwloc32 = Dep P95_HWLOC32_DIR 'c:\hwloc-win32-build-2.8.0'
    curl64 = Dep P95_CURL64_DIR '\curl-8.5.0-win64-mingw';   curl32 = Dep P95_CURL32_DIR '\curl-8.5.0-win32-mingw'
}

# Security files
$secFiles = 'security.h', 'security.c', 'secure5.c'
foreach ($f in $secFiles) {
    if (-not (Test-Path (Join-Path $SecurityDir $f))) { throw "$f not found in $SecurityDir (use -SecurityDir)" }
}
if (-not (Select-String -Path (Join-Path $SecurityDir 'security.h') -Pattern 'SECURITY_MODULES_PRESENT' -Quiet)) {
    throw "security.h in $SecurityDir does not define SECURITY_MODULES_PRESENT"
}

# Fresh source
Step "Cloning $Ref into $WorkDir"
if (Test-Path $WorkDir) { Remove-Item -Recurse -Force $WorkDir }
git clone -q $repo $WorkDir
if ($LASTEXITCODE -ne 0) { throw 'git clone failed' }
git -C $WorkDir checkout -q $Ref
if ($LASTEXITCODE -ne 0) { throw "git checkout $Ref failed" }
$commit = git -C $WorkDir rev-parse --short HEAD
foreach ($f in $secFiles) { Copy-Item (Join-Path $SecurityDir $f) $WorkDir }

$h = Get-Content (Join-Path $WorkDir 'commonc.h') -Raw
$v = [regex]::Match($h, '#define\s+VERSION\s+"(\d+)\.(\d+)"')
$b = [regex]::Match($h, '#define\s+BUILD_NUM\s+"(\d+)"')
if (-not $v.Success -or -not $b.Success) { throw 'VERSION or BUILD_NUM not found in commonc.h' }
# A release tag (v31.7b2 for version 31.7 build 2) must name the version it builds
$t = [regex]::Match($Ref, '^v(\d+)\.(\d+)b(\d+)$')
if ($t.Success -and ([int]$t.Groups[1].Value -ne [int]$v.Groups[1].Value -or [int]$t.Groups[2].Value -ne [int]$v.Groups[2].Value -or [int]$t.Groups[3].Value -ne [int]$b.Groups[1].Value)) {
    throw "Tag $Ref does not match commonc.h, which says version $($v.Groups[1].Value).$($v.Groups[2].Value) build $($b.Groups[1].Value)"
}
$ver = '{0}{1:D2}b{2:D2}' -f $v.Groups[1].Value, [int]$v.Groups[2].Value, [int]$b.Groups[1].Value
if (-not $OutDir) { $OutDir = Join-Path (Split-Path -Parent $repo) "prime95-release-$ver" }
Write-Host "Version $($v.Groups[1].Value).$($v.Groups[2].Value) build $($b.Groups[1].Value) ($ver), commit $commit"

try {
    # Build
    foreach ($s in @('gwnum\c64.bat', 'gwnum\c.bat', 'prime95\c64.bat', 'prime95\c.bat')) {
        Step "Building $s"
        Push-Location (Join-Path $WorkDir (Split-Path -Parent $s))
        cmd /c call (Split-Path -Leaf $s)
        $rc = $LASTEXITCODE
        Pop-Location
        if ($rc -ne 0) { throw "$s failed" }
    }
    $msbuild = & "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe" -latest -requires Microsoft.Component.MSBuild -find 'MSBuild\**\Bin\amd64\MSBuild.exe' | Select-Object -First 1
    foreach ($p in 'x64', 'Win32') {
        Step "Building prime95.exe ($p)"
        & $msbuild (Join-Path $WorkDir 'prime95\prime95.sln') -nologo -v:minimal -m -p:Configuration=Release "-p:Platform=$p"
        if ($LASTEXITCODE -ne 0) { throw "prime95.sln Release|$p failed" }
    }
    foreach ($s in 'c64.bat', 'c.bat') {
        Step "Building the Windows service (winnt\$s)"
        Push-Location (Join-Path $WorkDir 'winnt')
        cmd /c call $s
        $rc = $LASTEXITCODE
        Pop-Location
        if ($rc -ne 0) { throw "winnt\$s failed" }
    }

    # The executables, and a check that each one is a trusted build
    $exes = @{
        win64 = 'prime95\x64\Release\prime95.exe'; win32 = 'prime95\Release\prime95.exe'
        'win64.service' = 'winnt\ntprime64.exe'; 'win32.service' = 'winnt\ntprime.exe'
    }
    foreach ($e in $exes.Values) {
        $contents = [Text.Encoding]::ASCII.GetString([IO.File]::ReadAllBytes((Join-Path $WorkDir $e)))
        if ($contents.Contains('Untrusted Prime95')) { throw "$e is an untrusted build: the security files were not used" }
    }

    # Package
    Step "Packaging into $OutDir"
    New-Item -ItemType Directory -Force $OutDir | Out-Null
    $text = 'license.txt', 'readme.txt', 'stress.txt', 'undoc.txt', 'whatsnew.txt' | ForEach-Object { Join-Path $WorkDir "prime95\$_" }
    $dlls = @{
        64 = (Join-Path $dirs.curl64 'bin\libcurl-x64.dll'), (Join-Path $dirs.gmp64 'libgmp-gw1.dll'), (Join-Path $dirs.hwloc64 'bin\libhwloc-15.dll')
        32 = (Join-Path $dirs.curl32 'bin\libcurl.dll'), (Join-Path $dirs.gmp32 'libgmp-gw1.dll'), (Join-Path $dirs.hwloc32 'bin\libhwloc-15.dll')
    }
    $zips = @()
    foreach ($kind in 'win64', 'win32', 'win64.service', 'win32.service') {
        $bits = if ($kind -like 'win64*') { 64 } else { 32 }
        $files = @(Join-Path $WorkDir $exes[$kind]) + $dlls[$bits] + $text
        if ($kind -like '*.service') { $files += Join-Path $WorkDir 'winnt\readme.service.txt' }
        $zipPath = Join-Path $OutDir "p95v$ver.$kind.zip"
        if (Test-Path $zipPath) { Remove-Item $zipPath }
        $zip = [IO.Compression.ZipFile]::Open($zipPath, 'Create')
        try {
            foreach ($f in ($files | Sort-Object { Split-Path -Leaf $_ })) {
                if (-not (Test-Path $f)) { throw "Missing: $f" }
                [void][IO.Compression.ZipFileExtensions]::CreateEntryFromFile($zip, $f, (Split-Path -Leaf $f), 'Optimal')
            }
        } finally { $zip.Dispose() }
        $zips += $zipPath
        Write-Host ("  {0}: {1} files, {2:N1} MB" -f (Split-Path -Leaf $zipPath), $files.Count, ((Get-Item $zipPath).Length / 1MB))
    }
    $sums = $zips | ForEach-Object { '{0}  {1}' -f (Get-FileHash $_ -Algorithm SHA256).Hash.ToLower(), (Split-Path -Leaf $_) }
    [IO.File]::WriteAllText((Join-Path $OutDir 'SHA256SUMS'), (($sums -join "`n") + "`n"))
    Write-Host "Done: $ver from commit $commit"
} finally {
    foreach ($f in $secFiles) { Remove-Item -ErrorAction SilentlyContinue (Join-Path $WorkDir $f) }
}
