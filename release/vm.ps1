# Build an official Unix release of mprime on another machine over SSH, such
# as a Linux or FreeBSD virtual machine or a Mac, and fetch the tarball.
#
# It exports -Ref (a tag, branch or commit; default HEAD) of this repository,
# copies it and the private security files from -SecurityDir to the build
# machine, runs release/unix.sh there, and copies the tarball back into -OutDir
# (default: prime95-release-<version> next to this repository, where
# release\windows.ps1 puts the Windows zips), adding it to SHA256SUMS.  The
# security files are deleted from the build machine at the end.
#
# The build machine needs what release/unix.sh needs, and SSH key login.
# -Setup is run there first, by sh whatever the login shell is; by default it
# puts ~/tools/bin (uasm, objconv and on FreeBSD gmake) on PATH and enables a
# Red Hat devtoolset compiler if one is installed.  -UnixArgs are passed to
# unix.sh, for example '-b /opt/boost_1_76_0'.
#
#   powershell -ExecutionPolicy Bypass -File release\vm.ps1 -Target george@127.0.0.1 -Port 2222 -Ref v31.7b2

param(
    [Parameter(Mandatory)][string]$Target,
    [int]$Port = 22,
    [string]$Ref = 'HEAD',
    [string]$SecurityDir,
    [string]$OutDir,
    [string]$Setup = 'PATH=$HOME/tools/bin:$PATH; for f in /opt/rh/devtoolset-*/enable; do [ -f "$f" ] && . "$f"; done',
    [string]$UnixArgs = ''
)
$ErrorActionPreference = 'Stop'

$repo = Split-Path -Parent $PSScriptRoot
if (-not $SecurityDir) { $SecurityDir = Join-Path (Split-Path -Parent $repo) 'prime95-security' }
$secFiles = 'security.h', 'security.c', 'secure5.c'
foreach ($f in $secFiles) {
    if (-not (Test-Path (Join-Path $SecurityDir $f))) { throw "$f not found in $SecurityDir (use -SecurityDir)" }
}
$remote = 'prime95-release-work'    # relative to the home directory on the build machine
$ssh = @('-o', 'BatchMode=yes', '-p', $Port, $Target)
$scp = @('-o', 'BatchMode=yes', '-q', '-P', $Port)
function Step($text) { Write-Host "== $text" }
function Remote($cmd) { & ssh @ssh $cmd; if ($LASTEXITCODE -ne 0) { throw "Failed on ${Target}: $cmd" } }

$commit = git -C $repo rev-parse --short "$Ref^{commit}"
if ($LASTEXITCODE -ne 0) { throw "Unknown ref $Ref" }
$h = (git -C $repo show "${Ref}:commonc.h") -join "`n"
$v = [regex]::Match($h, '#define\s+VERSION\s+"(\d+)\.(\d+)"')
$b = [regex]::Match($h, '#define\s+BUILD_NUM\s+"(\d+)"')
if (-not $v.Success -or -not $b.Success) { throw 'VERSION or BUILD_NUM not found in commonc.h' }
$ver = '{0}{1:D2}b{2:D2}' -f $v.Groups[1].Value, [int]$v.Groups[2].Value, [int]$b.Groups[1].Value
if (-not $OutDir) { $OutDir = Join-Path (Split-Path -Parent $repo) "prime95-release-$ver" }
Write-Host "Version $($v.Groups[1].Value).$($v.Groups[2].Value) build $($b.Groups[1].Value) ($ver), commit $commit, on $Target port $Port"

# The source, with Unix line endings whatever this checkout uses
$tar = Join-Path $env:TEMP "prime95-$commit.tar"
$runSh = Join-Path $env:TEMP "prime95-$commit-run.sh"
Step "Exporting $Ref"
git -C $repo -c core.autocrlf=false -c core.eol=lf archive --format=tar --prefix=prime95/ -o $tar $commit
if ($LASTEXITCODE -ne 0) { throw 'git archive failed' }
# What runs on the build machine, as a script so that no remote shell parses it
$run = @(
    "cd `"`$HOME/$remote`" && tar -xf src.tar && rm src.tar || exit 1"
    $Setup
    "exec sh prime95/release/unix.sh -s `"`$HOME/$remote/sec`" -o `"`$HOME/$remote/out`" $UnixArgs"
)
[IO.File]::WriteAllText($runSh, (($run -join "`n") + "`n"))

try {
    Step "Copying the source and security files to $Target"
    Remote "rm -rf $remote && mkdir -p $remote/sec && chmod 700 $remote/sec"
    & scp @scp $tar "${Target}:$remote/src.tar"
    if ($LASTEXITCODE -ne 0) { throw 'scp of the source failed' }
    & scp @scp $runSh "${Target}:$remote/run.sh"
    if ($LASTEXITCODE -ne 0) { throw 'scp of run.sh failed' }
    & scp @scp ($secFiles | ForEach-Object { Join-Path $SecurityDir $_ }) "${Target}:$remote/sec/"
    if ($LASTEXITCODE -ne 0) { throw 'scp of the security files failed' }

    Step "Building on $Target"
    Remote "sh $remote/run.sh"

    Step "Fetching into $OutDir"
    New-Item -ItemType Directory -Force $OutDir | Out-Null
    $remoteSums = @{}
    foreach ($line in @(& ssh @ssh "cat $remote/out/SHA256SUMS")) {
        if ($line -match '^([0-9a-f]{64})  (\S+)$') { $remoteSums[$Matches[2]] = $Matches[1] }
    }
    $names = @($remoteSums.Keys | Sort-Object)
    if (-not $names) { throw "No tarball was built on $Target" }
    foreach ($n in $names) {
        & scp @scp "${Target}:$remote/out/$n" (Join-Path $OutDir $n)
        if ($LASTEXITCODE -ne 0) { throw "scp of $n failed" }
        if ((Get-FileHash (Join-Path $OutDir $n) -Algorithm SHA256).Hash.ToLower() -ne $remoteSums[$n]) { throw "$n was corrupted in transfer" }
        Write-Host "  $n"
    }

    # Add them to SHA256SUMS, replacing older entries of the same name
    $sumsPath = Join-Path $OutDir 'SHA256SUMS'
    $sums = @()
    if (Test-Path $sumsPath) { $sums = @(Get-Content $sumsPath | Where-Object { $_ -and ($_ -split '  ', 2)[1] -notin $names }) }
    $sums += $names | ForEach-Object { '{0}  {1}' -f (Get-FileHash (Join-Path $OutDir $_) -Algorithm SHA256).Hash.ToLower(), $_ }
    $sums = @($sums) | Sort-Object { ($_ -split '  ', 2)[1] }
    [IO.File]::WriteAllText($sumsPath, (($sums -join "`n") + "`n"))
    Write-Host "Done: $ver from commit $commit"
} finally {
    Remove-Item -ErrorAction SilentlyContinue $tar, $runSh
    # unix.sh deletes its copies from the tree; this deletes the originals
    & ssh @ssh "rm -rf $remote/sec $remote/prime95/security.h $remote/prime95/security.c $remote/prime95/secure5.c"
}
