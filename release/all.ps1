# Build every official release of prime95 and mprime into one directory.
#
# Runs release\windows.ps1 here and release\vm.ps1 for each Unix build
# machine below, all at the same time, and leaves in -OutDir (default:
# prime95-release-<version> next to this repository):
#   p95v<version>.win64.zip, .win32.zip, .win64.service.zip, .win32.service.zip
#   p95v<version>.linux64.tar.gz, .linux32.tar.gz, .FreeBSD64.tar.gz, .MacOSX.tar.gz
#   SHA256SUMS
# Build logs go to -LogDir (default: %TEMP%\prime95-release-logs).
#
# VirtualBox build machines that are powered off are started headless, and
# shut down again at the end unless -KeepVMs is given.  -Only picks some of the
# builds, for example -Only windows,linux64 (the four Windows zips are one build,
# named windows).
#
#   powershell -ExecutionPolicy Bypass -File release\all.ps1 -Ref v31.7b2

param(
    [string]$Ref = 'HEAD',
    [string]$SecurityDir,
    [string]$OutDir,
    [string]$LogDir = (Join-Path $env:TEMP 'prime95-release-logs'),
    [string[]]$Only,
    [switch]$KeepVMs
)
$ErrorActionPreference = 'Stop'

# The Unix build machines: SSH target and port, and the VirtualBox VM if it is one
$machines = [ordered]@{
    linux64   = @{ Target = 'george@127.0.0.1'; Port = 2222; VM = 'MSI CentOS7 64-bit' }
    linux32   = @{ Target = 'george@127.0.0.1'; Port = 2223; VM = 'MSI CentOS7 32-bit' }
    FreeBSD64 = @{ Target = 'george@127.0.0.1'; Port = 2224; VM = 'MSI FreeBSD12 64-bit' }
    MacOSX    = @{ Target = 'george@192.168.0.39'; Port = 22; VM = $null }
}
$vbox = "$env:ProgramFiles\Oracle\VirtualBox\VBoxManage.exe"

$repo = Split-Path -Parent $PSScriptRoot
if (-not $SecurityDir) { $SecurityDir = Join-Path (Split-Path -Parent $repo) 'prime95-security' }
$builds = @('windows') + $machines.Keys
if ($Only) {
    foreach ($o in $Only) { if ($o -notin $builds) { throw "Unknown build $o (builds: $($builds -join ', '))" } }
    $builds = $builds | Where-Object { $_ -in $Only }
}
function Step($text) { Write-Host "== $text" }

# The version, to name the output directory and the expected files
$commit = git -C $repo rev-parse --short "$Ref^{commit}"
if ($LASTEXITCODE -ne 0) { throw "Unknown ref $Ref" }
$h = (git -C $repo show "${Ref}:commonc.h") -join "`n"
$v = [regex]::Match($h, '#define\s+VERSION\s+"(\d+)\.(\d+)"')
$b = [regex]::Match($h, '#define\s+BUILD_NUM\s+"(\d+)"')
if (-not $v.Success -or -not $b.Success) { throw 'VERSION or BUILD_NUM not found in commonc.h' }
$ver = '{0}{1:D2}b{2:D2}' -f $v.Groups[1].Value, [int]$v.Groups[2].Value, [int]$b.Groups[1].Value
if (-not $OutDir) { $OutDir = Join-Path (Split-Path -Parent $repo) "prime95-release-$ver" }
New-Item -ItemType Directory -Force $OutDir, $LogDir | Out-Null
$OutDir = (Resolve-Path $OutDir).Path
Write-Host "Version $ver from commit ${commit}: $($builds -join ', ')"
Write-Host "Output in $OutDir, logs in $LogDir"

$started = @()
$procs = [ordered]@{}
try {
    # Start the build machines that are VirtualBox VMs and powered off
    foreach ($name in $builds | Where-Object { $_ -ne 'windows' -and $machines[$_].VM }) {
        $m = $machines[$name]
        $state = (& $vbox showvminfo $m.VM --machinereadable | Select-String '^VMState=').ToString()
        if ($state -ne 'VMState="running"') {
            Step "Starting VM $($m.VM)"
            & $vbox startvm $m.VM --type headless | Out-Null
            if ($LASTEXITCODE -ne 0) { throw "Could not start VM $($m.VM)" }
            $started += $m.VM
        }
    }
    foreach ($name in $builds | Where-Object { $_ -ne 'windows' }) {
        $m = $machines[$name]
        $ok = $false
        for ($i = 0; $i -lt 36 -and -not $ok; $i++) {
            # (Windows PowerShell turns a native command's redirected stderr
            # into an error that 'Stop' would throw, so this runs with 'Continue')
            $ok = & { $ErrorActionPreference = 'Continue'; & ssh -o BatchMode=yes -o ConnectTimeout=5 -p $m.Port $m.Target true 2>$null; $LASTEXITCODE -eq 0 }
            if (-not $ok) { Start-Sleep 5 }
        }
        if (-not $ok) { throw "Cannot reach $name ($($m.Target) port $($m.Port)) over SSH" }
    }

    # Start every build.  Each Unix build gets its own output directory, so
    # that they do not write SHA256SUMS at the same time.
    foreach ($name in $builds) {
        $log = Join-Path $LogDir "$name.log"
        if ($name -eq 'windows') {
            $script = 'windows.ps1'
            $args_ = @('-Ref', $Ref, '-SecurityDir', $SecurityDir, '-OutDir', (Join-Path $OutDir ".$name"))
        } else {
            $m = $machines[$name]
            $script = 'vm.ps1'
            $args_ = @('-Target', $m.Target, '-Port', $m.Port, '-Ref', $Ref, '-SecurityDir', $SecurityDir, '-OutDir', (Join-Path $OutDir ".$name"))
        }
        Step "Starting the $name build"
        $quoted = $args_ | ForEach-Object { '"{0}"' -f $_ }
        $procs[$name] = Start-Process powershell -PassThru -NoNewWindow -RedirectStandardOutput $log -RedirectStandardError "$log.err" `
            -ArgumentList (@('-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', "`"$(Join-Path $PSScriptRoot $script)`"") + $quoted)
        # .NET only keeps the exit code of a process whose handle was opened while it ran
        [void]$procs[$name].Handle
    }

    # Wait for them all
    $failed = @()
    foreach ($name in $procs.Keys) {
        $p = $procs[$name]
        $p.WaitForExit()
        $log = Join-Path $LogDir "$name.log"
        Get-Content "$log.err" -ErrorAction SilentlyContinue | Add-Content $log
        Remove-Item -ErrorAction SilentlyContinue "$log.err"
        if ($p.ExitCode -eq 0) {
            Step "$name build done"
        } else {
            Step "$name build FAILED (exit $($p.ExitCode)); the end of ${log}:"
            Get-Content $log -Tail 15 | ForEach-Object { "   $_" }
            $failed += $name
        }
    }

    # Gather the results into one directory with one SHA256SUMS
    $kinds = @()
    foreach ($name in $builds) {
        if ($name -eq 'windows') { $kinds += 'win64.zip', 'win32.zip', 'win64.service.zip', 'win32.service.zip' }
        else { $kinds += "$name.tar.gz" }
    }
    foreach ($name in $builds | Where-Object { $_ -notin $failed }) {
        $dir = Join-Path $OutDir ".$name"
        Get-ChildItem $dir -File | Where-Object { $_.Name -ne 'SHA256SUMS' } | Move-Item -Destination $OutDir -Force
        Remove-Item -Recurse -Force $dir
    }
    $files = $kinds | ForEach-Object { "p95v$ver.$_" } | Where-Object { Test-Path (Join-Path $OutDir $_) }
    $sumsPath = Join-Path $OutDir 'SHA256SUMS'
    $sums = @()
    if (Test-Path $sumsPath) { $sums = @(Get-Content $sumsPath | Where-Object { $_ -and ($_ -split '  ', 2)[1] -notin $files }) }
    $sums += $files | ForEach-Object { '{0}  {1}' -f (Get-FileHash (Join-Path $OutDir $_) -Algorithm SHA256).Hash.ToLower(), $_ }
    $sums = @($sums) | Sort-Object { ($_ -split '  ', 2)[1] }
    if ($sums) { [IO.File]::WriteAllText($sumsPath, (($sums -join "`n") + "`n")) }

    Step "Results in $OutDir"
    foreach ($k in $kinds) {
        $f = Join-Path $OutDir "p95v$ver.$k"
        if (Test-Path $f) { Write-Host ("  p95v{0}.{1,-20} {2,6:N1} MB" -f $ver, $k, ((Get-Item $f).Length / 1MB)) }
        else { Write-Host "  p95v$ver.$k  MISSING" }
    }
    if ($failed) { throw "Failed: $($failed -join ', ') (logs in $LogDir)" }
    Write-Host "Done: $ver from commit $commit"
} finally {
    foreach ($p in $procs.Values) { if (-not $p.HasExited) { $p.Kill() } }
    if (-not $KeepVMs) {
        foreach ($vm in $started) {
            Step "Shutting down VM $vm"
            & $vbox controlvm $vm acpipowerbutton | Out-Null
        }
    }
}
