setlocal

REM Set up the Visual Studio command line environment (any edition, 2019 or later).
REM Skipped if already running in a x86 Developer Command Prompt.

if "%VSCMD_ARG_TGT_ARCH%"=="x86" goto vsready
set VSDIR=
for /f "usebackq delims=" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set VSDIR=%%i
if not defined VSDIR (echo Visual Studio with the C++ tools was not found & exit /b 1)
call "%VSDIR%\VC\Auxiliary\Build\vcvarsall.bat" x86
:vsready

REM Without the private security files (official builds only), create empty
REM placeholders.  Otherwise #include "security.h" finds the Windows SDK's
REM security.h instead.

if not exist ..\security.h type nul > ..\security.h
if not exist ..\secure5.c type nul > ..\secure5.c

nmake /f compile
endlocal
