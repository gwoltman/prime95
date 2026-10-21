setlocal

REM Set up the Visual Studio command line environment (any edition, 2019 or later).
REM Skipped if already running in a x86 Developer Command Prompt.

if "%VSCMD_ARG_TGT_ARCH%"=="x86" goto vsready
set VSDIR=
for /f "usebackq delims=" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set VSDIR=%%i
if not defined VSDIR (echo Visual Studio with the C++ tools was not found & exit /b 1)
call "%VSDIR%\VC\Auxiliary\Build\vcvarsall.bat" x86
:vsready

REM Directory containing uasm.  Set UASM_DIR to use a different version.

if not defined UASM_DIR set UASM_DIR=c:\uasm_248
@set path=%UASM_DIR%;%path%

REM Create the output directories the makefile writes to

for %%d in (release\polymult debug\polymult linux macosx) do if not exist %%d mkdir %%d

nmake /f compile
endlocal
