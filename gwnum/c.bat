setlocal
if exist "c:\program files (x86)" goto x64
REM if exist "c:\program files\microsoft visual studio 10.0" goto x32_10
if exist "c:\program files\microsoft visual studio 8" goto x32_8
call "c:\program files\microsoft visual studio 9.0\vc\bin\vcvars32"
goto makeit
:x32_10
call "c:\program files\microsoft visual studio 10.0\vc\bin\vcvars32"
goto makeit
:x32_8
call "c:\program files\microsoft visual studio 8\vc\bin\vcvars32"
goto makeit
:x64
REM if exist "c:\program files (x86)\microsoft visual studio 10.0" goto x64_10
if exist "c:\program files (x86)\microsoft visual studio 8" goto x64_8
call "c:\program files (x86)\microsoft visual studio 9.0\vc\bin\vcvars32"
goto makeit
:x64_10
call "c:\program files (x86)\microsoft visual studio 10.0\vc\bin\vcvars32"
goto makeit
:x64_8
call "c:\program files (x86)\microsoft visual studio 8\vc\bin\vcvars32"
:makeit
@set path=y:\downloads\hjwasm216_x86;%path%
nmake /f compile
endlocal
