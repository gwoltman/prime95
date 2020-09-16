setlocal
REM if exist "c:\program files (x86)" goto x64
REM if exist "c:\program files\microsoft visual studio 10.0" goto x32_10
REM if exist "c:\program files\microsoft visual studio 8" goto x32_8
REM call "c:\program files\microsoft visual studio 9.0\vc\bin\vcvars32"
REM goto makeit
REM :x32_10
REM call "c:\program files\microsoft visual studio 10.0\vc\bin\vcvars32"
REM goto makeit
REM :x32_8
REM call "c:\program files\microsoft visual studio 8\vc\bin\vcvars32"
REM goto makeit
REM :x64
REM if exist "c:\program files (x86)\microsoft visual studio 10.0" goto x64_10
REM if exist "c:\program files (x86)\microsoft visual studio 8" goto x64_8
REM call "c:\program files (x86)\microsoft visual studio 9.0\vc\bin\vcvars32"
REM goto makeit
REM :x64_10
REM call "c:\program files (x86)\microsoft visual studio 10.0\vc\bin\vcvars32"
REM goto makeit
REM :x64_8
call "c:\program files (x86)\microsoft visual studio\2019\Professional\VC\Auxiliary\build\vcvarsall" x86
:makeit
@set path=c:\uasm_249;%path%
nmake /f compile
endlocal
