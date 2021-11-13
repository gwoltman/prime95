setlocal
REM if exist "c:\program files (x86)" goto x64
REM if exist "c:\program files\microsoft visual studio 10.0" goto x32_10
REM if exist "c:\program files\microsoft visual studio 8" goto x32_8
REM call "c:\program files\microsoft visual studio 9.0\vc\bin\x86_amd64\vcvarsx86_amd64"
REM goto makeit
REM :x32_10
REM call "c:\program files\microsoft visual studio 10.0\vc\bin\x86_amd64\vcvarsx86_amd64"
REM goto makeit
REM :x32_8
REM call "c:\program files\microsoft visual studio 8\vc\bin\x86_amd64\vcvarsx86_amd64"
REM goto makeit
REM :x64
REM if exist "c:\program files (x86)\microsoft visual studio 10.0" goto x64_10
REM if exist "c:\program files\microsoft visual studio 8" goto x64_8
REM call "c:\program files\microsoft visual studio 9.0\vc\bin\amd64\vcvarsamd64"
REM goto makeit
REM :x64_10
REM call "c:\program files\microsoft visual studio 10.0\vc\bin\amd64\vcvarsamd64"
REM goto makeit
REM :x64_8
REM call "c:\program files\microsoft visual studio 8\vc\bin\amd64\vcvarsamd64"
call "c:\program files (x86)\microsoft visual studio\2019\Professional\VC\Auxiliary\build\vcvarsall" x64
:makeit
@set path=c:\uasm_252;%path%
nmake /f compil64
endlocal
