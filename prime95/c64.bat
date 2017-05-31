setlocal
if exist "c:\program files (x86)" goto x64
REM if exist "c:\program files\microsoft visual studio 10.0" goto x32_10
if exist "c:\program files\microsoft visual studio 8" goto x32_8
call "c:\program files\microsoft visual studio 9.0\vc\bin\x86_amd64\vcvarsx86_amd64"
goto makeit
:x32_10
call "c:\program files\microsoft visual studio 10.0\vc\bin\x86_amd64\vcvarsx86_amd64"
goto makeit
:x32_8
call "c:\program files\microsoft visual studio 8\vc\bin\x86_amd64\vcvarsx86_amd64"
goto makeit
:x64
REM if exist "c:\program files (x86)\microsoft visual studio 10.0" goto x64_10
if exist "c:\program files\microsoft visual studio 8" goto x64_8
call "c:\program files\microsoft visual studio 9.0\vc\bin\amd64\vcvarsamd64"
goto makeit
:x64_10
call "c:\program files\microsoft visual studio 10.0\vc\bin\amd64\vcvarsamd64"
goto makeit
:x64_8
call "c:\program files\microsoft visual studio 8\vc\bin\amd64\vcvarsamd64"
:makeit
@set path=y:\downloads\hjwasm230_x86;%path%
nmake /f compil64
endlocal
