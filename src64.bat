REM Create and transfer a tar file of the GIMPS linux x86-64 sources

del src64

mkdir foobar
mkdir foobar\gwnum
mkdir foobar\linux

copy prime95\linux64\factor64.o foobar\linux
copy prime95\license.txt foobar\linux
copy prime95\stress.txt foobar\linux

xcopy gwnum\*.h foobar\gwnum
xcopy gwnum\*.c foobar\gwnum
xcopy gwnum\*.cpp foobar\gwnum
xcopy gwnum\*.ld foobar\gwnum
xcopy gwnum\linux64\*.a foobar\gwnum
copy gwnum\make64 foobar\gwnum\makefile

xcopy linux\*.h foobar\linux
xcopy linux\*.c foobar\linux
copy linux\releas64 foobar\linux\release
copy linux\make64 foobar\linux\makefile
xcopy linux\*.txt foobar\linux

tar cvf src64 *.h *.c qd/*.*

cd foobar\linux
call nocr release
call nocr makefile
call nocr whatsnew.txt
call nocr readme.txt
call nocr undoc.txt
call nocr license.txt
call nocr stress.txt
cd ..\gwnum
call nocr makefile
cd ..

tar rvf ..\src64 gwnum linux

del /Q gwnum\*
del /Q linux\*
rmdir gwnum
rmdir linux
cd ..
rmdir foobar

ftp -i -n -s:\rbin\src64 192.168.0.8
