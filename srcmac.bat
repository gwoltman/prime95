REM Create and transfer a tar file of the GIMPS Mac OS X sources

del srcmac

mkdir foobar
mkdir foobar\gwnum
mkdir foobar\macosx

copy prime95\macosx\factor32.o foobar\macosx
copy prime95\license.txt foobar\macosx
copy prime95\stress.txt foobar\macosx

xcopy gwnum\*.h foobar\gwnum
xcopy gwnum\*.c foobar\gwnum
xcopy gwnum\*.cpp foobar\gwnum
xcopy gwnum\*.ld foobar\gwnum
xcopy gwnum\macosx\*.* foobar\gwnum
copy gwnum\makemac foobar\gwnum\makefile

xcopy linux\*.h foobar\macosx
xcopy linux\*.c foobar\macosx
copy linux\release foobar\macosx\release
copy linux\makemac foobar\macosx\makefile
xcopy linux\*.txt foobar\macosx

tar cvf srcmac *.h *.c qd\*.*

cd foobar\macosx
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

tar rvf ..\srcmac gwnum macosx

del /Q gwnum\*
del /Q macosx\*
rmdir gwnum
rmdir macosx
cd ..
rmdir foobar

ftp -i -n -s:\rbin\srcmac 192.168.0.2
del srcmac
