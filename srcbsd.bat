REM Create a tar file of the GIMPS sources to put on a floppy
REM and take to a FreeBSD machine.

REM The Linux version must be built before the FreeBSD version.  The
REM Linux object files in ELF format (*.o) must have been copied back
REM to this computer.

del srcx
copy prime95\license.txt linux
copy prime95\stress.txt linux
cd linux
mkdir foo
copy release foo
copy makebsd foo
copy *.txt foo
copy menu.o foo
copy prime.o foo
call nocr release
call nocr makebsd
call nocr whatsnew.txt
call nocr readme.txt
call nocr undoc.txt
call nocr license.txt
call nocr stress.txt
del menu.o
del prime.o
cd ..
tar cvf srcx *.h *.c linux\*.c linux\*.h linux\*.txt linux\make* linux\*.o linux\rel*
REM tar rvf srcx linuxprp\*.c linuxprp\*.h linuxprp\*.txt linuxprp\make* linuxprp\rel*
cd linux
copy foo\release
copy foo\makebsd
copy foo\*.txt
copy foo\menu.o
copy foo\prime.o
del foo\release
del foo\makebsd
del foo\*.txt
del foo\menu.o
del foo\prime.o
rmdir foo
cd ..
del linux\license.txt
del linux\stress.txt
