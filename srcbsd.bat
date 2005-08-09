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
call nocr release
call nocr makebsd
call nocr whatsnew.txt
call nocr readme.txt
call nocr undoc.txt
call nocr license.txt
call nocr stress.txt
cd ..
cd gwnum
mkdir foo
copy makefile foo
copy makebsd foo
call nocr makefile
call nocr makebsd
cd ..
tar cvf srcx *.h *.c *.cpp linux\*.c linux\*.h linux\*.txt linux\make* linux\*.o linux\rel*
tar rvf srcx linuxprp\*.c linuxprp\*.h linuxprp\*.txt linuxprp\make* linuxprp\rel*
tar rvf srcx gwnum\*.h gwnum\*.c gwnum\*.cpp gwnum\make* gwnum\*.ld gwnum\*.o qd\*.*
cd linux
copy foo\release
copy foo\makebsd
copy foo\*.txt
del foo\release
del foo\makebsd
del foo\*.txt
rmdir foo
cd ..
cd gwnum
copy foo\makefile
copy foo\makebsd
del foo\makefile
del foo\makebsd
rmdir foo
cd ..
del linux\license.txt
del linux\stress.txt
