REM Create a tar file of the PRP linux sources to put on a floppy
REM and take to a linux machine.
xcopy prime95\*.obj linuxprp
tar cvzf srcy *.h *.c linuxprp\*.c linuxprp\*.h linuxprp\*.txt linuxprp\makefile linuxprp\*.obj linuxprp\release
del linuxprp\*.obj
