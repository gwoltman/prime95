# Makefile for Linux mprime
# by Jeremy Buhler (jbuhler@cs.washington.edu)
#
# NOTE: To convert the included object files to ELF, your binutils 
# installation must be able to read the coff-i386 object format.  
#
# Assuming your binutils are dynamically linked, you can add COFF support
# by configuring binutils with the switch '--enable-targets=i386-go32'
# and rebuilding libbfd.so .

CC = gcc
CFLAGS = -I.. -O2 -march=i486 -malign-double -DSTPCHK

CPP = g++
CPPFLAGS = -I.. -I../qd -O2 -march=i486 -malign-double -DSTPCHK

LFLAGS = -Wl,-M
LIBS   = -lm
LIBS2  = -lm -static

COFFOBJS = mult.o mult1.o mult2.o mult2a.o mult3.o mult3a.o mult4.o mult2p.o mult2ap.o mult3p.o mult3ap.o mult4p.o mult1aux.o mult2aux.o xmult1.o xmult1ax.o xmult2.o xmult2a.o xmult2ax.o xmult3.o xmult3a.o xmult3ax.o ecmhelp.o cpuid.o
FACTOROBJ = factor64.o
DUMMYOBJS = dummy0.o dummy4.o dummy8.o dummy12.o dummy16.o dummy20.o dummy24.o dummy28.o dummyt0.o dummyt4.o dummyt8.o dummyt12.o dummyt16.o dummyt20.o dummyt24.o dummyt28.o
LINUXOBJS = prime.o menu.o
COMMONOBJS = cpuid.o gwnum.o gwutil.o gwdbldbl.o
EXE      = mprime
EXE2	 = sprime

# insert proper dummy file to get mult.o and factor64.o's data aligned
# on 32 byte boundaries.  Use dummyNN.o to align mult.o and dummytNN.o
# to align factor64.o

DUMMY1   = dummy28.o
DUMMY2   = dummyt0.o
DUMMY3   = dummy28.o
DUMMY4   = dummyt12.o

#########################################################################

$(EXE): security $(LINUXOBJS) $(COMMONOBJS) $(COFFOBJS) $(FACTOROBJ) $(DUMMYOBJS)
	$(CC) $(LFLAGS) -o $(EXE) $(LINUXOBJS) $(COMMONOBJS) $(DUMMY1) $(COFFOBJS) $(DUMMY2) $(FACTOROBJ) $(LIBS)

$(EXE2): security $(LINUXOBJS) $(COMMONOBJS) $(COFFOBJS) $(FACTOROBJ) $(DUMMYOBJS)
	$(CC) $(LFLAGS) -o $(EXE2) $(LINUXOBJS) $(COMMONOBJS) $(DUMMY3) $(COFFOBJS) $(DUMMY4) $(FACTOROBJ) $(LIBS2)

security:
	[ ! -e ../security.h ] && touch ../security.h || true
	[ ! -e ../security.c ] && touch ../security.c || true

clean:
	rm -f $(EXE) $(EXE2) $(LINUXOBJS) $(COFFOBJS) $(DUMMYOBJS)

$(COFFOBJS) $(FACTOROBJ) $(DUMMYOBJS): %.o : %.obj
	objcopy -Oelf32-i386 -Icoff-i386 $< $@

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
