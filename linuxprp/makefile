# Makefile for Linux PRP
#
# Assuming your binutils are dynamically linked, you can add COFF support
# by configuring binutils with the switch '--enable-targets=i386-go32'
# and rebuilding libbfd.so .

CC = gcc
CFLAGS = -I.. -I../gwnum -O2 -march=i486 -malign-double

CPP = g++
CPPFLAGS = -I.. -I../gwnum -O2 -march=i486 -malign-double

LFLAGS = -Wl,-M
LIBS   = ../gwnum/gwnum.lib -lm -lstdc++
LIBS2  = ../gwnum/gwnum.lib -lm -lstdc++ -static

LINUXOBJS = prime.o menu.o
EXE      = prp
EXE2	 = sprp
DUMMY1   = ../linux/dummy24.o
DUMMY2   = ../linux/dummy24.o

#########################################################################

$(EXE): $(LINUXOBJS)
	$(CC) $(LFLAGS) -o $(EXE) $(LINUXOBJS) $(DUMMY1) $(LIBS)

$(EXE2): $(LINUXOBJS)
	$(CC) $(LFLAGS) -o $(EXE2) $(LINUXOBJS) $(DUMMY2) $(LIBS2)

clean:
	rm -f $(EXE) $(EXE2) $(LINUXOBJS)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
