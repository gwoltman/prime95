# Makefile for OS/2 mprime
# to compile with Watcom 11c / OpenWatcom
# use with GNUMAKE
# (c) Max Alekseyev

CC = wcc386
CPP = wpp386
CFLAGS = /ox /d0 /I..;..\linux;..\gwnum;..\qd /5s /D__IBMC__ /DTCPV40HDRS
LINK = wlink
LFLAGS = LIBPath ..\prime95;..\linux;..\gwnum op m import DosQuerySysState DOSCALLS.368

# debug version
# CFLAGS = /ox /d3 /I.. /5s /D__IBMC__ /D_DEBUG
# LFLAGS = LIBPath ..\prime95 op m debug all import DosQuerySysState DOSCALLS.368

COFFOBJS = cpuidhlp.obj,factor64.obj,gianthlp.obj,mult.obj,mult1.obj,mult1aux.obj,mult2.obj,mult2a.obj,mult2ap.obj,mult2aux.obj,mult2p.obj,mult3.obj,mult3a.obj,mult3ap.obj,mult3p.obj,mult4.obj,mult4p.obj,xmult1.obj,xmult1ax.obj,xmult2.obj,xmult2a.obj,xmult2ax.obj,xmult3.obj,xmult3a.obj,xmult3ax.obj
DUMMYOBJS = dummy4.obj,dummy8.obj,dummy12.obj,dummy16.obj,dummy20.obj,dummy24.obj,dummy28.obj
PRIMEOBJS = prime.obj,menu.obj

# PRIMESRC = prime.c,menu.c
# LIBS = tcpip32.lib
LIBS = so32dll.lib,tcp32dll.lib
EXE = primeos2.exe

$(EXE): prime.obj menu.obj gwdbldbl.obj	# $(PRIMEOBJS),$(COFFOBJS),$(DUMMYOBJS)
	$(LINK) $(LFLAGS) file prime.obj,menu.obj,gwdbldbl.obj name $(EXE) libfile $(COFFOBJS) library $(LIBS) 

prime.obj: ..\linux\prime.c
	$(CC) $(CFLAGS) $<

menu.obj: ..\linux\menu.c
	$(CC) $(CFLAGS) $<

gwdbldbl.obj: ..\gwnum\gwdbldbl.cpp
	$(CPP) $(CFLAGS) $<

clean:
	del $(EXE) *.obj
