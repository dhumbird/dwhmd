# SGI makefile for md series 2
# (c) 1999 Cam Abrams & UC Regents

BINDIR=/Users/humbird/md/bin

C_SRCS = args.c atom.c bicubic.c cfg.c chem.c clust.c cryst.c \
	 dblmat.c domd.c ion.c point.c push.c tbt_sicf.c tricubic.c

C_HDRS = $(C_SRCS:.c=.h) units.h sicf_params.h tricof.h genforce.h

C_R3D	= cfg2r3d.c colormap.c r3d_utils.c
C_R3DH	= colormap.h r3d_utils.h

OBJS    = $(C_SRCS:.c=.o)
OR3D    = $(C_R3D:.c=.o)

CC = gcc
RM = rm -f

D = 0
B3 = 1
V = 1.0
MXA=2000
CFLAGS_GENERAL = -O3 -D BUILD_NUMBER=`date +%e%m%y` -D DIAG=$(D) \
	 -D VERSION_NUMBER=$(V) -DTHREE_BODY=$(B3) -DTERSOFF_NEIGHBORS \
	 -DMXA=$(MXA) --fast-math \
	-funroll-loops 

CFLAGS = $(CFLAGS_GENERAL)

LFLAGS =  
LIBS	= -lm

md2:	$(OBJS) main.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) main.o $(LIBS)

package:	$(C_SRCS) main.c $(C_HDRS) $(C_R3D) $(C_R3DH)\
		 README manual.ps aSi.cfg makefile.sgi makefile.linux
	tar cf md2.tar $?
	gzip md2.tar

cfg.o:   cfg.h
tbt_sicf.o: cfg.h tbt_sicf.h sicf_params.h

cfgfix:	$(OBJS) cfgfix.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) cfgfix.o $(LIBS)

debug:	$(OBJS) debug.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) debug.o $(LIBS)

cgmin:	$(OBJS) cfg_cgmin.o  main_cgmin.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) cfg_cgmin.o main_cgmin.o $(LIBS)
	mv $@ $(BINDIR)

sat:	$(OBJS) cfg_sat.o  cfg_cgmin.o main_sat.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) cfg_sat.o cfg_cgmin.o main_sat.o $(LIBS)
	mv $@ $(BINDIR)

install: md2
	mv $? $(BINDIR)

cfg2r3d: $(OBJS) $(OR3D)
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(OR3D) $(LIBS)
	mv $@ $(BINDIR)

cfginfo: $(OBJS) $(OBJS) cfginfo.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) cfginfo.o $(LIBS)
	mv $@ $(BINDIR)

cfgrdf: $(OBJS) $(OBJS) cfgrdf.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) cfgrdf.o $(LIBS)
	mv $@ $(BINDIR)

nbrof: $(OBJS) $(OBJS) nbrof.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) nbrof.o $(LIBS)
	mv $@ $(BINDIR)

siff: $(OBJS) $(OBJS) siff.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) siff.o $(LIBS)
	mv $@ $(BINDIR)

fsif: $(OBJS) $(OBJS) fsif.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) fsif.o $(LIBS)
	mv $@ $(BINDIR)

clean:
	$(RM) *.o
