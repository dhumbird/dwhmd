# linux makefile for md series 2
# (c) 1999 Cam Abrams & UC Regents

C_SRCS = args.c atom.c bicubic.c cfg.c chem.c clust.c cryst.c dblmat.c domd.c \
         point.c push.c tbt_sicf.c tricubic.c ion.c

OBJS    = $(C_SRCS:.c=.o)
CC = /usr/bin/gcc
RM = rm -f

D = 0
B3 = 1
V = 1.0
MXA=1000
CFLAGS_GENERAL = -O3 -D BUILD_NUMBER=692001 -D DIAG=$(D) \
	 -D VERSION_NUMBER=$(V) -DTHREE_BODY=$(B3) -DMXA=$(MXA) -ffast-math \
	-march=pentium4 -mcpu=pentium4 -mmmx -msse -msse2 -funroll-loops \
	-mfpmath=sse

CFLAGS = $(CFLAGS_GENERAL)

LFLAGS = 
LIBS	= -lm

md2:	$(OBJS) main.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) main.o $(LIBS)

cfginfo: $(OBJS) cfginfo.o
	$(CC) $(LFLAGS) -o $@ $(OBJS) cfginfo.o $(LIBS)

cfg.o:   cfg.h
tbt_sicf.o: cfg.h tbt_sicf.h

cfg2r3d: $(OBJS) $(OR3D)
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(OR3D) $(LIBS)

clean:  
	$(RM) *.o

install: md2
	mv -f md2 /home/humbird/bin/cfamd2
