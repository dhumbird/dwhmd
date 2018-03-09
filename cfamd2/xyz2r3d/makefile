# general makefile for xyz2r3d
# (c) 1999 Cam Abrams & UC Regents

C_SRCS = args.c colormap.c dblmat.c et.c point.c r3d_utils.c xyz2r3d.c

OBJS    = $(C_SRCS:.c=.o)

CC = cc
RM = rm -f

D = 0
B3 = 1
V = 1.0
CFLAGS_GENERAL = -O3 -64 -D BUILD_NUMBER=`date +%e%m%y` -D DIAG=$(D) \
	 -D VERSION_NUMBER=$(V) -Xcpluscomm

CFLAGS = $(CFLAGS_GENERAL)

LFLAGS = -64 
LIBS	= -lm

xyz2r3d: $(OBJS)
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	$(RM) *.o
