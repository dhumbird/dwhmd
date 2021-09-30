
BINDIR=/home/dhumbird/md/bin

SRCS = i_atom.cpp i_config.cpp i_chemtab.cpp svector.cpp strfun.cpp \
sysfun.cpp mathfun.cpp cfginfo.cpp diag.cpp energy.cpp uptake.cpp xps.cpp \
depth.cpp traj.cpp

OBJS = $(SRCS:.cpp=.o)

default: cfginfo

CXX = g++
RM = rm -f
host:=$(shell hostname)

CPPFLAGS = -O3 -ffast-math -fpermissive

LIBS	= 

all: cfginfo

cfginfo: $(OBJS) 
	$(CXX) -o $@ $(OBJS) $(LIBS)

clean:  
	$(RM) *.o *~ cfginfo

install: cfginfo
	mv cfginfo $(BINDIR)

