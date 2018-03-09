
BINDIR=/home/humbird/md/bin

SRCS = i_atom.cpp i_config.cpp i_chemtab.cpp svector.cpp strfun.cpp \
sysfun.cpp mathfun.cpp cfginfo.cpp diag.cpp energy.cpp uptake.cpp xps.cpp \
depth.cpp traj.cpp

OBJS = $(SRCS:.cpp=.o)

default: cfginfo

CXX = g++
RM = rm -f
host:=$(shell hostname)

CPPFLAGS = -O3 -ffast-math -march=pentium4 -mcpu=pentium4 -funroll-loops -mfpmath=sse
ifeq ($(host),bedrock)
CPPFLAGS =  -O3 -ffast-math -march=pentium4 -mcpu=pentium4 -funroll-loops
endif

LIBS	= 

all: cfginfo

cfginfo: $(OBJS) 
	$(CXX) -o $@ $(OBJS) $(LIBS)

clean:  
	$(RM) *.o *~ cfginfo

install: cfginfo
	mv cfginfo $(BINDIR)

