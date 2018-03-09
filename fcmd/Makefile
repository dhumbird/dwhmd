
BINDIR=/home/dhumbird/bin

SRCS = c_io.cpp bicubic.cpp tricubic.cpp\
c_force.cpp c_vv.cpp chemtab.cpp atom.cpp nbr.cpp c_md.cpp c_base.cpp\
sysfun.cpp svector.cpp strfun.cpp errors.cpp mathfun.cpp c_atom.cpp\
c_cryst.cpp c_ion.cpp c_lcell.cpp c_sput.cpp c_uonly.cpp c_pass.cpp\
md_main.cpp m_run.cpp m_manip.cpp m_init.cpp m_ion.cpp m_bomb.cpp\
m_mix.cpp m_cgmin.cpp m_edep.cpp

OBJS = $(SRCS:.cpp=.o)

default: fcmd
atom.o: atom.h atom.cpp nbr.h
svector.o: svector.h svector.cpp
chemtab.o: chemtab.h chemtab.cpp

CXX = g++
RM = rm -f
host:=$(shell hostname)

CPPFLAGS = -O3 -ffast-math -fpermissive


LIBS	= 

fcmd: 	$(OBJS)  
	$(CXX) -o $@ $(OBJS) $(LIBS)

clean:  
	$(RM) *.o *~ fcmd

install: fcmd
	mv fcmd $(BINDIR)

