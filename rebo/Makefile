
BINDIR=/home/humbird/md/bin

CSRCS = c_io.cpp bicubic.cpp tricubic.cpp\
c_force.cpp c_vv.cpp chemtab.cpp atom.cpp nbr.cpp c_md.cpp c_base.cpp\
sysfun.cpp svector.cpp strfun.cpp errors.cpp mathfun.cpp c_atom.cpp\
c_uonly.cpp c_cryst.cpp c_pass.cpp c_ion.cpp \
c_lcell.cpp c_sput.cpp

MSRCS = md_main.cpp m_run.cpp m_manip.cpp m_init.cpp m_sat.cpp m_ion.cpp\
m_bomb.cpp m_grow.cpp

COBJS = $(CSRCS:.cpp=.o)
MOBJS = $(MSRCS:.cpp=.o)
default: rebo
atom.o: atom.h atom.cpp chemtab.o nbr.h
svector.o: svector.h svector.cpp
chemtab.o: chemtab.h chemtab.cpp

CXX = g++
RM = rm -f

CPPFLAGS = -O3 -ffast-math -march=pentium4 -mcpu=pentium4 -mmmx -msse -msse2 -funroll-loops -mfpmath=sse
LIBS	= 

all: rebo pdb 

rebo: 	$(COBJS) $(MOBJS) 
	$(CXX) -o $@ $(COBJS) $(MOBJS) $(LIBS)

pdb:	$(COBJS) pdbdump.o
	$(CXX) -o pdbdump $(COBJS) pdbdump.o $(LIBS)

clean:  
	$(RM) *.o *~ rebo pdbdump 

install: rebo pdb
	mv rebo pdbdump $(BINDIR)

