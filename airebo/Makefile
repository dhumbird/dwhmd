
BINDIR=/home/humbird/md/bin

CSRCS = c_io.cpp bicubic.cpp tricubic.cpp\
c_force.cpp c_vv.cpp chemtab.cpp atom.cpp nbr.cpp c_md.cpp c_base.cpp\
sysfun.cpp svector.cpp strfun.cpp errors.cpp mathfun.cpp c_atom.cpp\
c_uonly.cpp c_cryst.cpp c_pass.cpp c_ion.cpp \
c_lcell.cpp c_sput.cpp ljnbr.cpp

MSRCS = md_main.cpp m_run.cpp m_manip.cpp m_init.cpp m_sat.cpp m_ion.cpp\
m_bomb.cpp m_grow.cpp

ISRCS = info/cfginfo.cpp info/uptake.cpp info/diag.cpp info/energy.cpp

COBJS = $(CSRCS:.cpp=.o)
MOBJS = $(MSRCS:.cpp=.o)
IOBJS = $(ISRCS:.cpp=.o)

default: rebo
atom.o: atom.h atom.cpp chemtab.o nbr.h
svector.o: svector.h svector.cpp
chemtab.o: chemtab.h chemtab.cpp

CXX = /usr/local/bin/g++
RM = rm -f

CPPFLAGS = -O3 -ffast-math -march=pentium4 -mcpu=pentium4 -mmmx -msse -msse2 -funroll-loops -mfpmath=sse
LIBS	= 

all: rebo pdb cfginfo

rebo: 	$(COBJS) $(MOBJS) 
	$(CXX) -o $@ $(COBJS) $(MOBJS) $(LIBS)

pdb:	$(COBJS) pdbdump.o
	$(CXX) -o pdbdump $(COBJS) pdbdump.o $(LIBS)

cfginfo: $(IOBJS) $(COBJS) 
	$(CXX) -o info/$@ $(COBJS) $(IOBJS) $(LIBS)

clean:  
	$(RM) *.o *~ rebo pdbdump info/*.o info/*~

install: rebo pdb cfginfo
	mv rebo pdbdump info/cfginfo $(BINDIR)

