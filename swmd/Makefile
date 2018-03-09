
BINDIR=/home/humbird/md/bin


CFG_SRCS = c_sort.cpp c_base.cpp c_io.cpp c_cryst.cpp c_lcell.cpp c_md.cpp \
c_force.cpp c_ion.cpp c_atom.cpp c_vv.cpp c_sput.cpp c_si.cpp c_f.cpp \
c_mo.cpp c_sif.cpp c_pass.cpp c_top.cpp c_uonly.cpp chemtab.cpp particle.cpp \
svector.cpp sysfun.cpp strfun.cpp errors.cpp
CFG_OBJS = $(CFG_SRCS:.cpp=.o)

R3D_SRCS= g/rgb.cpp g/r3d_util.cpp g/colors.cpp
R3D_OBJS= $(R3D_SRCS:.cpp=.o)

INF_SRCS =info/velocity.cpp info/position.cpp info/energy.cpp info/ndens.cpp \
info/diffusion.cpp info/diag.cpp info/amorph.cpp  \
info/rdf.cpp info/uptake.cpp 
INF_OBJS = $(INF_SRCS:.cpp=.o)

PDB_SRCS=g/pdbdump.cpp
PDB_OBJS=$(PDB_SRCS:.cpp=.o)

MAIN_SRCS = md_main.cpp m_init.cpp m_sort.cpp m_ion.cpp m_bomb.cpp \
m_run.cpp m_manip.cpp m_sput.cpp  m_sat.cpp m_cgmin.cpp
MAIN_OBJS = $(MAIN_SRCS:.cpp=.o)

particle.o: particle.h particle.cpp chemtab.o
edep.o: edep.cpp edep.h consts.h
svector.o: svector.h svector.cpp
mkelip.o: mkelip.cpp strfun.o
chemtab.o: chemtab.h chemtab.cpp
g/cfg2r3d.o: svector.o

CXX = g++
RM = rm -f
host:=$(shell hostname)

CPPFLAGS = -O3 -ffast-math
ifeq ($(host),raptor)
CPPFLAGS =  -O3 -ffast-math -march=pentium4 -mcpu=pentium4 -funroll-loops -mfpmath=sse
endif
ifeq ($(host),ihml)
CPPFLAGS = -O3 -ffast-math -march=athlon-mp -mcpu=athlon-mp -funroll-loops -mfpmath=sse
endif
ifeq ($(host),bedrock)
CPPFLAGS =  -O3 -ffast-math -march=pentium4 -mcpu=pentium4 -funroll-loops
endif

LIBS	= -lm

all: swmd cfginfo cfg2r3d pdbdump

swmd: 	$(MAIN_OBJS) $(CFG_OBJS) 
	$(CXX) -o $@ $(CFG_OBJS) $(MAIN_OBJS) $(LIBS)

cfginfo: $(INF_OBJS) $(CFG_OBJS) info/cfginfo.o
	$(CXX) -o info/swcfginfo $(INF_OBJS) $(CFG_OBJS) info/cfginfo.o \
	$(LIBS)

cfg2r3d: $(R3D_OBJS) $(CFG_OBJS) g/cfg2r3d.o
	$(CXX) -o g/swcfg2r3d $(R3D_OBJS) $(CFG_OBJS) g/cfg2r3d.o $(LIBS)

pdbdump: $(PDB_OBJS) $(CFG_OBJS) 
	$(CXX) -o g/swpdbdump $(PDB_OBJS) $(CFG_OBJS) $(LIBS)
clean:  
	$(RM) *.o *~ core g/*.o g/*~ g/core swmd g/swcfg2r3d info/swcfginfo \
        info/core info/*.o info/*~

mdclean:
	$(RM) *.o *~ core 

install: swmd 
	mv swmd $(BINDIR)
