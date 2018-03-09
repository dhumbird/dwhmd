
BINDIR=/home/humbird/md/bin


CFG_SRCS = c_sort.cpp c_base.cpp c_io.cpp c_cryst.cpp c_lcell.cpp c_md.cpp \
c_force.cpp c_ion.cpp c_atom.cpp c_vv.cpp c_sput.cpp c_si.cpp c_f.cpp \
c_mo.cpp c_sif.cpp c_cl.cpp c_sicl.cpp c_pass.cpp c_top.cpp c_uonly.cpp
CFG_OBJS = $(CFG_SRCS:.cpp=.o)

MD_SRCS=chemtab.cpp particle.cpp svector.cpp sysfun.cpp strfun.cpp errors.cpp \
m_sat.cpp
MD_OBJS    = $(MD_SRCS:.cpp=.o)

R3D_SRCS=particle.cpp chemtab.cpp svector.cpp sysfun.cpp strfun.cpp \
errors.cpp g/rgb.cpp g/r3d_util.cpp g/colors.cpp
R3D_OBJS= $(R3D_SRCS:.cpp=.o)

INF_SRCS = chemtab.cpp particle.cpp svector.cpp sysfun.cpp strfun.cpp \
errors.cpp info/velocity.cpp info/position.cpp info/energy.cpp info/ndens.cpp \
info/diffusion.cpp info/diag.cpp info/amorph.cpp info/halogen.cpp \
info/rdf.cpp info/uptake.cpp info/dissoc.cpp info/sif3stat.cpp

INF_OBJS = $(INF_SRCS:.cpp=.o)

MAIN_SRCS = md_main.cpp m_init.cpp m_sort.cpp m_ion.cpp m_bomb.cpp m_edep.cpp \
m_run.cpp m_manip.cpp m_sput.cpp 

MAIN_OBJS = $(MAIN_SRCS:.cpp=.o)

particle.o: particle.h particle.cpp chemtab.o
edep.o: edep.cpp edep.h consts.h
svector.o: svector.h svector.cpp
mkelip.o: mkelip.cpp strfun.o
chemtab.o: chemtab.h chemtab.cpp
g/cfg2r3d.o: svector.o

CC = g++
RM = rm -f

CPPFLAGS = -O3 -ffast-math

LIBS	= -lm

all: md cfginfo cfg2r3d

md: 	$(MAIN_OBJS) $(CFG_OBJS) $(MD_OBJS) 
	$(CC) -o $@ $(MD_OBJS) $(CFG_OBJS) $(MAIN_OBJS) $(LIBS)

cfginfo: $(INF_OBJS) $(CFG_OBJS) info/cfginfo.o
	$(CC) -o info/$@ $(INF_OBJS) $(CFG_OBJS) info/cfginfo.o \
	$(LIBS)

cfg2r3d: $(R3D_OBJS) $(CFG_OBJS) g/cfg2r3d.o
	$(CC) -o g/$@ $(R3D_OBJS) $(CFG_OBJS) g/cfg2r3d.o $(LIBS)

clean:  
	$(RM) *.o *~ core g/*.o g/*~ g/core md g/cfg2r3d info/cfginfo \
	info/core info/*.o info/*~

mdclean:
	$(RM) *.o *~ core 

install: md cfginfo cfg2r3d 
	mv md info/cfginfo g/cfg2r3d $(BINDIR)

circtest: $(CFG_OBJS) $(MD_OBJS) circtest.o
	$(CC) -o circ $(CFG_OBJS) $(MD_OBJS) circtest.o $(LIBS)


