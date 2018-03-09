
BINDIR=/usr/home/humbird/md/bin


CFG_SRCS = c_base.cpp c_io.cpp c_cryst.cpp c_lcell.cpp c_md.cpp c_force.cpp \
c_ion.cpp c_atom.cpp c_info.cpp c_vv.cpp c_sorte.cpp c_sortr.cpp c_sortix.cpp \
subcell.cpp
CFG_OBJS = $(CFG_SRCS:.cpp=.o)

MD_SRCS=chemtab.cpp particle.cpp svector.cpp sysfun.cpp strfun.cpp \
mathfun.cpp bombard.cpp errors.cpp edep.cpp 
MD_OBJS    = $(MD_SRCS:.cpp=.o)

R3D_SRCS=particle.cpp chemtab.cpp svector.cpp sysfun.cpp strfun.cpp \
mathfun.cpp errors.cpp g/rgb.cpp g/r3d_util.cpp
R3D_OBJS= $(R3D_SRCS:.cpp=.o)

INF_SRCS = chemtab.cpp particle.cpp svector.cpp sysfun.cpp strfun.cpp \
mathfun.cpp errors.cpp info/velocity.cpp info/position.cpp info/energy.cpp \
info/rdf.cpp info/ndens.cpp
INF_OBJS = $(INF_SRCS:.cpp=.o)


particle.o: particle.h particle.cpp chemtab.o
edep.o: edep.cpp edep.h consts.h
svector.o: svector.h svector.cpp
mkelip.o: mkelip.cpp strfun.o
chemtab.o: chemtab.h chemtab.cpp
g/cfg2r3d.o: svector.o

CC = g++
RM = rm -f

CPPFLAGS = -O3 -ffast-math

LIBS	= -lm -lpthread

all: md cfginfo cfg2r3d

md: 	$(MD_OBJS) $(CFG_OBJS) md_main.o
	$(CC) -o $@ $(MD_OBJS) $(CFG_OBJS) md_main.o $(LIBS)

cfginfo: $(INF_OBJS) $(CFG_OBJS) info/cfginfo.o
	$(CC) -o info/$@ $(INF_OBJS) $(CFG_OBJS) info/cfginfo.o \
	$(LIBS)

cfg2r3d: $(R3D_OBJS) $(CFG_OBJS) g/cfg2r3d.o
	$(CC) -o g/$@ $(R3D_OBJS) $(CFG_OBJS) g/cfg2r3d.o $(LIBS)

clean:  
	$(RM) *.o *~ core g/*.o g/*~ g/core md g/cfg2r3d info/cfginfo \
	info/core info/*.o info/*~

install: md cfginfo cfg2r3d 
	mv md info/cfginfo g/cfg2r3d $(BINDIR)


