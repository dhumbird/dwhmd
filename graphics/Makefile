
BINDIR=/Users/humbird/bin


SRCS= rgb.cpp r3dfun.cpp rcfg.cpp svector.cpp miscfun.cpp \
ratom.cpp
OBJS= $(SRCS:.cpp=.o)

PDB_SRCS=rcfg.cpp ratom.cpp svector.cpp miscfun.cpp
PDB_OBJS=$(PDB_SRCS:.cpp=.o)

CXX = g++
RM = rm -f

CPPFLAGS = -O3 -ffast-math

LIBS	= -lm

all: cfg2r3d pdbdump

cfg2r3d: $(OBJS) cfg2r3d.o
	$(CXX) -o cfg2r3d $(OBJS) cfg2r3d.o $(LIBS)

pdbdump: $(PDB_OBJS) pdbdump.o
	$(CXX) -o pdbdump $(PDB_OBJS) pdbdump.o $(LIBS)
clean:  
	$(RM) *.o *~ core 

install: all
	mv cfg2r3d pdbdump $(BINDIR)
