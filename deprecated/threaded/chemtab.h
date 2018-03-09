//chemtab.h
//physical properties

#ifndef CHEMTAB_H
#define CHEMTAB_H

#include <map>
#include <sstream>
#include "sysfun.h"

//define some hard-coded defaults

//*****silicon*******
#define M_SI 28.086
#define RC_SI 14.22
#define TB_SI 1
//*****argon*********
#define M_AR 39.948
#define RC_AR 25
#define TB_AR 0

typedef struct cheminfo{
  float mass;
  float cutoff;
  bool three_body;
};

void ChemDef();


#endif
