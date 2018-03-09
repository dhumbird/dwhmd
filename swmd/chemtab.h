//chemtab.h
//physical properties

#ifndef CHEMTAB_H
#define CHEMTAB_H
using namespace std;
#include <map>
#include <sstream>
#include "sysfun.h"

//define some hard-coded defaults

//*****silicon (14)*******
#define M_14 28.086
#define RC_14 14.22
#define TB_14 1
//*****argon (18)*********
#define M_18 39.948
#define RC_18 25
#define TB_18 0
//*****helium (2)*********
#define M_2 4.0026
#define RC_2 25
#define TB_2 0
//*****fluorine (9)*******
#define M_9 18.998403
#define RC_9 19.104
#define TB_9 1
//*****neon (10)**********
#define M_10 20.17
#define RC_10 25
#define TB_10 0
//*****krypton (36)*******
#define M_36 83.8
#define RC_36 25
#define TB_36 0

typedef struct cheminfo{
  float mass;
  float cutoff;
  bool three_body;
};

void ChemDef();


#endif
