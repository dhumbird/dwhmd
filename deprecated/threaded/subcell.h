#ifndef SUBCELL_H
#define SUBCELL_H

#include <set>
#include "particle.h"
#include "consts.h"
#include "mathfun.h"

//this is a class defining the subcells used in the linked-cell method.

typedef set<particle*>::iterator SPI;

class subcell{
 public:
  double u;
  double Lx, Ly, Lz;
  //  double rij, rik, rjk, cosOjik, cosOijk, cosOikj;
  //double rija, rjka, rika, d3;
  //double hjik, hjik2, hijk, hijk2, hikj, hikj2;
  set<subcell*>::iterator c;
  //svector Rij, Rik, Rjk, Fij;	
  //svector Ruij, Ruik, Rujk;	
  set<particle*> atom;
  set<subcell*> next;
  //set<particle*> trips_list;
  SPI i, j, k, trip;
  subcell(){u=0;}
  set<subcell*>::iterator nbegin(){return next.begin();}
  set<subcell*>::iterator nend(){return next.end();}
  SPI abegin(){return atom.begin();}
  SPI aend(){return atom.end();}
  void clear(){atom.clear();}
  void insert(particle* p){atom.insert(p);}
  void ninsert(subcell* p){next.insert(p);}
  //static void *th_s_ForceEval(void *);
  //void s_ForceEval();
  static void *th_s_ReNeighbor(void *);
  void s_ReNeighbor();
  //double TwoBody(particle*, particle*);
  //double ThreeBody(particle*, particle*, particle*);
};


#endif

