//particle.h
//class listing for particle class

#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <set>
#include <iostream>
#include "svector.h"
#include "chemtab.h"

class particle{

 public:

  //printed members
  svector R; //position
  svector V; //velocity
  svector O; //origin
  int ix;
  short id;
  bool is_fixed;
  int xpass, ypass;

  //non-printed members
  svector P, F;
  void * my_cell;
  double m;
  double u;
  set<particle*> nlist;
  double rc; //squared cutoff, sq. ang
  bool three_body;
  void * p_clust;

  particle();
  void SetProps();
  void PosUpdate(double&, double&, double&, double&);
  void ts_init() {u=0; F.clear(); p_clust=NULL;}
  double Ek() {return 0.5*m*V.sqmag();}
  friend ostream& operator << (ostream&, particle&);
  friend istream& operator >> (istream&, particle&);
};

#endif
