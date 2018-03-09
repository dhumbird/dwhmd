#ifndef CONFIG_H
#define CONFIG_H

#ifdef __linux__
#  define _REENTRANT
#  define _P __P
#endif

#include <vector>
#include <string>
#include <fstream>
#include <set>
#include <algorithm>
#include <cstdio>
#include <pthread.h>
#include "particle.h"
#include "subcell.h"
#include "consts.h"
#include "strfun.h"
#include "mathfun.h"

typedef vector<particle>::iterator VPI;
typedef set<particle*>::iterator SPI;


//************************************************************************
class config{

 public:
  vector<particle>* atom;
  VPI begin, end, i, j, p;
  SPI i0, j0, k0, trip;
  set<subcell*>::iterator c;
  set<particle*> trips_list;
  vector<subcell> cells;
  vector<subcell>::iterator cells_it;
  set<subcell*> cellset1, cellset2;
  double t, Lx, Ly, Lz, T, dt, dt2, dtsq2, u, ek;
  double ek_av, u_av, d1, Lcx, Lcy, Lcz, sLz, sLz2, Lz2;
  float  bhb_dt;
  int N, NFixed, Nsput, a, x, y, z, Nx, Ny, Nz, Nc;
  string name;
  bool chem_is_def;
  short mat; //base chemical of the crystal
  pthread_t threads[2];
  int *status;

  //******inline functions**********
  ~config(){delete atom;}
  void Setdt(double del_t){dt=del_t; dt2=del_t/2; dtsq2=dt*dt2;}
  void PrintE(){printf("%f\t%e\t%e\t%e\t%f\n", t, u, ek, u+ek, T);}
  void RePartition(){cells.clear(); Partition();}
  //******in c_atom.cpp*************
  bool erase(VPI);
  void CopyAtom(VPI);
  VPI atomix(int);
  //******in c_base.cpp*************
  config();
  config(string);
  void resize(int);
  void init();
  //******in c_cryst.cpp*************
  void MakeCryst(float, short, short, short, float, bool, bool, short, int,
		 string);
  //******in c_force.cpp*************
  double TwoBody(particle*, particle*);
  double ThreeBody(particle*, particle*, particle*);
  void ForceEval(bool);
  //******in c_info.cpp**************
  double MaxZ();
  double MinZ();
  //******in c_io.cpp****************
  void Load(string);
  void Dump(bool);
  void SetProps();
  //******in c_ion.cpp***************
  void AddIon(float, short, float, bool, float, bool, float, int, bool, bool,
	      bool, bool, bool, float);
  void DelIon(bool);
  //******in c_lcell.cpp*************
  void Partition();
  void Allocate();
  //******in c_md.cpp****************
  void dtOptimize();
  void ReNeighbor();
  void TimeStepInit();
  void Thermo(float, float, int, double, bool, bool, bool);
  void Run(double, int, double, bool, bool, bool);
  bool CheckSput(VPI);
  //******in c_vv.cpp****************
  void FirstVV();
  void OFirstVV();
  void SecondVV();
  //******in c_sorte.cpp************
  void Ek_sort_asc();
  void Ek_sort_desc();
  void F_sort_asc();
  void F_sort_desc();
  //******in c_sortr.cpp************
  void Rz_sort_asc();
  void Rz_sort_desc();
  //void Rx_sort_asc();
  //void Rx_sort_desc();
  //void Ry_sort_asc();
  //void Ry_sort_desc();
  //******in c_sortix.cpp***********
  void ix_sort_asc();
  void ix_sort_desc();
};


#endif



