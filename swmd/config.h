#ifndef CONFIG_H
#define CONFIG_H

using namespace std;

#include <vector>
#include <set>
#include "particle.h"

typedef vector<particle>::iterator VPI;
typedef set<particle*>::iterator SPI;

#include <string>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include "consts.h"
#include "mathfun.h"
#include "strfun.h"
#include "c_subcell.h"
#include "c_cluster.h"

//************************************************************************
class config{
 private:
  svector Rij, Rik, Rjk, Fij;	
  svector Ruij, Ruik, Rujk;	
  double rij, rik, rjk, cosOjik, cosOijk, cosOikj;
  double cosOjik2, cosOijk2, cosOikj2;
  double rija, rjka, rika, rija2, rjka2, rika2;
  double hjik, hjik2, hijk, hijk2, hikj, hikj2;
  double expa1jik, expa1ijk, expa1ikj;
  double expa2jik, expa2ijk, expa2ikj;
  
 public:
  vector<particle>* atom;
  VPI begin, end, i, j;
  SPI i0, j0, k0, trip;
  set<subcell*>::iterator c;
  set<particle*> trips_list;
  vector<subcell> cells;
  vector<subcell>::iterator cells_it;
  vector<cluster> cfg_clusts;
  double t, Lx, Ly, Lz, T, dt, dt2, dtsq2, u, ek;
  double Lcx, Lcy, Lcz, Lz_full, the_top;
  float  bhb_dt;
  void * this_cell;
  int N, NFixed, Nmax, Nx, Ny, Nz, Nc;
  string name;
  short mat; 
  set<int> sput_list;
  int inert;

  //******inline functions**********
  config(){init();}
  config(string scfg){init(); Load(scfg); Partition(); ReNeighbor();}
  ~config(){delete atom;}
  void Setdt(double del_t){dt=del_t; dt2=del_t/2; dtsq2=dt*dt2;}
  void PrintE(){printf("%f\t%8.3f\t%8.3f\t%8.3f\t%6.2f\n", t, u, ek, u+ek, T);}
  void TimeStepInit(){u=0; ek=0; for (i=begin; i<end; i++) i->ts_init();}
  double MaxZ(){
    double z=-1000; for (VPI w=begin; w<end; w++) z=w->R.z>?z; return z;}
  double MinZ(){
    double z=1000; for (VPI w=begin; w<end; w++) z=w->R.z<?z; return z;}
  void SetProps(){for (i=begin; i<end; i++) i->SetProps();}
  subcell* WhichCell(svector& R){
    if (R.z > the_top) Partition();
    return &cells[(int)((R.x + Lx/2)/Lcx)%Nx +
		 Nx*(((int)(fabs((R.y - Ly/2)/Lcy))%Ny) +
		     Ny*((int)(fabs((R.z - the_top)/Lcz))%Nz))];
  }
  //******in c_atom.cpp*************
  bool erase(VPI);
  VPI append(short);
  VPI atomix(int);
  void resize(int);
  //******in c_base.cpp*************
  void init();
  void reset();
  //******in c_cryst.cpp*************
  void MakeCryst(float, short, short, short, float, short, int, string);
  void AddFixed(float);
  void CheckDepth();
  //******in c_force.cpp*************
  void ForceEval();
  double TwoBody(particle*, particle*);
  double ThreeBody(particle*, particle*, particle*);
  double u_TwoBody(particle*, particle*);
  double u_ThreeBody(particle*, particle*, particle*);
  //******in c_si.cpp****************
  double SiSi(particle*, particle*);
  double SiSiSi(particle*, particle*, particle*);
  //******in c_f.cpp*****************
  double FF(particle*, particle*);
  double FFF(particle*, particle*, particle*);
  //******in c_sif.cpp***************
  double FSi(particle*, particle*);
  double FSiSi(particle*, particle*, particle*);
  double FFSi(particle*, particle*, particle*);
  //******in c_mo.cpp****************
  double HeSi(particle*, particle*);
  double NeSi(particle*, particle*);
  double ArSi(particle*, particle*);
  double KrSi(particle*, particle*);
  double HeF(particle*, particle*);
  double NeF(particle*, particle*);
  double ArF(particle*, particle*);
  //******in c_io.cpp****************
  void Load(string);
  void Dump();
  //******in c_ion.cpp***************
  VPI AddIon(float, int, float, bool, float, bool, float, int, bool);
  void DelIon(bool);
  //******in c_lcell.cpp*************
  void Partition();
  //******in c_md.cpp****************
  void ReNeighbor();
  void Thermo(float, float, int, double, bool, bool, bool);
  void Run(double, int, double, bool, bool, int, bool, bool,bool);
  void RunQuenched(double, int, double, bool, bool, int, bool, bool,
		   bool,float);
  void dtOptimize();
  //******in c_sput.cpp**************
  void Visit(particle*, cluster*);
  void BoundNeighbors();
  void CheckSput();
  void Sputter(int);
  void ClustSput(cluster*);
  void Sweep();
  //******in c_vv.cpp****************
  void FirstVV();
  void SecondVV();
  //******in c_sort.cpp*************
  void Ek_sort_asc();
  void Ek_sort_desc();
  void Rz_sort_asc();
  void Rz_sort_desc();
  void ix_sort_asc();
  void ix_sort_desc();
  //******in c_pass.cpp*************
  double FindWell(short,svector*,float,int);
  double U_on_i(VPI);
  //******in c_top.cpp**************
  double AvgTop();
  //******in c_uonly.cpp************
  double u_SiSi();
  double u_SiSiSi();
  double u_FF();
  double u_FFF();
  double u_FSi();
  double u_FSiSi(particle*, particle*, particle*);
  double u_FFSi(particle*, particle*, particle*);
};


#endif



