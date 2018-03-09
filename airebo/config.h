#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <set>
#include "atom.h"

typedef vector<atom>::iterator VAI;
typedef set<atom*>::iterator SAI;

#include <string>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include "consts.h"
#include "mathfun.h"
#include "strfun.h"
#include "sysfun.h"
#include "tricubic.h"
#include <sstream>
#include "c_subcell.h"
#include "c_cluster.h"

//************************************************************************
class config{

 private:
  double b1, b2, b3, beta1, beta2, beta3, aa, qq, alpha;
  svector Fij;

  double dlam1, dlam2, S, Sprime; 
  double VA_ij, VR_ij, dVA_ij, dVR_ij;

  void * this_cell;

 public:
  vector<atom>* atom_list;
  //atom* j; atom* k; atom* l; atom* ii;
  VAI begin, end;//, i, vai_j;
  double t, Lx, Ly, Lz, T, dt, dt2, dtsq2, u, ek;
  int N, NFixed, Nmax, Nx, Ny, Nz, Nc;
  string name;
  set<int> surf_list, sput_list;
  double Lcx, Lcy, Lcz, Lz_full, the_top;
  set<subcell*>::iterator c; 
  vector<subcell> cells;
  vector<subcell>::iterator vscit;
  SAI sai_i, sai_j;
  subcell* cells_it;
  int inert;
  vector<cluster> cfg_clusts;

  //******inline functions**********
  config(){init();}
  config(string scfg){init(); Load(scfg); Partition(); ReNeighbor();}
  ~config(){delete atom_list;}
  void Setdt(double del_t){dt=del_t; dt2=del_t/2; dtsq2=dt*dt2;}
  void PrintE(){printf("%12.7f %10.5f %10.5f %12.7f %6.2f\n",t,u,ek,u+ek,T);}
  void TimeStepInit(){u=0; ek=0; for (VAI i=begin; i<end; i++) i->ts_init();}
  void SetProps(){for (VAI i=begin; i<end; i++) i->SetProps();}
  double MinZ(){
    double z=1000; for (VAI w=begin; w<end; w++) z=mfmin(w->R.z,z); return z;}
  double MaxZ(){
    double z=-1000; for (VAI w=begin; w<end; w++) z=mfmax(w->R.z,z); return z;}
  subcell* WhichCell(svector& R){
    if (R.z > the_top) Partition();
    return &cells[(int)((R.x + Lx/2)/Lcx)%Nx +
                 Nx*(((int)(fabs((R.y - Ly/2)/Lcy))%Ny) +
                     Ny*((int)(fabs((R.z - the_top)/Lcz))%Nz))];
  }
  //******in c_atom.cpp************
  void resize(int);
  VAI atomix(int);
  bool erase(VAI);
  VAI append(short);
  //******in c_lcell.cpp************
  void Partition();
  //******in c_base.cpp*************
  void reset();
  void init();
  //******in c_cryst.cpp************
  void MakeCryst(float, short, short, short, float, short, int, string);
  void Maxwell(float,int);
  //******in c_force.cpp*************
  void ForceEval();
  double BondOrder(nbr*, double, bool);
  //******in c_io.cpp****************
  void Load(string);
  void Dump();
  //******in c_md.cpp****************
  void dtOptimize();
  void ReNeighbor();
  void Run(double, int, double, bool, bool, int, bool, bool,bool);
  void RunQuenched(double, int, double, bool, bool, int, bool, bool, 
		   bool, float);
  void Thermo(float, float, int, double, bool, bool, bool,double);
  //******in c_vv.cpp****************
  void FirstVV();
  void SecondVV();
  //******in c_uonly.cpp************
  double Uonly();
  //*****in c_pass.cpp**************
  double U_on_i(atom*);
  void SurfAcc(short);
  void Passivate(short);
  void UPath(VAI);
  //*****in c_ion.cpp***************
  VAI AddIon(float, string, float, bool, float, bool, float, int, bool);
  void DelIon(bool);
  //*****in c_sput.cpp**************
  void Visit(atom*, cluster*);
  void CheckSput();
  void ClustSput(cluster*);
  void Sputter(int);
};


#endif



