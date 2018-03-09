#include "config.h"
#define SI_F_av 1.6
#define SI_F_max 1.81
#define SI_SI_max 2.77
#define SI_rad 1.2
#define F_rad 0.9

double sisi=0;
double sisisi=0;
double ff=0; double fff=0;
double fsi=0; double ffsi=0; double fsisi=0;

//***********************************************************************
double config::U_on_i(VPI s){
  if ((this_cell=WhichCell(s->R)) != s->my_cell){
    ((subcell*) s->my_cell)->erase(&(*s));
    ((subcell*) this_cell)->insert(&(*s));
  }
  subcell * scp=WhichCell(s->R);
  for (j0=s->nlist.begin(); j0!=s->nlist.end(); j0++)
    (*j0)->nlist.erase(&(*s));
  
  s->nlist.clear();
  for (c=scp->nbegin(); c!=scp->nend(); c++)
    for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
      if (&(*s)!=*j0)
	if ((s->R - (*j0)->R).minsqmag(Lx, Ly) < (s->rc >? (*j0)->rc)){
	  s->nlist.insert(*j0);
	  (*j0)->nlist.insert(&(*s));
	}
  double d=0;
  for(j0=s->nlist.begin(); j0!=s->nlist.end(); j0++){
    d+=u_TwoBody(&(*s),*j0);
    if (s->three_body){
      trips_list.clear();
      for (trip=(*j0)->nlist.begin(); trip!=(*j0)->nlist.end(); trip++)
	if (*trip!=&(*s)) trips_list.insert(*trip);
      trip=j0;
      for (trip++; trip!=s->nlist.end(); trip++)
	trips_list.insert(*trip);
      for(trip=s->nlist.begin(); trip!=j0; trip++)
	trips_list.erase(*trip);
      for(k0=trips_list.begin(); k0!=trips_list.end(); k0++)
	d+=u_ThreeBody(&(*s),*j0, *k0);
    }
  }
  cerr<<"Si-Si "<<sisi<<" \tF-F "<<ff<<endl;
  cerr<<"Si-Si-Si "<<sisisi<<" \tF-F-F "<<fff<<endl;
  cerr<<"Si-F "<<fsi<<" \tF-Si-Si "<<fsisi<<" \tF-F-Si "<<ffsi<<endl;
  return d;
}
//***********************************************************************
double config::FindWell(short id, svector* SV, float To, int xx){
   
  VPI newatom=append(id);
  double u_min=1000;
  VPI s=atomix(xx);
  svector P;
  float bl=SI_F_av;
  float bl_SQ=bl*bl;
  //find true minimum around this atom.
  for (float phi=0; phi<=180; phi+=3){
    float phi_rad=phi/180*PI;
    float sin_phi_rad=sin(phi_rad);
    float increm=PI/sin_phi_rad/60.0;
    for (float theta=0; phi!=0 && phi!=180 && theta<2*PI; theta+=increm){
      P.set(bl*sin_phi_rad*cos(theta), bl*sin_phi_rad*sin(theta),
	    bl*cos(phi_rad));
      P+=s->R;
      P.minimg(Lx,Ly);
      newatom->R=P;
      double uu=U_on_i(newatom);
      if (uu< u_min){
	*SV=P;
	u_min=uu;
      }
    }
  }
  // cerr<<"Findwell on "<<xx<<" returned "<<u_min<<" at "<<*SV<<endl;
  erase(newatom);
  Nmax--;
  TimeStepInit();
  return u_min;
}
