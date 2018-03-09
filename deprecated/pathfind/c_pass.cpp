#include "config.h"
#define SI_CL 2.1
#define SI_F 1.81
#define SI_SI 2.77
#define SI_F_sq 3.2761
#define SI_SI_sq 7.6729
#define VDW_SI 2.1
#define VDW_F 1.3
#define VDW_CL 1.65

//***********************************************************************
double config::U_on_i(VPI s, float max){
  if ((this_cell=WhichCell(s->R)) != s->my_cell){
        ((subcell*) s->my_cell)->erase(s);
        ((subcell*) this_cell)->insert(s);
  }
  subcell * scp=WhichCell(s->R);
  for (j0=s->nlist.begin(); j0!=s->nlist.end(); j0++){
    (*j0)->nlist.erase(s);
  }
  s->nlist.clear();
  for (c=scp->nbegin(); c!=scp->nend(); c++)
    for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
      if (s!=*j0)
	if ((s->R - (*j0)->R).minsqmag(Lx, Ly) < 
	    (s->rc >? (*j0)->rc)){
	  s->nlist.insert(*j0);
	  (*j0)->nlist.insert(s);
	}
  double d=0;
  for(j0=s->nlist.begin(); d<max && j0!=s->nlist.end(); j0++){
    d+=u_TwoBody(s,*j0)/2;
    if (s->three_body){
      trips_list.clear();
      for (trip=(*j0)->nlist.begin(); trip!=(*j0)->nlist.end(); trip++)
	if (*trip!=s) trips_list.insert(*trip);
      trip=j0;
      for (trip++; trip!=s->nlist.end(); trip++)
	trips_list.insert(*trip);
      for(trip=s->nlist.begin(); trip!=j0; trip++)
	trips_list.erase(*trip);
      for(k0=trips_list.begin(); k0!=trips_list.end(); k0++)
	d+=u_ThreeBody(s, *j0, *k0)/3;
    }
  }
  if (d<max) return d;
  else return max;
}
//**************************************************************************
double config::FindWell(short id, svector* P){
  VPI newatom=append(id);
  float res=.5; //angstroms
  int nnx=(int)(Lx/res);
  int nny=(int)(Ly/res);
  int nnz=(int)(Lz_full/res);
  float lcx=Lx/nnx;
  float lcy=Ly/nny;
  float lcz=Lz_full/nnz;
  float grid[nnx][nny][nnz];
  double u_min=1000;
  for (int z=0; z<nnz; z++)
    for (int x=0; x<nnx; x++)
      for (int y=0; y<nny; y++){
        newatom->R.set(Lx/2-x*lcx, Ly/2-y*lcy, Lz_full-Lz/2-z*lcz);
	newatom->R.minimg(Lx,Ly);
	newatom->u=grid[x][y][z]=U_on_i(newatom, 1);
        newatom->u=U_on_i(newatom,1);
	if (newatom->u < u_min){
	 u_min=newatom->u;
	 //*P=newatom->R;
	 *P=svector((float)x, (float)y, (float)z);
        }
      }
  fstream fout("ugrid.out",ios::out);
  char outbuf[256];
  sprintf(outbuf,"%i\t%i\t%i\n", 8, 8, 8);
  fout<<outbuf;
  for (int z=0; z<8; z++){
    for (int y=4; y<12; y++){
      for (int x=20; x<28; x++){
	sprintf(outbuf, "%7.4f ",grid[x][y][z]);
	fout<<outbuf;
      }
      fout<<endl;
    }
    fout<<endl;
  }
  cout<<*P<<endl;
  erase(newatom);
  Nmax--;
  return u_min;
}
//*************************************************************************
void config::Passivate(short id){
//    set<int>::iterator sii;
//    set<int> pass_list;
//    for (sii=surf_list.begin(); sii!=surf_list.end(); sii++){
//      VPI s=atomix(*sii);
//      if (s->nlist.size()<4){
//        pass_list.insert(*sii);
//      }
//      else{
//        int coord=0;
//        for (SPI spi=s->nlist.begin(); spi!=s->nlist.end(); spi++){
//  	if ((*spi)->id==14){
//  	  if (((*spi)->R - s->R).minsqmag(Lx,Ly) < SI_SI_sq)
//  	    coord++;
//  	}
//  	else if ((*spi)->id==9){
//  	  if (((*spi)->R - s->R).minsqmag(Lx,Ly) < SI_F_sq)
//  	    coord++;
//  	}
//        }
//        if (coord < 4)
//  	pass_list.insert(*sii);
//      }
//    }
  
//    if (!pass_list.empty()){
//      int kept=pass_list.size();
//      cerr<<"* Passivating. Found "<<kept<<" sites...";
//      float bl=0;
//      svector P;
//      if (id==9) bl=SI_F;
//      //else if (id==17) bl=SI_CL;
//      else if (id==14) bl=SI_SI;
//      else {cerr<<"\nI shouldn't passivate with this species!\n"; exit(1);}
    
//      for (sii=pass_list.begin(); sii!=pass_list.end(); sii++){
//        VPI newatom=append(id);
//        VPI target=atomix(*sii);
      
//        newatom->R=newatom->P;
//        for (double n=0; n<0.10; n+=dt){
//  	if (n>0){
//  	  newatom->PosUpdate(dt,dtsq2,Lx,Ly);
//  	  newatom->V+=dt2*newatom->F/newatom->m;
//  	  if ((this_cell=WhichCell(newatom->R)) != newatom->my_cell){
//  	    ((subcell*) newatom->my_cell)->erase(newatom);
//  	    ((subcell*) this_cell)->insert(newatom);
//  	  }
//  	}
//  	TimeStepInit();
//  	ReNeighbor();
//  	if (U_on_i(newatom) == 0) break;
//  	newatom->V+=dt2*newatom->F/newatom->m;
//        }
      
//        //if it is less than 110% of the bond length away, keep it.
//        if ((newatom->R - target->R).minmag(Lx,Ly)<bl*1.1){
//  	newatom->V.clear();
//  	newatom->O=newatom->R;
//        }
//        else{
//  	kept--;
//  	erase(newatom);
//  	Nmax--;
//        }
//      }
//    }
//    Partition(); TimeStepInit(); ReNeighbor();
//    cerr<<" Kept "<<kept<<". Total atoms: "<<N<<endl;
}


