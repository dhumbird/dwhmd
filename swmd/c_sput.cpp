#include "config.h"
//**************************************************************************
void config::Visit(particle* w, cluster* cl){
  SPI spi_j;
  cl->insert(w);
  if (w->is_fixed) cl->cl_is_fixed=1;
  for (spi_j=w->nlist.begin(); spi_j!=w->nlist.end(); spi_j++)
    if((*spi_j)->p_clust==NULL)
      Visit(*spi_j, cl);
}
//***************************************************************************
void config::CheckSput(){
  cfg_clusts.clear();
  for (VPI s=begin; s<end; s++)
    if (s->p_clust==NULL){
      cfg_clusts.resize(cfg_clusts.size()+1);
      Visit(&(*s),&(*(cfg_clusts.end()-1)));
    }
  
  if (cfg_clusts.size()>1){
    vector<cluster>::iterator cl_it;
    for (cl_it=cfg_clusts.begin(); cl_it!=cfg_clusts.end(); cl_it++)
      if (!cl_it->cl_is_fixed)
	ClustSput(&(*cl_it));
  }
  if (!sput_list.empty()){
    set<int>::iterator sput;
    for (sput=sput_list.begin(); sput!=sput_list.end(); sput++)
      Sputter(*sput);
    sput_list.clear();
    Partition();
    ReNeighbor();
  }
}
//************************************************************************
void config::BoundNeighbors(){
  for (i=begin; i<end; i++) i->nlist.clear();
  for (cells_it=cells.begin(); cells_it!=cells.end(); cells_it++)
    for (i0=cells_it->abegin(); i0!=cells_it->aend(); i0++)  
      for (c=cells_it->nbegin(); c!=cells_it->nend(); c++)
        for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
          if (*i0 < *j0){
	    float rc=0;
	    if ((*i0)->id + (*j0)->id==28) rc=7.673;
	    else if ((*i0)->id + (*j0)->id==23) rc=3.276;
	    else if ((*i0)->id + (*j0)->id==18) rc=2;
	    else rc=(*i0)->rc >? (*j0)->rc;
	    if (((*i0)->R - (*j0)->R).minsqmag(Lx, Ly) < rc){ 
	      (*i0)->nlist.insert(*j0);
	      (*j0)->nlist.insert(*i0);
	    }
	  }
}
//**************************************************************************
void config::ClustSput(cluster* cl){
  //to project the candidate's path, normalize its velocity to
  //1 A/ps. Then use R+=V to move it 1 A along this line.
  svector R,V; double m=0; float r=0;
  int nSi, nF;
  nSi=nF=0;
  for (SPI s=cl->begin(); s!=cl->end(); s++){
    R+=(*s)->R * (*s)->m;
    r=(*(cl->begin()))->R.x - (*s)->R.x;
    if ((r/Lx) > 0.5) R.x += Lx * (*s)->m;
    if ((r/Lx) < -0.5) R.x -= Lx * (*s)->m;
    r=(*(cl->begin()))->R.y - (*s)->R.y;
    if ((r/Ly) > 0.5) R.y += Ly * (*s)->m;
    if ((r/Ly) < -0.5) R.y -= Ly * (*s)->m;
    V+=(*s)->m * (*s)->V;
    m+=(*s)->m;
  }
  V/=V.mag();
  R/=m;
  R.minimg(Lx, Ly);
  if (R.z < -abs(Lz/2) && V.z < 0)
    for (SPI s=cl->begin(); s!=cl->end(); s++)
      sput_list.insert((*s)->ix);
  else if (V.z > 0){
    bool got_neighbors=0;
    while (R.z < the_top && !got_neighbors){
      subcell* sc=WhichCell(R);
      for (c=sc->nbegin(); c!=sc->nend(); c++)
	for (SPI sj=(*c)->abegin(); sj!=(*c)->aend(); sj++){
	  if (!got_neighbors)
	    if (!binary_search(cl->begin(), cl->end(), *sj))
	      if ((R - (*sj)->R).minsqmag(Lx, Ly) < (*sj)->rc)
		got_neighbors=1;
	}
      R+=V; R.minimg(Lx, Ly);
    }
    if (!got_neighbors){
      for (SPI s=cl->begin(); s!=cl->end(); s++){
	if ((*s)->id==14) nSi++;
	if ((*s)->id==9) nF++;
	sput_list.insert((*s)->ix);
      }
      if (nSi>0 || nF>0){
	cerr<<"* Cluster of ";
	if (nSi>0){ 
	  cerr<<"Si";
	  if (nSi>1) cerr<<nSi;
	}
	if (nF>0){
	  cerr<<"F";
	  if (nF>1) cerr<<nF;
	}
	cerr<<" sputtered"<<endl;
      }
    }
  }
}
//**************************************************************************
void config::Sputter(int w){
  VPI s=atomix(w);
  if (w==inert){
    cerr<<"* Ion scattered. Time: "<<t<<". Exit V: "<<s->V
	<<"("<<s->Ek()<<" eV)"<<endl;
    inert=-1;
  }
  else{
    if (s->R.z < -abs(Lz/2))
      cerr<<"* Atom "<<w<<" ("<<s->id<<") pushed through bottom. Time: "<<t<<endl;
    else
      cerr<<"* Atom "<<w<<" ("<<s->id<<") sputtered from "
	  <<((s->P.is_zero())?s->R:s->P)<<" Time: "<<t<<endl;
  }
  erase(s);
}
//***************************************************************************
void config::Sweep(){
  cerr<<"* Sweeping...";
  BoundNeighbors();
  //encluster
  cfg_clusts.clear();
  for (VPI s=begin; s<end; s++)
    if (s->p_clust==NULL){
      cfg_clusts.resize(cfg_clusts.size()+1);
      Visit(&(*s),&(*(cfg_clusts.end()-1)));
    }
  if (cfg_clusts.size()>1){
    vector<cluster>::iterator cl_it;
    cerr<<"Found "<<cfg_clusts.size()-1;
    for (cl_it=cfg_clusts.begin(); cl_it!=cfg_clusts.end(); cl_it++)
      if (!cl_it->cl_is_fixed)
	ClustSput(&(*cl_it));
	
  }
  if (!sput_list.empty()){
    set<int>::iterator sput;
    for (sput=sput_list.begin(); sput!=sput_list.end(); sput++)
      Sputter(*sput);
    sput_list.clear();
    Partition();
  }
  cerr<<endl;
  ReNeighbor();
}





