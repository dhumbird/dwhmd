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
  for (i=begin; i<end; i++)
    if (i->p_clust==NULL){
      cfg_clusts.resize(cfg_clusts.size()+1);
      Visit(i,cfg_clusts.end()-1);
    }
  
  if (cfg_clusts.size()>1)
    for (cl_it=cfg_clusts.begin(); cl_it!=cfg_clusts.end(); cl_it++)
      if (!cl_it->cl_is_fixed)
	ClustSput(cl_it);

  if (!sput_list.empty()){
    for (sput=sput_list.begin(); sput!=sput_list.end(); sput++)
      Sputter(*sput);
    sput_list.clear();
    Partition();
    ReNeighbor();
  }
}
//**************************************************************************
void config::ClustSput(cluster* cl){
  //to project the candidate's path, normalize its velocity to
  //1 A/ps. Then use R+=V to move it 1 A along this line.
  svector R,V; double m=0; float r=0;
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
	for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++){
	  if (!got_neighbors)
	    if (!binary_search(cl->begin(), cl->end(), *j0))
	      if ((R - (*j0)->R).minsqmag(Lx, Ly) < (*j0)->rc)
		got_neighbors=1;
	}
      R+=V; R.minimg(Lx, Ly);
    }
    if (!got_neighbors)
      for (SPI s=cl->begin(); s!=cl->end(); s++)
	sput_list.insert((*s)->ix);
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
bool config::Sweep(){
  cerr<<"* Sweeping...";
  
  //encluster
  cfg_clusts.clear();
  for (i=begin; i<end; i++)
    if (i->p_clust==NULL){
      cfg_clusts.resize(cfg_clusts.size()+1);
      Visit(i,cfg_clusts.end()-1);
    }
  
  if (cfg_clusts.size()>1){
    cerr<<"Found "<<cfg_clusts.size()-1<<"...";
    for (cl_it=cfg_clusts.begin(); cl_it!=cfg_clusts.end(); cl_it++){
      if (!cl_it->cl_is_fixed){
	svector R,Ro,V;
	double m=0; float r=0;
	//check for that goddamn inert
	for (SPI s=cl_it->begin(); s!=cl_it->end(); s++)
	  if ((*s)->ix==inert) Sputter(inert);
	//construct center-of-mass position and velocity
	for (SPI s=cl_it->begin(); s!=cl_it->end(); s++){
	  R+=(*s)->m * (*s)->R;
	  r=(*(cl_it->begin()))->R.x - (*s)->R.x;
	  if ((r/Lx) > 0.5) R.x += Lx * (*s)->m;
	  if ((r/Lx) < -0.5) R.x -= Lx * (*s)->m;
	  r=(*(cl_it->begin()))->R.y - (*s)->R.y;
	  if ((r/Ly) > 0.5) R.x += Ly * (*s)->m;
	  if ((r/Ly) < -0.5) R.x -= Ly * (*s)->m;
	  V+=(*s)->m * (*s)->V;
	  m+=(*s)->m;
	}
	V/=V.mag();
	R/=m;
	R.minimg(Lx, Ly);
	Ro=R;
	bool got_neighbors=0;
	while (!got_neighbors){
	  subcell* sc=WhichCell(R);
	  for (c=sc->nbegin(); c!=sc->nend(); c++)
	    for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	      if (!got_neighbors)
		if (!binary_search(cl_it->begin(), cl_it->end(), *j0))
		  if ((R - (*j0)->R).minsqmag(Lx, Ly) < (*j0)->rc)
		    got_neighbors=1;
	  R+=V; R.minimg(Lx,Ly);
	}
	//put the cluster there and equilibrate.
	for (SPI s=cl_it->begin(); s!=cl_it->end(); s++){
	  (*s)->R+=(R-Ro);
	  (*s)->R.minimg(Lx,Ly);
	}
	Run(0.05, 0, 0, 0, 0, -1, 1, 0,0);
	return 0;
      }
    }
  }
  cerr<<"done."<<endl;
  return 1;
}





