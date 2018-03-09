#include "config.h"
//**************************************************************************
void config::Visit(atom* w, cluster* cl){
  VNI vj; atom* aj;
  cl->insert(w);
  if (w->is_fixed) cl->cl_is_fixed=1;
  for (vj=w->nlist.begin(); vj!=w->nlist.end(); vj++){
    aj=(atom*)vj->a2;
    if(aj->p_clust==NULL)
      Visit(aj, cl);
  }
}
//***************************************************************************
void config::CheckSput(){
  cfg_clusts.clear();
  for (VAI s=begin; s<end; s++)
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
// void config::BoundNeighbors(){
//   for (i=begin; i<end; i++) i->nlist.clear();
//   for (cells_it=cells.begin(); cells_it!=cells.end(); cells_it++)
//     for (i0=cells_it->abegin(); i0!=cells_it->aend(); i0++)  
//       for (c=cells_it->nbegin(); c!=cells_it->nend(); c++)
//         for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
//           if (*i0 < *j0){
// 	    float rc=0;
// 	    if ((*i0)->id + (*j0)->id==28) rc=7.673;
// 	    else if ((*i0)->id + (*j0)->id==23) rc=3.276;
// 	    else if ((*i0)->id + (*j0)->id==18) rc=2;
// 	    else rc=(*i0)->rc >? (*j0)->rc;
// 	    if (((*i0)->R - (*j0)->R).minsqmag(Lx, Ly) < rc){ 
// 	      (*i0)->nlist.insert(*j0);
// 	      (*j0)->nlist.insert(*i0);
// 	    }
// 	  }
// }
//**************************************************************************
void config::ClustSput(cluster* cl){
  //to project the candidate's path, normalize its velocity to
  //1 A/ps. Then use R+=V to move it 1 A along this line.
  svector sR,sV; double sm=0; float sr=0;
  int nC, nH;
  nC=nH=0;
  for (SAI s=cl->begin(); s!=cl->end(); s++){
    sR+=(*s)->R * (*s)->m;
    sr=(*(cl->begin()))->R.x - (*s)->R.x;
    if ((sr/Lx) > 0.5) sR.x += Lx * (*s)->m;
    if ((sr/Lx) < -0.5) sR.x -= Lx * (*s)->m;
    sr=(*(cl->begin()))->R.y - (*s)->R.y;
    if ((sr/Ly) > 0.5) sR.y += Ly * (*s)->m;
    if ((sr/Ly) < -0.5) sR.y -= Ly * (*s)->m;
    sV+=(*s)->m * (*s)->V;
    sm+=(*s)->m;
  }
  sV/=sV.mag();
  sR/=sm;
  sR.minimg(Lx, Ly);
  if (sR.z < -abs(Lz/2) && sV.z < 0)
    for (SAI s=cl->begin(); s!=cl->end(); s++)
      sput_list.insert((*s)->ix);
  else if (sV.z > 0){
    bool got_neighbors=0;
    while (sR.z < the_top && !got_neighbors){
      subcell* sc=WhichCell(sR);
      for (c=sc->nbegin(); c!=sc->nend(); c++)
	for (SAI sj=(*c)->abegin(); sj!=(*c)->aend(); sj++){
	  if (!got_neighbors)
	    if (!binary_search(cl->begin(), cl->end(), *sj))
	      if ((sR - (*sj)->R).minsqmag(Lx, Ly) < 4)
		got_neighbors=1;
	}
      sR+=sV; sR.minimg(Lx, Ly);
    }
    if (!got_neighbors){
      for (SAI s=cl->begin(); s!=cl->end(); s++){
	if ((*s)->id==6) nC++;
	if ((*s)->id==1) nH++;
	sput_list.insert((*s)->ix);
      }
      if (nC>0 || nH>0){
	cerr<<"* Cluster of ";
	if (nC>0){ 
	  cerr<<"C";
	  if (nC>1) cerr<<nC;
	}
	if (nH>0){
	  cerr<<"H";
	  if (nH>1) cerr<<nH;
	}
	cerr<<" sputtered"<<endl;
      }
    }
  }
}
//**************************************************************************
void config::Sputter(int w){
  VAI s=atomix(w);
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
// void config::Sweep(){
//   cerr<<"* Sweeping...";
//   BoundNeighbors();
//   //encluster
//   cfg_clusts.clear();
//   for (VPI s=begin; s<end; s++)
//     if (s->p_clust==NULL){
//       cfg_clusts.resize(cfg_clusts.size()+1);
//       Visit(s,cfg_clusts.end()-1);
//     }
//   if (cfg_clusts.size()>1){
//     vector<cluster>::iterator cl_it;
//     cerr<<"Found "<<cfg_clusts.size()-1;
//     for (cl_it=cfg_clusts.begin(); cl_it!=cfg_clusts.end(); cl_it++)
//       if (!cl_it->cl_is_fixed)
// 	ClustSput(cl_it);
	
//   }
//   if (!sput_list.empty()){
//     set<int>::iterator sput;
//     for (sput=sput_list.begin(); sput!=sput_list.end(); sput++)
//       Sputter(*sput);
//     sput_list.clear();
//     Partition();
//   }
//   cerr<<endl;
//   ReNeighbor();
// }





