#include "config.h"
extern map<short,double> RC_SQ;
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
    if (s->R.z < the_bottom) sput_list.insert(s->ix);
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
    //sputter backwards so I save as many indices as possible
    set<int>::reverse_iterator sput;
    for (sput=sput_list.rbegin(); sput!=sput_list.rend(); sput++)
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
  svector sR,sV; double rc;
  double rmin=1000;
  //center-of-mass velocity
  for (SAI s=cl->begin(); s!=cl->end(); s++){
    sV+=(*s)->m * (*s)->V;
    rmin=mfmin((*s)->R.z,rmin);
  }
  sV/=sV.mag();
  svector Rtemp;
  if (sV.z > 0){
    bool got_neighbors=0;
    while (rmin < the_top && !got_neighbors){
      for (SAI i=cl->begin(); i!=cl->end(); i++){
	for (VAI j=begin; j<end; j++){
	  if (!got_neighbors)
	    if (!binary_search(cl->begin(), cl->end(), &(*j))){
	      Rtemp=(*i)->R + sR; Rtemp.minimg(Lx,Ly);
	      if ((*i)->ix==inert) rc=36;
	      else rc=RC_SQ[(*i)->id+j->id];
	      if ((Rtemp - j->R).minsqmag(Lx, Ly) < rc)
		got_neighbors=1;
	    }
	}
      }
      sR+=sV;
      rmin+=sV.z;
    }
    if (!got_neighbors){
      FillSputList(cl,0);
    }
  }
}
//*********************************************************************
string config::ClusterName(cluster* cl){
  int nC, nF, nSi, nCl;
  nC=nF=nSi=nCl=0;
  for (SAI s=cl->begin(); s!=cl->end(); s++){
    if ((*s)->id==6) nC++;
    if ((*s)->id==9) nF++;
    if ((*s)->id==14) nSi++;
    if ((*s)->id==17) nCl++;
  }
  ostringstream name;
  if (nC>0 || nF>0 || nSi>0 || nCl>0){
    if (nSi>0){
      name<<"Si";
      if (nSi>1) name<<nSi;
    }
    if (nC>0){ 
      name<<"C";
      if (nC>1) name<<nC;
    }
    if (nF>0){
      name<<"F";
      if (nF>1) name<<nF;
    }
    if (nCl>0){
      name<<"Cl";
      if (nCl>1) name<<nCl;
    }
  }
  return name.str();
}
//**********************************************************************
void config::FillSputList(cluster* cl, bool sweep){
  double cek=0;
  svector sV; double mass=0;
  //center-of-mass velocity
  for (SAI s=cl->begin(); s!=cl->end(); s++){
    //cek+=(*s)->Ek();
    sV+=(*s)->m * (*s)->V;
    mass+=(*s)->m;
    sput_list.insert((*s)->ix);
  }
  sV/=mass;
  cek=0.5*mass*sV.sqmag();
  cerr<<"* Cluster of "<<ClusterName(cl);
  if (sweep) cerr<<" swept."<<endl;
  else cerr<<" sputtered. "<<sV<<" "<<cek<<" eV."<<endl;
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
    if (s->R.z < the_bottom)
      cerr<<"* Atom "<<w<<" ("<<s->id<<") pushed through bottom. Time: "<<t<<endl;
    else
      cerr<<"* Atom "<<w<<" ("<<s->id<<") sputtered from "
	  <<((s->P.is_zero())?s->R:s->P)<<" Time: "<<t<<endl;
  }
  erase(s);
}
//***************************************************************************
void config::BoundVisit(atom* w, cluster* cl){
  VNI vj; atom* aj;
  cl->insert(w);
  if (w->is_fixed) cl->cl_is_fixed=1;
  for (vj=w->nlist.begin(); vj!=w->nlist.end(); vj++){
    if (vj->f == 1){
      aj=(atom*)vj->a2;
      if(aj->p_clust==NULL)
	BoundVisit(aj, cl);
    }
  }
}
//***********************************************************************
void config::Sweep(){
  //cerr<<"* Sweeping...";
  cfg_clusts.clear();
  for (VAI s=begin; s<end; s++)
    if (s->p_clust==NULL){
      cfg_clusts.resize(cfg_clusts.size()+1);
      BoundVisit(&(*s),&(*(cfg_clusts.end()-1)));
    }
  if (cfg_clusts.size()>1){
    vector<cluster>::iterator cl_it;
    //cerr<<"Found "<<cfg_clusts.size()-1<<endl;
    for (cl_it=cfg_clusts.begin(); cl_it!=cfg_clusts.end(); cl_it++){
      if (!cl_it->cl_is_fixed){
	for (SAI s=cl_it->begin(); s!=cl_it->end(); s++){
	  for (VNI vj=(*s)->nlist.begin(); vj!=(*s)->nlist.end(); vj++){
	    atom* aj=(atom*)vj->a2;
	    if (!binary_search(cl_it->begin(), cl_it->end(), aj))
	      cl_it->fmax = mfmax(vj->f, cl_it->fmax);
	  }
	}
	if (AutoSweep(&(*cl_it)) || cl_it->fmax < 0.10){
//  	  for (SAI s=cl_it->begin(); s!=cl_it->end(); s++)
//   	    cerr<<(*s)->ix<<" "<<(*s)->R<<" ";
//   	  cerr<<cl_it->fmax<<endl;
	  svector Rcom;
	  float cl_rad=0; float cl_mass=0;
	  SAI anchor=cl_it->begin();
	  for (SAI s=cl_it->begin(); s!=cl_it->end(); s++){
	    for (SAI s2=cl_it->begin(); s2!=cl_it->end(); s2++){
	      float rad=((*s)->R - (*s2)->R).minsqmag(Lx,Ly);
	      cl_rad=mfmax(cl_rad, rad);
	    }
	    if (s==anchor)
 	      Rcom=(*s)->R*(*s)->m;
 	    else{
 	      float tr=(*s)->R.x - (*anchor)->R.x;
	      if (fabs(tr)<Lx/2)
		Rcom.x+=(*s)->R.x*(*s)->m;
	      else{
		if (tr<0) Rcom.x+=((*s)->R.x+Lx)*(*s)->m;
		else Rcom.x += ((*s)->R.x-Lx)*(*s)->m;
	      }
	      tr=(*s)->R.y - (*anchor)->R.y;
	      if (fabs(tr)<Ly/2)
		Rcom.y+=(*s)->R.y*(*s)->m;
	      else{
		if (tr<0) Rcom.y+=((*s)->R.y+Ly)*(*s)->m;
		else Rcom.y += ((*s)->R.y-Ly)*(*s)->m;
	      }
	      Rcom.z+=(*s)->R.z*(*s)->m;
	    }
 	    cl_mass+=(*s)->m;
	  }
	  Rcom/=cl_mass;
	  Rcom.minimg(Lx,Ly);
	  cl_rad=sqrt(cl_rad)/2;
	  cl_rad*=cl_rad;
	  //line of sight path?
	  bool flag=1;
	  if (the_top - Rcom.z < 15){
	    float sin_phi_rad,phi_rad,increm; svector P, Rtemp;
	    for (float phi=0; flag && phi<=90; phi+=3){
	      phi_rad=phi/180*PI;
	      sin_phi_rad=sin(phi_rad);
	      increm=PI/sin_phi_rad/90.0;
	      for (float theta=0; flag && phi!=0 && phi!=180 && theta<2*PI; 
		   theta+=increm){
		P.set(1, theta, phi_rad);
		P.Sph2Rec();
		Rtemp=Rcom;
		bool got_neighbors=0;
		while (Rtemp.z < the_top && !got_neighbors){
		  for (SAI s=cl_it->begin(); s!=cl_it->end(); s++){
		    for (VAI vj=begin; vj<end; vj++){
		      if (!got_neighbors)
			if (!binary_search(cl_it->begin(), cl_it->end(), &(*vj))){
			  if ((Rtemp - vj->R).minsqmag(Lx, Ly) < cl_rad){
			    got_neighbors=1;
			    vj=end;
			    s=cl_it->end();
			  }
			}
		    }
		  }
		  Rtemp+=P;
		}
		if (!got_neighbors) flag=0;
	      }
	    }
	    if (!flag) FillSputList(&(*cl_it), 1);
	  }
	}
      }
    }
    if (!sput_list.empty()){
      set<int>::iterator sput;
      for (sput=sput_list.begin(); sput!=sput_list.end(); sput++)
	Sputter(*sput);
      sput_list.clear();
      Partition();
    }
  }
  //else cerr<<"done"<<endl;
  ReNeighbor();
}
//******************************************************************
bool config::AutoSweep(cluster* cl){
  string name=ClusterName(cl);
  if (name=="SiF4") return 1;
  else if (name=="Si2F6") return 1;
  else if (name=="Si3F8") return 1;
  else if (name=="SiCl4") return 1;
  else if (name=="Si2Cl6") return 1;
  else if (name=="Si3Cl8") return 1;
  else if (name=="CF4") return 1;
  else if (name=="C2F6") return 1;
  else return 0;
}
