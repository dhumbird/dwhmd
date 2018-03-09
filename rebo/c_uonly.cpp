#include "config.h"
extern map<short,double> Q, A, ALPHA, B1, B2, B3, BETA1, BETA2, BETA3;
//**********************************************************************
double config::Uonly(){
  double uu=0;
  for (i=begin; i<end; i++){
    for (bond_ij=i->nlist.begin(); bond_ij!=i->nlist.end(); bond_ij++){
      //establish bond pointers
      if ((j=(atom*)bond_ij->a2) > &(*i)){
	bond_ji=j->FindNbr(&(*i));
	scr_Rhat.clear();
	iNconj=jNconj=0;
	//calculate VA_ij, VR_ij;
	type=bond_ij->type;
	r=bond_ij->r;
	VA_ij=bond_ij->f*(B1[type]*exp(-BETA1[type]*r)
			  +B2[type]*exp(-BETA2[type]*r)
			  +B3[type]*exp(-BETA3[type]*r));
	
	VR_ij=bond_ij->f*(1+Q[type]/r)*A[type]*exp(-ALPHA[type]*r);
	double NH_ij = i->NH - (j->id==1)*bond_ij->f;
	double NC_ij = i->NC - (j->id==6)*bond_ij->f;
	double NH_ji = j->NH - (i->id==1)*bond_ij->f;
	double NC_ji = j->NC - (i->id==6)*bond_ij->f;
	double Nt_ij = NH_ij + NC_ij;
	double Nt_ji = NH_ji + NC_ji;
	//get P_ij, P_ji
	P_ij=P_ji=dHP_ij=dHP_ji=dCP_ij=dCP_ji=0;
	if (i->id==6){
	  if (j->id==1){ 
	    Pch_bicubicint(NH_ij, NC_ij, &P_ij, &dHP_ij,&dCP_ij);
	  }
	  if (j->id==6){
	    Pcc_bicubicint(NH_ij, NC_ij, &P_ij, &dHP_ij,&dCP_ij);
	    Pcc_bicubicint(NH_ji, NC_ji, &P_ji, &dHP_ji,&dCP_ji);
	  }
	}
	else {
	  if (j->id==6){
	    Pch_bicubicint(NH_ji, NC_ji, &P_ji, &dHP_ji,&dCP_ji);
	  }
	}
	bsp_ij=P_ij; bsp_ji=P_ji;
	
	if (Nt_ij < 3.2) S=1;
	else if (Nt_ij < 3.7){
	  xik=2*(Nt_ij-3.2);
	  S = 1 - pow(xik,3)*(6*xik*xik-15*xik+10);
	}
	else S=0;
	
	//***************for k neighbor of i*************************
	for (bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
	  if ((k=(atom*)bond_ik->a2)!=j){
	    if ((bond_jk=j->FindNbr(k))==j->nlist.end()){
	      svector R=(j->R - k->R); R.minimg(Lx,Ly);
	      rjk=R.mag();
	      R/=rjk;
	      scr_Rhat.push_back(R);
	    }
	    else rjk=bond_jk->r;
	    //*******b-sigma-pi calculations************
	    if (i->id == 1)
	      el=exp(4.0*((k->id==1)*RHH + (k->id==6)*RCH - bond_ik->r
			  -(j->id==1)*RHH - (j->id==6)*RCH + bond_ij->r));
	    else el = 1.0;
	    cosO = bond_ij->Rhat ^ bond_ik->Rhat;
	    Gpoly(cosO, i->id, &G, &Gprime, &gam, &gamprime);
	    if (gam!=0) g = G+S*(gam-G);
	    else g=G;
	    bsp_ij += g*bond_ik->f*el;
	    //*********conjugation*****************
	    if (k->id==6){
	      xik=k->Nt-bond_ik->f;
	      if (xik<2) iNconj+=bond_ik->f;
	      else if (xik<3) iNconj += bond_ik->f * 0.5*(1+cos(xik));
	    }
	  }
	}
	if (Nt_ji < 3.2) S=1;
	else if (Nt_ji < 3.7){
	  xik=2*(Nt_ji-3.2);
	  S = 1 - pow(xik,3)*(6*xik*xik-15*xik+10);
	}
	else S=0;
	//***************for k neighbor of j*************************
	for (bond_jk=j->nlist.begin(); bond_jk!=j->nlist.end(); bond_jk++){
	  if ((k=(atom*)bond_jk->a2)!=&(*i)){
	    if ((bond_ik=i->FindNbr(k))==i->nlist.end()){
	      svector R=(i->R - k->R); R.minimg(Lx,Ly);
	      rik=R.mag();
	      R/=rik;
	      scr_Rhat.push_back(R);
	    }
	    else rik=bond_ik->r;
	    //*********b-sigma-pi calculation*************
	    if (j->id == 1){
	      el=exp(4.0*((k->id==1)*RHH + (k->id==6)*RCH - bond_jk->r
			  -(i->id==1)*RHH - (i->id==6)*RCH + bond_ij->r));
	    }
	    else{
	      el = 1.0;
	    }
	    cosO = bond_ji->Rhat ^ bond_jk->Rhat;
	    Gpoly(cosO, j->id, &G, &Gprime, &gam, &gamprime);
	    if (gam!=0) g = G+S*(gam-G);
	    else g=G;
	    bsp_ji += g * bond_ij->f * el;
	    //*********conjugation************
	    if (k->id==6){
	      xik=k->Nt-bond_jk->f;
	      if (xik<2) jNconj+=bond_jk->f;
	      else if (xik<3) jNconj+=bond_jk->f*0.5*(1+cos(xik));
	    }
	  }
	}
	
	//finalize b-sigma-pi
	bsp_ij = 1.0 / sqrt(1 + bsp_ij);
	bsp_ji = 1.0 / sqrt(1 + bsp_ji);
	
	//***************conjugation************************
	Nconj_ij=1+iNconj*iNconj+jNconj*jNconj;
	if (type==0)
	  Fcc_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &F_ij, &dFi_ij, &dFj_ij, &dFc_ij);
	else if (type==2)
	  Fch_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &F_ij, &dFi_ij, &dFj_ij, &dFc_ij); 
	else if (type==1)
	  Fhh_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &F_ij, &dFi_ij, &dFj_ij, &dFc_ij); 
	else 
	  F_ij=dFi_ij=dFj_ij=dFc_ij=0;
	
	//***************dihedral*****************************
	bDH_ij=0;
	if (bond_ij->type==0){
	  Tcc_tricubicint(Nt_ij, Nt_ji, Nconj_ij,
			  &T_ij, &dTi_ij, &dTj_ij, &dTc_ij);
	  Rji=bond_ji->r*bond_ji->Rhat;
	  for (bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
	    if ((k=(atom*)bond_ik->a2) != j){
	      Rik=bond_ik->r*bond_ik->Rhat;
	      C_ikji = Rik * Rji;
	      c_ikji = C_ikji.mag();
	      for (bond_jl=j->nlist.begin(); bond_jl!=j->nlist.end(); 
		   bond_jl++){
		if ((l=(atom*)bond_jl->a2) != &(*i)){
		  
		  bond_lj=l->FindNbr(j);
		  Rlj=bond_lj->r*bond_lj->Rhat;
		  C_jilj = Rji * Rlj;
		  c_jilj = C_jilj.mag();
		  cosO = (1/c_jilj/c_ikji)*(C_ikji ^ C_jilj);

		  switch (k->id + l->id){
		  case 12: e_tors=eCCCC; break;
		  case 7:  e_tors=eHCCC; break;
		  case 2:  e_tors=eHCCH; break;
		  }
		  uu+= bond_ik->f * bond_ji->f * bond_jl->f * 
		    (tors_fac * e_tors * pow(cosO/2.0,10) - 0.1*e_tors);
		  bDH_ij+=(1-cosO*cosO)*bond_ik->w * bond_jl->w;
		}
	      }
	    }
	  }
	  bDH_ij *= T_ij;
	}
      	
	//compute bbar_ij
	bbar_ij=(bsp_ij + bsp_ji + F_ij + bDH_ij)/2;
	uu+=(VR_ij-bbar_ij*VA_ij);
      }
    }
  }
  return uu;
}
