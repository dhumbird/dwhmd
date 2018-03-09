#include "config.h"
extern map<short,double> Q, A, ALPHA, B1, B2, B3, BETA1, BETA2, BETA3;
double fij, f1ij, tt;
//**********************************************************************
void config::ForceEval(){
  for (i=begin; i<end; i++){
    if (i->ix == inert){
      for (bond_ij=i->nlist.begin(); bond_ij!=i->nlist.end(); bond_ij++){
	j=(atom*)bond_ij->a2;
	r=bond_ij->r;
	if (j->id == 6){
	  moa=MO_AR_C_a;
	  moeps=MO_AR_C_EPS;
	}
	else if (j->id == 1){
	  moa=MO_AR_H_a;
	  moeps=MO_AR_H_EPS;
	}
	rija=r/moa;
	Fij=moeps/r*(0.35*SpExp(-0.3*rija)*(1/r+0.3/moa)
		     +0.55*SpExp(-1.2*rija)*(1/r+1.2/moa)
		     +0.1*SpExp(-6*rija)*(1/r+6/moa))*bond_ij->Rhat;
	i->F += Fij; j->F -= Fij;
	u+= moeps/r*(0.35*SpExp(-0.3*rija)
		     +0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
      }
    }
    else{
      for (bond_ij=i->nlist.begin(); bond_ij!=i->nlist.end(); bond_ij++){
	//establish bond pointers
	if ((j=(atom*)bond_ij->a2) > &(*i)){
	  bond_ji=j->FindNbr(&(*i));
	  //clear storage
	  iFv_ij=0; iFv_ik.clear(); iFv_jk.clear(); iFv_ik2.clear();
	  iFv_kl.clear();
	  jFv_ij=0; jFv_ik.clear(); jFv_jk.clear(); jFv_jk2.clear();
	  jFv_kl.clear();
	  scr_Rhat.clear();
	  sumdg_i=sumdg_j=0;
	  iNconj=jNconj=0;
	  //calculate VA_ij, VR_ij;
	  type=bond_ij->type;
	  r=bond_ij->r;
	  fij=bond_ij->f; f1ij=bond_ij->fprime;
	  
	  b1=B1[type]; b2=B2[type]; b3=B3[type];
	  beta1=BETA1[type]; beta2=BETA2[type]; beta3=BETA3[type];
	  VA_ij=b1*exp(-beta1*r) + b2*exp(-beta2*r) + b3*exp(-beta3*r);
	  dVA_ij=b1*beta1*exp(-beta1*r) + b2*beta2*exp(-beta2*r)
	    +b3*beta3*exp(-beta3*r);
	  dVA_ij=VA_ij*f1ij-fij*dVA_ij;
	  VA_ij*=fij;
	  qq=Q[type]/r; aa=A[type]*SpExp(-ALPHA[type]*r);
	  VR_ij=fij*(1+qq)*aa;
	  dVR_ij=aa*(f1ij*(1+qq)-fij*(qq/r+ALPHA[type]*(1+qq)));
   
	  double NH_ij = i->NH - (j->id==1)*bond_ij->f;
	  double NC_ij = i->NC - (j->id==6)*bond_ij->f;
	  double NH_ji = j->NH - (i->id==1)*bond_ij->f;
	  double NC_ji = j->NC - (i->id==6)*bond_ij->f;
	  double Nt_ij = NH_ij + NC_ij;
	  double Nt_ji = NH_ji + NC_ji;

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

	  PolySwitch(2*(Nt_ij-3.2), &S, &Sprime);
	  Sprime*=2;
 	  //***************for k neighbor of i*************************
	  for (bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
	    if ((k=(atom*)bond_ik->a2)!=j){
	      if ((bond_jk=j->FindNbr(k))==j->nlist.end()){
		svector R=(j->R - k->R); R.minimg(Lx,Ly);
		rjk=R.mag();
		R/=rjk;
		scr_Rhat.push_back(R);
	      }
	      else
		rjk=bond_jk->r;
	      //*******b-sigma-pi calculations************
	      if (i->id == 1){
		el=exp(4.0*((k->id==1)*RHH + (k->id==6)*RCH - bond_ik->r
			    -(j->id==1)*RHH - (j->id==6)*RCH + r));
		dlam1=4.0; dlam2=-4.0; 
		elf=bond_ik->f*el;
	      }
	      else{
		el = 1.0; dlam1=dlam2=0;
		elf=bond_ik->f;
	      }
	      cosO = bond_ij->Rhat ^ bond_ik->Rhat;
	      Gpoly(cosO, i->id, &G, &Gprime, &gam, &gamprime);
  
	      if (gam!=0){
		g = G+S*(gam-G);
		g1=Gprime+S*(gamprime-Gprime);
		sumdg_i+=elf*(gam-G);
	      }
	      else{
		g=G; g1=Gprime;
	      }
	      
	      iFv_ij+=(elf*(g1*(1/bond_ik->r - cosO/r) + g*dlam1));
	      iFv_ik.push_back(bond_ik->fprime*
			       (g*el +(k->id==6)*dCP_ij + (k->id==1)*dHP_ij)
			       +elf*(g1*(1/r - cosO/bond_ik->r) + g*dlam2));
	      iFv_jk.push_back(elf*(-g1*rjk/bond_ik->r/r));
	      bsp_ij += g*elf;
	      //*********conjugation*****************
	      if (k->id==6){
		PolySwitch(k->Nt - bond_ik->f - 2, &Yik, &Y1ik);
		iNconj += bond_ik->f * Yik;
		iFv_ik2.push_back(bond_ik->fprime * Yik);
		iFv_kl.push_back(bond_ik->f*Y1ik);
	      }
	      else{
		iFv_ik2.push_back(0); iFv_kl.push_back(0);
	      }
	    }
	  }
	  sumdg_i*=Sprime;
	  
	  PolySwitch(2*(Nt_ji-3.2), &S, &Sprime);
	  Sprime*=2;

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
			    -(i->id==1)*RHH - (i->id==6)*RCH + r));
		dlam1=4.0; dlam2=-4.0; 
		elf=bond_jk->f*el;
	      }
	      else{
		el = 1.0; dlam1=dlam2=0;
		elf=bond_jk->f;
	      }
	      cosO = bond_ji->Rhat ^ bond_jk->Rhat;
	      Gpoly(cosO, j->id, &G, &Gprime, &gam, &gamprime);
	      
	      if (gam!=0){
		g = G+S*(gam-G);
		g1=Gprime+S*(gamprime-Gprime);
		sumdg_j+=elf*(gam-G);
	      }
	      else{
		g=G; g1=Gprime;
	      }
	      
	      jFv_ij+=(elf*(g1*(1/bond_jk->r - cosO/r) + g*dlam1));
	      jFv_ik.push_back(elf*(-g1*rik/bond_jk->r/r));
	      jFv_jk.push_back(bond_jk->fprime
			       *(g*el + (k->id==6)*dCP_ji + (k->id==1)*dHP_ji)
			       +elf*(g1*(1/r - cosO/bond_jk->r) + g*dlam2));
	      bsp_ji += g * elf;
	      //*********conjugation************
	      if (k->id==6){
		PolySwitch(k->Nt - bond_jk->f - 2, &Yik, &Y1ik);
		jNconj += bond_jk->f * Yik;
		jFv_jk2.push_back(bond_jk->fprime * Yik);
		jFv_kl.push_back(bond_jk->f*Y1ik);
	      }
	      else{
		jFv_jk2.push_back(0); jFv_kl.push_back(0);
	      }
	    }
	  }
	  sumdg_j*=Sprime;

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
	  if (type==0){
	    Tcc_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &T_ij, &dTi_ij,
			    &dTj_ij, &dTc_ij);
	    Sprime=0.5*VA_ij*T_ij;
	    Rji=r*bond_ji->Rhat;
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
		    
		    Fij=C_ikji*Rji; //used as temp. vars
		    R=C_jilj*Rji;
		      
		    dcosO_ik = cosO/c_ikji/c_ikji*Fij-1/c_jilj/c_ikji*R;
		    
		    dcosO_ji = cosO/c_jilj/c_jilj*(C_jilj*Rlj)
		      - cosO/c_ikji/c_ikji*(C_ikji*Rik)
		      - 1/c_ikji/c_jilj*(C_ikji * Rlj + Rik * C_jilj);
		    
		    dcosO_lj = -cosO/c_jilj/c_jilj*R 
		      + 1/c_ikji/c_jilj*Fij;

		    //***AIREBO torsional potential
		    switch (k->id + l->id){
		    case 12: e_tors=eCCCC; break;
		    case 7:  e_tors=eHCCC; break;
		    case 2:  e_tors=eHCCH; break;
		    }
		    v_tors = tors_fac * e_tors
		      * pow(cosO/2.0,10) - 0.1*e_tors;
		    dv_tors = tors_fac * e_tors * 5 * pow(cosO/2.0,9)
		      *(bond_ik->f*bond_ji->f*bond_jl->f);
		    u+= bond_ik->f * bond_ji->f * bond_jl->f * v_tors;
		    
		    k->F+=dcosO_ik*dv_tors;
		    i->F-=dcosO_ik*dv_tors;
		    i->F+=dcosO_ji*dv_tors;
		    j->F-=dcosO_ji*dv_tors;
		    l->F-=dcosO_lj*dv_tors;
		    j->F+=dcosO_lj*dv_tors;
		    
		    if (bond_ik->fprime){
		      Fij=(v_tors*bond_jl->f*bond_ji->f*bond_ik->fprime)
			*bond_ik->Rhat;
		      i->F-=Fij; k->F+=Fij;
		    }
		    if (bond_ji->fprime){
		      Fij=(v_tors*bond_jl->f*bond_ik->f*bond_ji->fprime)
			*bond_ji->Rhat;
		      i->F+=Fij; j->F-=Fij;
		    }
		    if (bond_jl->fprime){
		      Fij=(v_tors*bond_ji->f*bond_ik->f*bond_jl->fprime)
			*bond_jl->Rhat;
		      j->F-=Fij; l->F+=Fij;
		    }
		    //***REBO standard torsional
		    if (T_ij){
		      S=cosO*2*Sprime*bond_jl->w*bond_ik->w;
		      k->F+=dcosO_ik*S;
		      i->F-=dcosO_ik*S;
		      i->F+=dcosO_ji*S;
		      j->F-=dcosO_ji*S;
		      l->F-=dcosO_lj*S;
		      j->F+=dcosO_lj*S;
		      
		      S = (1-cosO*cosO);
		      bDH_ij+=S*bond_ik->w * bond_jl->w;
		      if (bond_jl->wprime){
			Fij=Sprime*bond_ik->w*S*bond_jl->wprime
			  *bond_jl->Rhat;
			l->F-=Fij; j->F+=Fij;
		      }
		      if (bond_ik->wprime){
			Fij=Sprime*bond_jl->w*S*bond_ik->wprime
			  *bond_ik->Rhat;
			k->F-=Fij; i->F+=Fij;
		      }
		    }
		  }
		}
	      }
	    }
	    //add terms from the conjugation; they will get handled later.
	    if (T_ij){
	      dFi_ij += dTi_ij * bDH_ij;
	      dFj_ij += dTj_ij * bDH_ij;
	      dFc_ij += dTc_ij * bDH_ij;
	      bDH_ij *= T_ij;
	    }
	  }
	  
	  //compute bbar_ij
	  bbar_ij=(bsp_ij + bsp_ji + F_ij + bDH_ij)/2;
	  u+=(VR_ij-bbar_ij*VA_ij);
	  
	  cnt=0;
	  scrcount=0;

	  bsp_ij=-0.25*pow(bsp_ij,3)*VA_ij; //bsp becomes temp variable
	  bsp_ji=-0.25*pow(bsp_ji,3)*VA_ij;

	  Fij=(-dVR_ij+bbar_ij*dVA_ij+bsp_ij*iFv_ij+bsp_ji*jFv_ij)
	    *bond_ij->Rhat;
	  i->F+=Fij; j->F-=Fij;

	  for (bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
	    if ((k=(atom*)bond_ik->a2)!=j){
	      if ((bond_jk=j->FindNbr(k))==j->nlist.end())
		R=scr_Rhat[scrcount++];
	      else R=bond_jk->Rhat;
	      force_ik=bsp_ij*iFv_ik[cnt];
	      force_jk=bsp_ij*iFv_jk[cnt];
	      if (bond_ik->fprime){
		force_ik+=bsp_ij*bond_ik->fprime*sumdg_i;
		if (dFi_ij) force_ik+=0.5*VA_ij*bond_ik->fprime*dFi_ij;
	      }
	      if (k->id==6 && dFc_ij){
		Sprime = VA_ij * dFc_ij * iNconj ;
		if (iFv_ik2[cnt]) force_ik+=Sprime*iFv_ik2[cnt];
		if (iFv_kl[cnt]){
		  Sprime*=iFv_kl[cnt];
		  for (bond_kl=k->nlist.begin(); bond_kl!=k->nlist.end();
		       bond_kl++)
		    if (bond_kl->fprime && (l=(atom*)bond_kl->a2)!=&(*i)){
		      Fij=(Sprime*bond_kl->fprime)*bond_kl->Rhat;
		      k->F += Fij; l->F -= Fij;
		    }
		}
	      }
	      Fij = force_ik*bond_ik->Rhat;
	      i->F += Fij; k->F -= Fij;
	      Fij = force_jk*R;
	      j->F += Fij; k->F -= Fij;
	      cnt++;
	    }
	  }
	  cnt=0;
	  for (bond_jk=j->nlist.begin(); bond_jk!=j->nlist.end(); bond_jk++){
	    if ((k=(atom*)bond_jk->a2)!=&(*i)){
	      if ((bond_ik=i->FindNbr(k))==i->nlist.end())
		R=scr_Rhat[scrcount++];
	      else R=bond_ik->Rhat;
	      force_jk=bsp_ji*jFv_jk[cnt];
	      force_ik=bsp_ji*jFv_ik[cnt];
	      if (bond_jk->fprime){
		force_jk+=bsp_ji*bond_jk->fprime*sumdg_j;
		if (dFj_ij) force_jk+=0.5*VA_ij*bond_jk->fprime*dFj_ij;
	      }
	      if (k->id==6 && dFc_ij){
		Sprime=VA_ij*dFc_ij*jNconj;
		if (jFv_jk2[cnt]) force_jk+=Sprime*jFv_jk2[cnt];
		if (jFv_kl[cnt]){
		  Sprime*=jFv_kl[cnt];
		  for (bond_kl=k->nlist.begin(); bond_kl!=k->nlist.end();
		       bond_kl++)
		    if (bond_kl->fprime && (l=(atom*)bond_kl->a2)!=j){
		      Fij = (Sprime*bond_kl->fprime)*bond_kl->Rhat;
		      k->F += Fij; l->F -= Fij;
		    }
		}
	      }
	      Fij = force_ik*R;
	      i->F += Fij; k->F -= Fij;
	      Fij = force_jk*bond_jk->Rhat;
	      j->F += Fij; k->F -= Fij;
	      cnt++;
	    }
	  }
	}
      }
    }
  }
}
