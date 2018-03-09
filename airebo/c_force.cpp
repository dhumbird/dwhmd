#include "config.h"
extern map<short,double> Q, A, ALPHA, B1, B2, B3, BETA1, BETA2, BETA3,
LJ_EPS, LJ_SIGMA, B_MIN, B_MAX, REBO_F_MIN, V_LJ_RC;

//**********************************************************************
void config::ForceEval(){
  for (VAI i_ = begin; i_ < end; i_++){
    if (i_->ix == inert){
      for (VNI bond_ij=i_->nlist.begin(); bond_ij!=i_->nlist.end(); bond_ij++){
	atom* j=(atom*)bond_ij->a2;
	double r=bond_ij->r;
	double moa, moeps;
	if (j->id == 6){
	  moa=MO_AR_C_a;
	  moeps=MO_AR_C_EPS;
	}
	else if (j->id == 1){
	  moa=MO_AR_H_a;
	  moeps=MO_AR_H_EPS;
	}
	double rija=r/moa;
	Fij=moeps/r*(0.35*SpExp(-0.3*rija)*(1/r+0.3/moa)
		     +0.55*SpExp(-1.2*rija)*(1/r+1.2/moa)
		     +0.1*SpExp(-6*rija)*(1/r+6/moa))*bond_ij->Rhat;
	i_->F += Fij; j->F -= Fij;
	u+= moeps/r*(0.35*SpExp(-0.3*rija)
		     +0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
      }
    }
    else{
      double V_ij, dV_ij, Sb;
      //LENNARD-JONES
      for (VNJI ij=i_->ljnlist.begin(); ij!=i_->ljnlist.end(); ij++){
	atom* j = (atom*)ij->a2;
	ij->PreComp();
	V_ij= LJ_EPS[ij->type]*
	  (pow(LJ_SIGMA[ij->type]/ij->r_true,12)
	   -pow(LJ_SIGMA[ij->type]/ij->r_true,6)) - V_LJ_RC[ij->type];
	dV_ij = LJ_EPS[ij->type]*
	  (-12*pow(LJ_SIGMA[ij->type]/ij->r_true,12)/ij->r_true
	   +6*pow(LJ_SIGMA[ij->type]/ij->r_true,6)/ij->r_true);
	   
	if (ij->f_LJ){
	  ij->r = REBO_F_MIN[ij->type];
	  ij->f = ij->w = 1; 
	  ij->fprime = ij->wprime = 0;
	  Sb=BondOrder(&(*ij), -ij->C*V_ij*ij->f_LJ, 1);
	}
	else 
	  Sb=0;

	u += ij->C * V_ij * (1 - ij->f_LJ + ij->f_LJ * Sb);
	Fij = - ij->C * ((dV_ij*(1 + ij->f_LJ * (Sb-1)))
			 + V_ij * (Sb-1) * ij->fprime_LJ)* ij->Rhat;
 	i_->F+=Fij; j->F-=Fij;
	
	V_ij*=(1 + ij->f_LJ * (Sb-1));
	if (ij->n3){
	  Fij=V_ij* 
	    ij->n1->fprime * ij->n2->f * ij->n3->f * ij->n1->Rhat;
	  ((atom*)ij->n1->a1)->F+=Fij; ((atom*)ij->n1->a2)->F-=Fij;
	  Fij=V_ij*
	    ij->n1->f * ij->n2->fprime * ij->n3->f * ij->n2->Rhat;
	  ((atom*)ij->n2->a1)->F+=Fij; ((atom*)ij->n2->a2)->F-=Fij;
	  Fij=V_ij*
	    ij->n1->f * ij->n2->f * ij->n3->fprime * ij->n3->Rhat;
	  ((atom*)ij->n3->a1)->F+=Fij; ((atom*)ij->n3->a2)->F-=Fij;
	}
	else if (ij->n2){
	  Fij=V_ij* 
	    ij->n1->fprime * ij->n2->f * ij->n1->Rhat;
	  ((atom*)ij->n1->a1)->F+=Fij; ((atom*)ij->n1->a2)->F-=Fij;
	  Fij=V_ij* 
	    ij->n1->f * ij->n2->fprime * ij->n2->Rhat;
	  ((atom*)ij->n2->a1)->F+=Fij; ((atom*)ij->n2->a2)->F-=Fij;
	}
	else if (ij->n1){
	  Fij=V_ij* 
	    ij->n1->fprime * ij->n1->Rhat;
	  ((atom*)ij->n1->a1)->F+=Fij; ((atom*)ij->n1->a2)->F-=Fij;
	}
      }
    
    //NON - LENNARD-JONES
      for (VNI ij=i_->nlist.begin(); ij!=i_->nlist.end(); ij++){
	atom* j = (atom*)ij->a2;
	if (j > &(*i_)){
	  int type=ij->type;
	  double r=ij->r;
	  b1=B1[type]; b2=B2[type]; b3=B3[type];
	  beta1=BETA1[type]; beta2=BETA2[type]; beta3=BETA3[type];
	  VA_ij=b1*exp(-beta1*r) + b2*exp(-beta2*r) + b3*exp(-beta3*r);
	  dVA_ij=b1*beta1*exp(-beta1*r) + b2*beta2*exp(-beta2*r)
	    +b3*beta3*exp(-beta3*r);
	  dVA_ij = VA_ij * ij->fprime - ij->f * dVA_ij;
	  VA_ij *= ij->f;
	  double qq=Q[type]/r; 
	  double aa=A[type] * exp(-ALPHA[type] * r);
	  VR_ij = ij->f * (1+qq) * aa;
	  dVR_ij= aa * (ij->fprime*(1+qq)-ij->f*(qq/r+ALPHA[type]*(1+qq)));
	  double b_ij=BondOrder(&(*ij), VA_ij, 0);
	  u+=VR_ij - b_ij*VA_ij;
	  Fij=(-dVR_ij+b_ij*dVA_ij)*ij->Rhat;
	  i_->F+=Fij; j->F-=Fij;
	}
      }
    }
  }
}
//******************************************************
double config::BondOrder(nbr* bond_ij, double Pre, bool spoof){
  svector Ril, Rjk, Rij, Rjl, Rik, R;
  double ril, rjk;
  double sumdg_i,sumdg_j,iNconj,jNconj,bsp_ij,bsp_ji;
  vector<double> i_cosO, j_cosO, dctjk, dctij, dctik, dctji, dctil, dctjl;
  vector<double> iFv_ik, iFv_jk, iFv_ik2, iFv_kl;
  vector<double> jFv_il, jFv_jl, jFv_jl2, jFv_lm;
  vector<svector> Ril_save, Rjk_save;
  double iFv_ij, jFv_ij, force_ik, force_jk, force_il, force_jl;
  double cosO, el, elf, g, g1, xik, Yik, Y1ik, F_ij, dFi_ij, dFj_ij, dFc_ij;
  double G, Gprime, gam, gamprime, bbar_ij;
  double Nconj_ij, bDH_ij;
  int icnt, jcnt;
  double T_ij, dTi_ij, dTj_ij, dTc_ij;
  double e_tors, v_tors, dv_tors;
  map<atom*, svector> tmpF; 
  svector C_jik, C_ijl;

  iFv_ij=jFv_ij=0; 
  sumdg_i=sumdg_j=iNconj=jNconj=bsp_ij=bsp_ji=0;

  atom* i=(atom*)bond_ij->a1;
  atom* j=(atom*)bond_ij->a2;
  double rij_true; svector Rij_true; nbr* real_ij=NULL;

  if (spoof){
    if (i->FindNbr(j)!=i->nlist.end()){
      real_ij = &(*(i->FindNbr(j)));
    }
    Rij_true = (i->R - j->R); Rij_true.minimg(Lx,Ly);
    rij_true = Rij_true.mag();
  }
  
  double NH_ij, NC_ij, Nt_ij, NH_ji, NC_ji, Nt_ji;
  if (!spoof){
    NH_ij = i->NH - (j->id==1)*bond_ij->f;
    NC_ij = i->NC - (j->id==6)*bond_ij->f;
    NH_ji = j->NH - (i->id==1)*bond_ij->f;
    NC_ji = j->NC - (i->id==6)*bond_ij->f; 
  }
  else{
    if (real_ij){
      //The REAL bond_ij->f was put in by precomputations and must be removed.
      //bond_ij->f at this point is just 1
      //Nt_ij(spoof) must equal Nt_ij(real)
      NH_ij = i->NH - (j->id==1)*(real_ij->f);
      NC_ij = i->NC - (j->id==6)*(real_ij->f);
      NH_ji = j->NH - (i->id==1)*(real_ij->f);
      NC_ji = j->NC - (i->id==6)*(real_ij->f); 
    }
    else{
      NH_ij = i->NH;
      NC_ij = i->NC;
      NH_ji = j->NH;
      NC_ji = j->NC;
    }
  }
  Nt_ij = NH_ij + NC_ij;
  Nt_ji = NH_ji + NC_ji;

  PolySwitch(2*(Nt_ij-3.2), &S, &Sprime);
  Sprime*=2;
  //***************for k neighbor of i*************************
  for (VNI bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
    atom* k=(atom*)bond_ik->a2;
    if (k!=j){
      Rjk=bond_ik->r*bond_ik->Rhat - bond_ij->r*bond_ij->Rhat;
      Rjk.minimg(Lx,Ly);
      rjk=Rjk.mag();
      Rjk/=rjk;
      Rjk_save.push_back(Rjk);
      //*******b-sigma-pi calculations************
      if (i->id == 1){
	el=exp(4.0*((k->id==1)*RHH + (k->id==6)*RCH - bond_ik->r
		    -(j->id==1)*RHH - (j->id==6)*RCH + bond_ij->r));
	dlam1=4.0; dlam2=-4.0; 
	elf=bond_ik->f*el;
      }
      else{
	el = 1.0; dlam1=dlam2=0;
	elf=bond_ik->f;
      }
	
      cosO = bond_ij->Rhat ^ bond_ik->Rhat;
      Gpoly(cosO, i->id, &G, &Gprime, &gam, &gamprime);
      i_cosO.push_back(cosO);

      dctjk.push_back(-rjk/bond_ij->r/bond_ik->r);
      dctik.push_back(1/bond_ij->r - cosO/bond_ik->r);
      if (spoof) dctij.push_back(0);
      else dctij.push_back(1/bond_ik->r - cosO/bond_ij->r);

      if (gam!=0){
	g = G+S*(gam-G);
	g1=Gprime+S*(gamprime-Gprime);
	sumdg_i+=elf*(gam-G);
      }
      else{
	g=G; g1=Gprime;
      }

      if (!spoof) iFv_ij+=(elf*(g1*dctij.back() + g*dlam1));
      iFv_ik.push_back(bond_ik->fprime*g*el+elf*(g1*dctik.back() + g*dlam2));
      iFv_jk.push_back(elf*g1*dctjk.back());
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
  for (VNI bond_jl=j->nlist.begin(); bond_jl!=j->nlist.end(); bond_jl++){
    atom* l=(atom*)bond_jl->a2;
    if (l!=i){
      Ril=bond_ij->r*bond_ij->Rhat + bond_jl->r*bond_jl->Rhat;
      Ril.minimg(Lx,Ly);
      ril=Ril.mag();
      Ril/=ril;
      Ril_save.push_back(Ril);
      //*********b-sigma-pi calculation*************
      if (j->id == 1){
	el=exp(4.0*((l->id==1)*RHH + (l->id==6)*RCH - bond_jl->r
		    -(i->id==1)*RHH - (i->id==6)*RCH + bond_ij->r));
	dlam1=4.0; dlam2=-4.0; 
	elf=bond_jl->f*el;
      }
      else{
	el = 1.0; dlam1=dlam2=0;
	elf=bond_jl->f;
      }
      
      cosO = (-1*bond_ij->Rhat) ^ bond_jl->Rhat;
      Gpoly(cosO, j->id, &G, &Gprime, &gam, &gamprime);
      j_cosO.push_back(cosO);
      
      dctjl.push_back(1/bond_ij->r - cosO/bond_jl->r);
      dctil.push_back(-ril/bond_jl->r/bond_ij->r);
      if (spoof) dctji.push_back(0);
      else dctji.push_back(1/bond_jl->r - cosO/bond_ij->r);

      if (gam!=0){
	g = G+S*(gam-G);
	g1=Gprime+S*(gamprime-Gprime);
	sumdg_j+=elf*(gam-G);
      }
      else{
	g=G; g1=Gprime;
      }
      if (!spoof) jFv_ij+=elf*(g1*dctji.back() + g*dlam1);
      jFv_il.push_back(elf*g1*dctil.back());
      jFv_jl.push_back(bond_jl->fprime*g*el +elf*(g1*dctjl.back() + g*dlam2));
      bsp_ji += g * elf;
      //*********conjugation************
      if (l->id==6){
	PolySwitch(l->Nt - bond_jl->f - 2, &Yik, &Y1ik);
	jNconj += bond_jl->f * Yik;
	jFv_jl2.push_back(bond_jl->fprime * Yik);
	jFv_lm.push_back(bond_jl->f*Y1ik);
      }
      else{
	jFv_jl2.push_back(0); jFv_lm.push_back(0);
      }
    }
  }
  sumdg_j*=Sprime;
  
  
  double P_ij,P_ji,dHP_ij,dHP_ji,dCP_ij,dCP_ji;
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
  bsp_ij+=P_ij; bsp_ji+=P_ji;

  bsp_ij = 1.0 / sqrt(1 + bsp_ij);
  bsp_ji = 1.0 / sqrt(1 + bsp_ji);
	  
  //***************conjugation************************
  Nconj_ij=1+iNconj*iNconj+jNconj*jNconj;  
  if (bond_ij->type==0)
    Fcc_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &F_ij, &dFi_ij, &dFj_ij, &dFc_ij);
  else if (bond_ij->type==2)
    Fch_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &F_ij, &dFi_ij, &dFj_ij, &dFc_ij); 
  else if (bond_ij->type==1)
    Fhh_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &F_ij, &dFi_ij, &dFj_ij, &dFc_ij); 
  else 
    F_ij=dFi_ij=dFj_ij=dFc_ij=0;
  
  //******************  dihedral*****************************
  bDH_ij=T_ij=dTi_ij=dTj_ij=dTc_ij=0;
  if (bond_ij->type==0){
    Tcc_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &T_ij, &dTi_ij, &dTj_ij, &dTc_ij);
    if (T_ij || !spoof){
      if (spoof){
	tmpF[i]=svector(0,0,0);
	tmpF[j]=svector(0,0,0);
	for (VNI bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
	  atom* k=(atom*)bond_ik->a2;
	  tmpF[k]=svector(0,0,0);
	  for (VNI bond_jl=j->nlist.begin(); bond_jl!=j->nlist.end();
	       bond_jl++){
	    atom* l=(atom*)bond_jl->a2;
	    tmpF[l]=svector(0,0,0);
	  }
	}
      }
      Rij = bond_ij->r * bond_ij->Rhat;
      icnt = -1;
      for (VNI bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
	atom* k=(atom*)bond_ik->a2;
	if (k != j){
	  icnt++;
	  Rik=bond_ik->r*bond_ik->Rhat;
	  Rjk=Rjk_save[icnt];
	  svector Crk = Rik * Rij;
	  double cosjik=i_cosO[icnt];
	  double sin2jik=1-cosjik*cosjik;
	  jcnt=-1;
	  for (VNI bond_jl=j->nlist.begin(); bond_jl!=j->nlist.end();
	       bond_jl++){
	    atom* l=(atom*)bond_jl->a2;
	    if (l != i){
	      jcnt++;
	      Rjl=bond_jl->r*bond_jl->Rhat;
	      Ril=Ril_save[jcnt];
	      svector Crl = Rij * Rjl;
	      double cosijl=j_cosO[jcnt];
	      double sin2ijl=1-cosijl*cosijl;
	      double cwnom=bond_ik->r*bond_jl->r
		*bond_ij->r*bond_ij->r*sqrt(sin2jik)*sqrt(sin2ijl);

	      double dt1dik=1/bond_ik->r - dctik[icnt]/sin2jik*cosjik;
	      double dt1djk=-dctjk[icnt]/sin2jik*cosjik;
	      double dt1djl=1/bond_jl->r - dctjl[jcnt]/sin2ijl*cosijl;
	      double dt1dil=-dctil[jcnt]/sin2ijl*cosijl;
	      double dt1dij=2/bond_ij->r-dctij[icnt]/sin2jik*cosjik
		-dctji[jcnt]/sin2ijl*cosijl;
	      
	      double cwnum=Crk ^ Crl;
	      cosO=cwnum/cwnom;
	      double bt=1-cosO*cosO;
	      
	      svector dt2dik(-Rij.z*Crl.y+Rij.y*Crl.z,
			     -Rij.x*Crl.z+Rij.z*Crl.x,
			     -Rij.y*Crl.x+Rij.x*Crl.y);
	      svector dt2djl(-Rij.y*Crk.z+Rij.z*Crk.y,
			     -Rij.z*Crk.x+Rij.x*Crk.z,
			     -Rij.x*Crk.y+Rij.y*Crk.x);
	      svector dt2dij(Rik.z*Crl.y-Rjl.z*Crk.y-Rik.y*Crl.z+Rjl.y*Crk.z,
			     Rik.x*Crl.z-Rjl.x*Crk.z-Rik.z*Crl.x+Rjl.z*Crk.x,
			     Rik.y*Crl.x-Rjl.y*Crk.x-Rik.x*Crl.y+Rjl.x*Crk.y);
	      if (T_ij){
		bDH_ij+=bond_ik->w * bond_jl->w*bt;
		double aa, aaa1, at2;
		if (!spoof){
		  aa=-Pre*cosO/cwnom*T_ij*bond_jl->w*bond_ik->w;
		  aaa1=Pre*0.5*bt*T_ij;
		  at2=aa*cwnum;
		  //i-j
		  Fij=-dt1dij*at2*bond_ij->Rhat + aa*dt2dij;
		  i->F+=Fij; j->F-=Fij;
		  //i-k
		  Fij=(-dt1dik*at2+aaa1*bond_jl->w*bond_ik->wprime)
		    *bond_ik->Rhat + aa*dt2dik;
		  i->F+=Fij; k->F-=Fij;
		  //j-l
		  Fij=(-dt1djl*at2+aaa1*bond_jl->wprime*bond_ik->w)
		    *bond_jl->Rhat + aa*dt2djl;
		  j->F+=Fij; l->F-=Fij;
		  //j-k
		  Fij=-dt1djk*at2*Rjk;
		  j->F+=Fij; k->F-=Fij;
		  //i-l
		  Fij=-dt1dil*at2*Ril;
		  i->F+=Fij; l->F-=Fij;
 		}
		else{
		  aa=-cosO/cwnom*T_ij*bond_jl->w*bond_ik->w;
		  aaa1=0.5*bt*T_ij;
		  at2=aa*cwnum;
		  //i-j is zero
		  //i-k
		  Fij=(-dt1dik*at2+aaa1*bond_jl->w*bond_ik->wprime)
		    *bond_ik->Rhat + aa*dt2dik;
		  tmpF[i]+=Fij; tmpF[k]-=Fij;
		  //j-l
		  Fij=(-dt1djl*at2+aaa1*bond_jl->wprime*bond_ik->w)
		    *bond_jl->Rhat + aa*dt2djl;
		  tmpF[j]+=Fij; tmpF[l]-=Fij;
		  //j-k
		  double ff=-dt1djk*at2;
		  Fij=ff*Rjk;
		  tmpF[i]+=Fij; tmpF[k]-=Fij;
		  Fij *= (bond_ij->r/rij_true);
		  tmpF[i] -= Fij; tmpF[j] += Fij;
		  ff *= -bond_ij->r/pow(rij_true,3)*(Rij_true^Rjk);
		  Fij=ff*Rij_true;
		  tmpF[i]-=Fij; tmpF[j]+=Fij;
		  //i-l
		  ff=-dt1dil*at2;
		  Fij=ff*Ril;
		  tmpF[j]+=Fij; tmpF[l]-=Fij;
		  Fij*=(bond_ij->r/rij_true);
		  tmpF[i] += Fij; tmpF[j] -= Fij;
		  ff *= -bond_ij->r/pow(rij_true,3)*(Rij_true^Ril);
		  Fij=ff*Rij_true;
		  tmpF[i] += Fij; tmpF[j] -= Fij;
		}
	      }
	      //***AIREBO torsional potential
	      if (!spoof){
		switch (k->id + l->id){
		case 12: e_tors=eCCCC; break;
		case 7:  e_tors=eHCCC; break;
		case 2:  e_tors=eHCCH; break;
		}
		v_tors = tors_fac * e_tors
		  * pow(cosO/2.0,10) - 0.1*e_tors;
		dv_tors = tors_fac * e_tors * 5 * pow(cosO/2.0,9)
		  *(bond_ik->f*bond_ij->f*bond_jl->f);
		u+= bond_ik->f * bond_ij->f * bond_jl->f * v_tors;
		
		double aa=-1/cwnom*dv_tors;
		double at2=aa*cwnum;

		//i-j
		Fij=(-dt1dij*at2-v_tors*bond_jl->f*bond_ik->f*bond_ij->fprime)
		  *bond_ij->Rhat + aa*dt2dij;
		i->F+=Fij; j->F-=Fij;
		//i-k
		Fij=(-dt1dik*at2-v_tors*bond_jl->f*bond_ij->f*bond_ik->fprime)
		  *bond_ik->Rhat + aa*dt2dik;
		i->F+=Fij; k->F-=Fij;
		//j-l
		Fij=(-dt1djl*at2-v_tors*bond_ij->f*bond_ik->f*bond_jl->fprime)
		  *bond_jl->Rhat + aa*dt2djl;
		j->F+=Fij; l->F-=Fij;	      
		//j-k
		Fij=-dt1djk*at2*Rjk;
		j->F+=Fij; k->F-=Fij;
		//i-l
		Fij=-dt1dil*at2*Ril;
		i->F+=Fij; l->F-=Fij;
	      }
	    }
	  }
	}
      }
      if (T_ij || dTi_ij || dTj_ij || dTc_ij){
	dFi_ij += dTi_ij * bDH_ij;
	dFj_ij += dTj_ij * bDH_ij;
	dFc_ij += dTc_ij * bDH_ij;
	bDH_ij *= T_ij;
      }
    }
  }
  bbar_ij=(bsp_ij + bsp_ji + F_ij + bDH_ij)/2.0;
  bsp_ij=-0.25*pow(bsp_ij,3)*Pre; //bsp becomes temp variable
  bsp_ji=-0.25*pow(bsp_ji,3)*Pre;

  if (spoof){ //if b_ij is too big, we can jump out here.
    PolySwitch((bbar_ij-B_MIN[bond_ij->type])/
	       (B_MAX[bond_ij->type]-B_MIN[bond_ij->type]), &S, &Sprime);
    if (Sprime==0){
      return S;
    }
    else{
      Pre *= Sprime / (B_MAX[bond_ij->type]-B_MIN[bond_ij->type]);
      if (T_ij){
	i->F+=tmpF[i] * Pre;
	j->F+=tmpF[j] * Pre;
	for (VNI bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
	  atom* k=(atom*)bond_ik->a2;
	  if (k != j)
	    k->F+=tmpF[k]*Pre;
	}
	for (VNI bond_jl=j->nlist.begin(); bond_jl!=j->nlist.end(); bond_jl++){
	  atom* l=(atom*)bond_jl->a2;
	  if (l != i)
	    l->F+=tmpF[l]*Pre;
	}
      }
    }
  }
  else{
    Fij=(bsp_ij*iFv_ij+bsp_ji*jFv_ij)*bond_ij->Rhat;
    i->F+=Fij; j->F-=Fij;
  }
  
  icnt=0;
  for (VNI bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
    atom* k=(atom*)bond_ik->a2;
    if (k != j){
      Rjk=Rjk_save[icnt];
      force_ik=bsp_ij*iFv_ik[icnt];
      force_jk=bsp_ij*iFv_jk[icnt];
      if (bond_ik->fprime){
	force_ik+=bsp_ij*bond_ik->fprime
	  *(sumdg_i + (k->id==6)*dCP_ij + (k->id==1)*dHP_ij);
	if (dFi_ij) force_ik+=0.5*Pre*bond_ik->fprime*dFi_ij;
      }
      if (k->id==6 && dFc_ij){
	Sprime = Pre * dFc_ij * iNconj ;
	if (iFv_ik2[icnt]) force_ik+=Sprime*iFv_ik2[icnt];
	if (iFv_kl[icnt]){
	  Sprime*=iFv_kl[icnt];
	  for (VNI bond_kl=k->nlist.begin(); bond_kl!=k->nlist.end(); 
	       bond_kl++)
	    if (bond_kl->fprime){
	      atom* l=(atom*)bond_kl->a2;
	      if (l != i){
		Fij=(Sprime*bond_kl->fprime)*bond_kl->Rhat;
		k->F += Fij; l->F -= Fij;
	      }
	    }
	}
      }
      Fij = force_ik*bond_ik->Rhat;
      i->F += Fij; k->F -= Fij;
      Fij = force_jk * Rjk;
      if (spoof){
	//cerr<<force_jk<<endl;
	i->F += Fij; k->F -= Fij;
	Fij *= (bond_ij->r/rij_true);
	i->F -= Fij; j->F += Fij;
	force_jk *= -bond_ij->r/pow(rij_true,3)*(Rij_true^Rjk);
	Fij=force_jk*Rij_true;
	i->F-=Fij; j->F+=Fij;
      }
      else{
	j->F += Fij; k->F -= Fij;
      }
      icnt++;
    }
  }
  jcnt=0;
  for (VNI bond_jl=j->nlist.begin(); bond_jl!=j->nlist.end(); bond_jl++){
    atom* l=(atom*)bond_jl->a2;
    if (l !=i){
      Ril=Ril_save[jcnt];
      force_jl=bsp_ji*jFv_jl[jcnt];
      force_il=bsp_ji*jFv_il[jcnt];
      if (bond_jl->fprime){
	force_jl+=bsp_ji*bond_jl->fprime
	  *(sumdg_j + (l->id==6)*dCP_ji + (l->id==1)*dHP_ji);
	if (dFj_ij) force_jl+=0.5*Pre*bond_jl->fprime*dFj_ij;
      }
      if (l->id==6 && dFc_ij){
	Sprime=Pre*dFc_ij*jNconj;
	if (jFv_jl2[jcnt]) force_jl+=Sprime*jFv_jl2[jcnt];
	if (jFv_lm[jcnt]){
	  Sprime*=jFv_lm[jcnt];
	  for (VNI bond_lm=l->nlist.begin(); bond_lm!=l->nlist.end();
	       bond_lm++)
	    if (bond_lm->fprime){
	      atom* m=(atom*)bond_lm->a2;
	      if (m!=j){
		Fij = (Sprime*bond_lm->fprime)*bond_lm->Rhat;
		l->F += Fij; m->F -= Fij;
	      }
	    }
	}
      }
      Fij = force_il*Ril;
      if (spoof){
	//cerr<<force_il<<endl;
	j->F += Fij; l->F -= Fij;
	Fij *= (bond_ij->r/rij_true);
	j->F -= Fij; i->F += Fij;
	force_il *= -bond_ij->r/pow(rij_true,3)*(Rij_true^Ril);
	Fij=force_il*Rij_true;
	j->F-=Fij; i->F+=Fij;
      }
      else{
	i->F += Fij; l->F -= Fij;
      }
      Fij = force_jl*bond_jl->Rhat;
      j->F += Fij; l->F -= Fij;
      jcnt++;
    }
  }
  if (spoof) return S;
  else return bbar_ij;
}
