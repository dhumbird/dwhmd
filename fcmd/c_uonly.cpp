#include "config.h"
extern map<short,double> A, B, MORSE_LAM, MORSE_MU,
SPL_RA, SPL_RB, MO_A, MO_EPS, MO_S, SPL_A, SPL_B, SPL_C;
//**********************************************************************
double config::Uonly(){
  double uu=0;
  double NF_ij, NC_ij, NSi_ij, Nt_ij, NF_ji, NC_ji, NSi_ji, Nt_ji;
  double P_ij, P_ji, dFP_ij, dFP_ji, dCP_ij, dCP_ji, dlam;
  double el, g, g1, xik, Yik, Y1ik, F_ij, dFi_ij, dFj_ij, dFc_ij;
  double eta_i, eta_j, delta_i, delta_j;

  for (VAI i_=begin; i_<end; i_++){
    //calculate VA_ij, VR_ij;
    for (VNI ij=i_->nlist.begin(); ij!=i_->nlist.end(); ij++){
      atom* j = (atom*)ij->a2;
      if (j > &(*i_)){
	scr_Rhat.clear();
	iNconj=jNconj=0;
	int type=i_->id+j->id;
	double r=ij->r;
	  	  
	VA_ij = ij->f * B[type] * SpExp(-MORSE_MU[type] * r);
		
	if (r < SPL_RA[type]){
	  double rija=r/MO_A[type];
	  VR_ij = MO_EPS[type]/r*(0.35*SpExp(-0.3*rija)
			   +0.55*SpExp(-1.2*rija)
			   +0.1*SpExp(-6*rija)) + MO_S[type];
	}
	else if (r < SPL_RB[type]){
	  VR_ij = SPL_C[type] + SpExp(SPL_A[type] * r + SPL_B[type]);
	}
	else{
	  VR_ij = ij->f*A[type] * SpExp(-MORSE_LAM[type] * r);
	}

	NF_ij  = i_->Nmap[9]  - (j->id==9)  * ij->f;
	NC_ij  = i_->Nmap[6]  - (j->id==6)  * ij->f;
	NSi_ij = i_->Nmap[14] - (j->id==14) * ij->f;
	NF_ji  = j->Nmap[9]  - (i_->id==9)  * ij->f;
	NC_ji  = j->Nmap[6]  - (i_->id==6)  * ij->f;
	NSi_ji = j->Nmap[14] - (i_->id==14) * ij->f;
	
	Nt_ij = NF_ij + NC_ij + NSi_ij;
	Nt_ji = NF_ji + NC_ji + NSi_ji;
	
	P_ij=P_ji=dFP_ij=dFP_ji=dCP_ij=dCP_ji=0;
  
	if (ij->type==12){
	  Pcc_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
	  Pcc_bicubicint(NF_ji, NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
	}
	else if (ij->type==15){
	  if (i_->id==6) Pcf_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
	  else Pcf_bicubicint(NF_ji, NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
	}
	bsp_ij=P_ij; bsp_ji=P_ji;
	
	//***************for k neighbor of i*************************
	for (VNI bond_ik=i_->nlist.begin(); bond_ik!=i_->nlist.end(); bond_ik++){
	  atom* k=(atom*)bond_ik->a2;
	  if (k!=j){
	    //*******b-sigma-pi calculations************
	    elambda(i_->id, j->id, k->id, ij->r, bond_ik->r, &el, &dlam);
	    G_angle(i_->id, j->id, k->id,
		    ij->Rhat ^ bond_ik->Rhat, &g, &g1);
	    bsp_ij += g*el*bond_ik->f;
	    //*********conjugation*****************
	    if (ij->type==12 && k->id==6){
	      PolySwitch(k->Nt - bond_ik->f - 2, &Yik, &Y1ik);
	      iNconj += bond_ik->f * Yik;
	    }
	  }
	}
	//***************for k neighbor of j*************************
	for (VNI bond_jk=j->nlist.begin(); bond_jk!=j->nlist.end(); bond_jk++){
	  atom* k=(atom*)bond_jk->a2;
	  if (k!=&(*i_)){
	    //*********b-sigma-pi calculation*************
	    elambda(j->id, i_->id, k->id, ij->r, bond_jk->r, &el, &dlam);
	    G_angle(j->id, i_->id, k->id, 
		    -(ij->Rhat ^ bond_jk->Rhat), &g, &g1);
	    bsp_ji += g * el*bond_jk->f;
	    //*********conjugation************
	    if (ij->type==12 && k->id==6){
	      PolySwitch(k->Nt - bond_jk->f - 2, &Yik, &Y1ik);
	      jNconj += bond_jk->f * Yik;
	    }
	  }
	}
	
	//***************conjugation************************
	if (ij->type==12){
	  Nconj_ij=1+iNconj+jNconj;
	  Fcc_tricubicint(Nt_ij, Nt_ji, Nconj_ij, 
			  &F_ij, &dFi_ij, &dFj_ij, &dFc_ij);
	}	  
	else 
	  F_ij=dFi_ij=dFj_ij=dFc_ij=0;
	
	set_eta(i_->id, j->id, &eta_i, &delta_i);
	set_eta(j->id, i_->id, &eta_j, &delta_j);
	bbar_ij = 0.5*(pow(1 + pow(bsp_ij, eta_i), -delta_i)
		       + pow(1 + pow(bsp_ji, eta_j), -delta_j)
		       +F_ij);
	uu+=(VR_ij-bbar_ij*VA_ij);
      }
    }
  }
  return uu;
}
