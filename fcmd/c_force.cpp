#include "config.h"
extern map<short,double> A, B, MORSE_LAM, MORSE_MU, RE,
SPL_RA, SPL_RB, MO_A, MO_EPS, MO_S, SPL_A, SPL_B, SPL_C;

//**********************************************************************
void config::ForceEval()
{
  for (VAI i_ = begin; i_ < end; i_++)
  {
    if (i_->ix == inert)
    {
      for (set<atom*>::iterator ij=ion_nbr.begin(); ij!=ion_nbr.end(); ij++)
      {
      	atom* j=*ij;
      	svector R=i_->R - j->R; R.minimg(Lx,Ly);
      	double r=R.mag();
      	int type = i_->id+j->id;
      	double moa=MO_A[type];
      	double moeps=MO_EPS[type];
      	double rija=r/moa;
      	Fij=moeps/r*(0.35*SpExp(-0.3*rija)*(1/r+0.3/moa)
      		     +0.55*SpExp(-1.2*rija)*(1/r+1.2/moa)
      		     +0.1*SpExp(-6*rija)*(1/r+6/moa))/r*R;
      	i_->F += Fij; j->F -= Fij;
      	u+= moeps/r*(0.35*SpExp(-0.3*rija)
      		     +0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
      }
    }
    else
    {
      for (VNI ij=i_->nlist.begin(); ij!=i_->nlist.end(); ij++)
      {
      	atom* j = (atom*)ij->a2;
      	if (j > &(*i_))
        {
    	    int type=i_->id+j->id;
    	    double r=ij->r;
          double aa = B[type] * SpExp(-MORSE_MU[type] * r);
          VA_ij = ij->f * aa;
          dVA_ij = aa * (ij->fprime - ij->f*MORSE_MU[type]);
      	  if (r < SPL_RA[type])
          {
      	    double moa=MO_A[type];
      	    double moeps=MO_EPS[type];
      	    double rija=r/moa;
      	    VR_ij = moeps/r*(0.35*SpExp(-0.3*rija)
      			     +0.55*SpExp(-1.2*rija)
      			     +0.1*SpExp(-6*rija)) + MO_S[type];
      	    dVR_ij = -moeps/r*(0.35*SpExp(-0.3*rija)*(1/r+0.3/moa)
      			      +0.55*SpExp(-1.2*rija)*(1/r+1.2/moa)
      			      +0.1*SpExp(-6*rija)*(1/r+6/moa));
      	  }
      	  else if (r < SPL_RB[type])
          {
      	    VR_ij = SPL_C[type] + SpExp(SPL_A[type] * r + SPL_B[type]);
      	    dVR_ij = SPL_A[type] * SpExp(SPL_A[type] * r + SPL_B[type]);
      	  }
      	  else
          {
      	    aa=A[type] * SpExp(-MORSE_LAM[type] * r);
      	    VR_ij = ij->f * aa;
      	    dVR_ij = aa * (ij->fprime - ij->f*MORSE_LAM[type]);
      	  }
          double b_ij=BondOrder(&(*ij), VA_ij);
          //cerr<<b_ij<<endl;
      	  u+=VR_ij - b_ij*VA_ij;
      	  Fij=(-dVR_ij+b_ij*dVA_ij)*ij->Rhat;
      	  i_->F+=Fij; j->F-=Fij;
      	}
      }
    }
  }
}
//******************************************************
double config::BondOrder(nbr* bond_ij, double Pre){
  svector Rik, Rjk, R;
  double rik, rjk;
  iFv_ij=0; iFv_ik.clear(); iFv_jk.clear(); iFv_ik2.clear();
  jFv_ij=0; jFv_ik.clear(); jFv_jk.clear(); jFv_jk2.clear();
  scr_Rhat.clear();
  iNconj=jNconj=0;
  
  atom* i=(atom*)bond_ij->a1;
  atom* j=(atom*)bond_ij->a2;
 
  double NF_ij, NC_ij, NSi_ij, Nt_ij, NF_ji, NC_ji, NSi_ji, Nt_ji;
  double NCl_ij, NCl_ji;
  double P_ij, P_ji, dFP_ij, dFP_ji, dCP_ij, dCP_ji, dlam;
  double cosO, el, elf, g, g1, xik, Yik, Y1ik, F_ij, dFi_ij, dFj_ij, dFc_ij;
  double eta_i, eta_j, delta_i, delta_j;

  NF_ij  = i->Nmap[9]  - (j->id==9)  * bond_ij->f;
  NC_ij  = i->Nmap[6]  - (j->id==6)  * bond_ij->f;
  NSi_ij = i->Nmap[14] - (j->id==14) * bond_ij->f;
  NCl_ij = i->Nmap[17] - (j->id==17) * bond_ij->f;
  NF_ji  = j->Nmap[9]  - (i->id==9)  * bond_ij->f;
  NC_ji  = j->Nmap[6]  - (i->id==6)  * bond_ij->f;
  NSi_ji = j->Nmap[14] - (i->id==14) * bond_ij->f;
  NCl_ji = j->Nmap[17] - (i->id==17) * bond_ij->f;

  Nt_ij = NF_ij + NC_ij + NSi_ij + NCl_ij;
  Nt_ji = NF_ji + NC_ji + NSi_ji + NCl_ji;

  P_ij=P_ji=dFP_ij=dFP_ji=dCP_ij=dCP_ji=0;
  
  if (bond_ij->type==12)
  {
    Pcc_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
    Pcc_bicubicint(NF_ji, NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
  }
  else if (bond_ij->type==15)
  {
    if (i->id==6) Pcf_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
    else Pcf_bicubicint(NF_ji, NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
  }
  else if (bond_ij->type==23)
  {
    if (i->id==14) Psif_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
    else Psif_bicubicint(NF_ji,  NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
  }
  else if (bond_ij->type==31)
  {
    if (i->id==14) Psicl_bicubicint(NCl_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
    else Psicl_bicubicint(NCl_ji,  NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
  }
  bsp_ij=P_ij; bsp_ji=P_ji;

  //***************for k neighbor of i*************************
  for (VNI bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++)
  {
    atom* k=(atom*)bond_ik->a2;
    if (k!=j)
    {
      VNI bond_jk=j->FindNbr(k);
      if (bond_jk!=j->nlist.end())
      {
      	rjk=bond_jk->r;
      	Rjk=bond_jk->Rhat;
      }
      else
      {
      	Rjk=j->R-k->R;
      	Rjk.minimg(Lx,Ly);
      	rjk=Rjk.mag();
      	Rjk/=rjk;
      }
      scr_Rhat.push_back(Rjk);
      //*******b-sigma-pi calculations************
      elambda(i->id, j->id, k->id, bond_ij->r, bond_ik->r, &el, &dlam);
      elf=el*bond_ik->f;
      cosO = bond_ij->Rhat ^ bond_ik->Rhat;
      G_angle(i->id, j->id, k->id, cosO, &g, &g1);
      //cerr<<cosO<<" "<<g<<endl;
      iFv_ij+=(elf*(g1*(1/bond_ik->r - cosO/bond_ij->r) + g*dlam));
      iFv_ik.push_back(bond_ik->fprime*
		       (g*el + ((k->id==9||k->id==17) ? dFP_ij : dCP_ij))
		       +elf*(g1*(1/bond_ij->r - cosO/bond_ik->r) - g*dlam));
      iFv_jk.push_back(elf*(-g1*rjk/bond_ik->r/bond_ij->r));
            
      bsp_ij += g*elf;
      //*********conjugation*****************
      if (bond_ij->type==12 && k->id==6)
      {
      	PolySwitch(k->Nt - bond_ik->f - 2, &Yik, &Y1ik);
      	iNconj += bond_ik->f * Yik;
      	iFv_ik2.push_back(bond_ik->fprime * Yik);
      	bond_ik->Fv_kl=bond_ik->f*Y1ik;
      }
      else
      {
      	iFv_ik2.push_back(0); 
      	bond_ik->Fv_kl=0;
      }
    }
  }
  //cerr<<bsp_ij<<" ";
  //***************for k neighbor of j*************************
  for (VNI bond_jk=j->nlist.begin(); bond_jk!=j->nlist.end(); bond_jk++)
  {
    atom* k=(atom*)bond_jk->a2;
    if (k!=i)
    {
      VNI bond_ik=i->FindNbr(k);
      if (bond_ik!=i->nlist.end())
      {
      	rik=bond_ik->r;
      	Rik=bond_ik->Rhat;
      }
      else
      {
      	Rik=i->R-k->R;
      	Rik.minimg(Lx,Ly);
      	rik=Rik.mag();
      	Rik/=rik;
      }
      scr_Rhat.push_back(Rik);
      //*********b-sigma-pi calculation*************
      elambda(j->id, i->id, k->id, bond_ij->r, bond_jk->r, &el, &dlam);
      elf=el*bond_jk->f;
      cosO = -(bond_ij->Rhat ^ bond_jk->Rhat);
      G_angle(j->id, i->id, k->id, cosO, &g, &g1);

      jFv_ij+=(elf*(g1*(1/bond_jk->r - cosO/bond_ij->r) + g*dlam));
      jFv_ik.push_back(elf*(-g1*rik/bond_jk->r/bond_ij->r));
      jFv_jk.push_back(bond_jk->fprime
		       *(g*el + ((k->id==9||k->id==17) ? dFP_ji : dCP_ji))
		       +elf*(g1*(1/bond_ij->r - cosO/bond_jk->r) - g*dlam));
      bsp_ji += g * elf;
      //*********conjugation************
      if ( bond_ij->type==12 && k->id==6)
      {
      	PolySwitch(k->Nt - bond_jk->f - 2, &Yik, &Y1ik);
      	jNconj += bond_jk->f * Yik;
      	jFv_jk2.push_back(bond_jk->fprime * Yik);
      	bond_jk->Fv_kl=bond_jk->f*Y1ik;
      }
      else
      {
      	jFv_jk2.push_back(0); 
      	bond_jk->Fv_kl=0;
      }
    }
  }
  //cerr<<bsp_ji<<endl;
  //***************conjugation************************
  if (bond_ij->type==12)
  {
    Nconj_ij=1+iNconj+jNconj;
    Fcc_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &F_ij, &dFi_ij, &dFj_ij, &dFc_ij);
  }
  else 
    F_ij=dFi_ij=dFj_ij=dFc_ij=0;
  //***********************************************
  double dtemp=bsp_ji;
  set_eta(i->id, j->id, &eta_i, &delta_i);
  set_eta(j->id, i->id, &eta_j, &delta_j);
  Pre*=0.5;
  bbar_ij = 0.5*(pow(1 + pow(bsp_ij, eta_i), -delta_i)
		 + pow(1 + pow(bsp_ji, eta_j), -delta_j)
		 +F_ij);

  //bsp becomes temp variable
  if (eta_i!=1)
  {
    if (bsp_ij <= 0) bsp_ij=0; //unphysical, but happens due to roundoff
    else
      bsp_ij=pow(1 + pow(bsp_ij, eta_i), -delta_i-1)*
	               eta_i*pow(bsp_ij, eta_i-1);
  }
  else
    bsp_ij=pow(1 + bsp_ij, -delta_i-1);
  bsp_ij*=-Pre*delta_i;
  
  if (eta_j!=1)
  {
    if (bsp_ji <= 0) bsp_ji=0;
    else
      bsp_ji=pow(1 + pow(bsp_ji, eta_j), -delta_j-1)*
	               eta_j*pow(bsp_ji, eta_j-1);
  }
  else  
    bsp_ji=pow(1 + bsp_ji, -delta_j-1);
  
  bsp_ji*=-Pre*delta_j;
  Fij=(bsp_ij*iFv_ij+bsp_ji*jFv_ij)*bond_ij->Rhat;
  i->F+=Fij; j->F-=Fij;
  cnt=0;
  scrcount=0;
  for (VNI bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++)
  {
    atom* k=(atom*)bond_ik->a2;
    if (k != j)
    {
      force_ik=bsp_ij*iFv_ik[cnt];
      if (bond_ik->fprime && dFi_ij) force_ik+=Pre*bond_ik->fprime*dFi_ij;
      if (dFc_ij && k->id==6)
      {
      	Sprime = Pre * dFc_ij;
      	if (iFv_ik2[cnt])
        {
          force_ik+=Sprime*iFv_ik2[cnt];
        }
    	if (bond_ik->Fv_kl)
      {
    	  Sprime*=bond_ik->Fv_kl;
    	  for (VNI bond_kl=k->nlist.begin(); bond_kl!=k->nlist.end(); bond_kl++)
    	    if (bond_kl->fprime)
          {
    	      atom* l=(atom*)bond_kl->a2;
    	      if (l != i)
            {
          		Fij=(Sprime*bond_kl->fprime)*bond_kl->Rhat;
          		k->F += Fij; l->F -= Fij;
            }
    	    }
        }
      }
      Fij = force_ik*bond_ik->Rhat;
      i->F += Fij; k->F -= Fij;
      
      Fij = bsp_ij * iFv_jk[cnt] * scr_Rhat[scrcount++];
      j->F += Fij; k->F -= Fij;
      cnt++;
    }
  }
  cnt=0;
  for (VNI bond_jk=j->nlist.begin(); bond_jk!=j->nlist.end(); bond_jk++)
  {
    atom* k=(atom*)bond_jk->a2;
    if (k != i)
    {
      force_jk=bsp_ji*jFv_jk[cnt];
      if (bond_jk->fprime && dFj_ij) 
        force_jk+=Pre*bond_jk->fprime*dFj_ij;
      if (dFc_ij && k->id==6)
      {
      	Sprime = Pre * dFc_ij;
      	if (jFv_jk2[cnt]) 
          force_jk+=Sprime*jFv_jk2[cnt];
      	if (bond_jk->Fv_kl)
        {
      	  Sprime*=bond_jk->Fv_kl;
          for (VNI bond_kl=k->nlist.begin(); bond_kl!=k->nlist.end(); bond_kl++)
            if (bond_kl->fprime)
            {
      	      atom* l=(atom*)bond_kl->a2;
      	      if (l!=j)
              {
            		Fij = (Sprime*bond_kl->fprime)*bond_kl->Rhat;
            		k->F += Fij; l->F -= Fij;
              }
            }
        }
      }
      Fij = bsp_ji*jFv_ik[cnt]*scr_Rhat[scrcount++];
      i->F += Fij; k->F -= Fij;
      Fij = force_jk*bond_jk->Rhat;
      j->F += Fij; k->F -= Fij;
      cnt++;
    }
  }
  return bbar_ij;
}
//*******************************************************************
void config::elambda(short i_id, short j_id, short k_id, 
		     double r_ij, double r_ik,
		     double* el, double* dlam){
  double alpha=0; double beta=1;
  if (i_id==9){
    alpha=TAN_alpha; beta=TAN_beta;
  }
  if (i_id==6 && j_id==9 && k_id==9){
    alpha=TAN_alpha; beta=TAN_beta;
  }
  if (i_id==14){
    if (j_id==14 && k_id==14){
      alpha=TER_alpha; beta=TER_beta;
    }
    if (j_id==9 || k_id==9){
      alpha=MUR_alpha; beta=MUR_beta;
    }
    if (j_id==17 || k_id==17){
      alpha=MUR_alpha; beta=MUR_beta;
    }
  }
  if (i_id==17){
    alpha=TAN_alpha; 
    beta=TAN_beta;
  }
  if (alpha){
    double n = (r_ij - RE[i_id + j_id]) - (r_ik - RE[i_id + k_id]);
    *el=SpExp(alpha*pow(n, beta));
    *dlam=alpha;
    if (beta!=1){
      *dlam*=beta*pow(n, beta-1);
    }
  }
  else{
    *el = 1.0; *dlam=0;
  }
}
//******************************************************************
void config::G_angle(short i_id, short j_id, short k_id, 
		     double cosO, double* g, double* g1)
{
  //cerr<<i_id<<" "<<j_id<<" "<<k_id<<endl;
  //cerr<<cosO<<" ";
  if (i_id==6)
  {
    if (j_id==14)
    { //Tersoff
      double c0 = -0.57058 - cosO;
      double c02 = c0*c0;
      *g = 1.5724e-07 * (1 + 7.5326e-9 - 1.4477264e-07 / (19.219456 + c02));
      *g1 = -4.55281e-14 * c0 / pow((19.219456 + c02),2);
    }
    else
    { //Brenner
      double c0 = -1 - cosO;
      double c02 = c0*c0;
      *g = 0.00020813 * (8890.8 - 108900 / (12.25 + c02));
      *g1 = -45.330714 * c0 / pow((12.25 + c02),2);
    }
  }
  else if (i_id==14)
  {
    if (j_id==9 || k_id==9)
    { //Murty
      *g = 0.0216 + 0.27*pow(-0.47-cosO, 2);
      *g1 = -0.54*(-0.47-cosO);
    }
    //corrected Jan 2018. Previously "if"
    else if (j_id==17 || k_id==17)
    //if (j_id==17 || k_id==17)
    { //Murty
      *g = 0.0216 + 0.27*pow(-0.47-cosO, 2);
      *g1 = -0.54*(-0.47-cosO);
    }
    else
    { //Tersoff
      *g = 0.16*pow(-0.59826-cosO, 2);
      *g1 = -0.32*(-0.59826-cosO);
    }
  }
  else if (i_id==9)
  { //Tanaka
    *g=4;
    *g1=0;
  }
  else if (i_id==17)
  { //Tanaka
    *g=4;
    *g1=0;
  }

}
//***********************************************************************
void config::set_eta(short i_id, short j_id, double* eta, double* delta){
  if (i_id==6){
    if (j_id==14){ //C-Tersoff
      *eta=TERC_eta; *delta=TERC_delta;
    }
    else{ //Brenner
      *eta=BRE_eta; *delta=BRE_delta;
    }
  }
  else if (i_id==14){
    if (j_id==9){
      *eta=MUR_eta; *delta=MUR_delta;
    }
    else if (j_id==17){
      *eta=MUR_eta; *delta=MUR_delta;
    }
    else{
      *eta=TERSi_eta; *delta=TERSi_delta;
    }
  }
  else if (i_id==9){
    *eta=TAN_eta; *delta=TAN_delta;
  }
  else if (i_id==17){
    *eta=TAN_eta; *delta=TAN_delta;
  }
}
