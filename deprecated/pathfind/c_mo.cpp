#include "config.h"

const double MO_AR_SI_a = 0.388874*pow(sqrt(14)+sqrt(18),-TWO_THIRDS);
const double MO_AR_F_a = 0.388874*pow(sqrt(9)+sqrt(18),-TWO_THIRDS);
const double MO_AR_CL_a = 0.388874*pow(sqrt(17)+sqrt(18),-TWO_THIRDS);

const double MO_NE_SI_a = 0.388874*pow(sqrt(14)+sqrt(10),-TWO_THIRDS);
const double MO_NE_F_a = 0.388874*pow(sqrt(9)+sqrt(10),-TWO_THIRDS);
const double MO_NE_CL_a = 0.388874*pow(sqrt(17)+sqrt(10),-TWO_THIRDS);

const double MO_HE_SI_a = 0.388874*pow(sqrt(14)+sqrt(2),-TWO_THIRDS);
const double MO_HE_F_a = 0.388874*pow(sqrt(9)+sqrt(2),-TWO_THIRDS);
const double MO_HE_CL_a = 0.388874*pow(sqrt(17)+sqrt(2),-TWO_THIRDS);

const double MO_KR_SI_a = 0.388874*pow(sqrt(14)+sqrt(36),-TWO_THIRDS);
const double MO_KR_F_a = 0.388874*pow(sqrt(9)+sqrt(36),-TWO_THIRDS);
const double MO_KR_CL_a = 0.388874*pow(sqrt(17)+sqrt(36),-TWO_THIRDS);

const double MO_AR_SI_EPS = 14.39965*14*18;
const double MO_AR_F_EPS = 14.39965*18*9;
const double MO_AR_CL_EPS = 14.39965*18*17;

const double MO_HE_SI_EPS = 14.39965*14*2;
const double MO_HE_F_EPS = 14.39965*2*9;
const double MO_HE_CL_EPS = 14.39965*2*17;

const double MO_NE_SI_EPS = 14.39965*10*14;
const double MO_NE_F_EPS = 14.39965*10*9;
const double MO_NE_CL_EPS = 14.39965*10*17;

const double MO_KR_SI_EPS = 14.39965*36*14;
const double MO_KR_F_EPS = 14.39965*36*9;
const double MO_KR_CL_EPS = 14.39965*36*17;

//*************************************************************************
double config::HeSi(particle* I, particle *J){
  rija=rij/MO_HE_SI_a;
  Fij=MO_HE_SI_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_HE_SI_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_HE_SI_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_HE_SI_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_HE_SI_EPS/rij*(0.35*SpExp(-0.3*rija)
			   +0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}
//*************************************************************************
double config::ArSi(particle* I, particle *J){
  rija=rij/MO_AR_SI_a;
  Fij=MO_AR_SI_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_AR_SI_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_AR_SI_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_AR_SI_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_AR_SI_EPS/rij*(0.35*SpExp(-0.3*rija)
			   
+0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}
//************************************************************************
double config::NeSi(particle* I, particle *J){
  rija=rij/MO_NE_SI_a;
  Fij=MO_NE_SI_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_NE_SI_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_NE_SI_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_NE_SI_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_NE_SI_EPS/rij*(0.35*SpExp(-0.3*rija)
			   
+0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}
//*************************************************************************
double config::HeF(particle* I, particle *J){
  rija=rij/MO_HE_F_a;
  Fij=MO_HE_F_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_HE_F_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_HE_F_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_HE_F_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_HE_F_EPS/rij*(0.35*SpExp(-0.3*rija)
			   +0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}
//*************************************************************************
double config::ArF(particle* I, particle *J){
  rija=rij/MO_AR_F_a;
  Fij=MO_AR_F_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_AR_F_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_AR_F_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_AR_F_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_AR_F_EPS/rij*(0.35*SpExp(-0.3*rija)
			   
+0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}
//************************************************************************
double config::NeF(particle* I, particle *J){
  rija=rij/MO_NE_F_a;
  Fij=MO_NE_F_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_NE_F_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_NE_F_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_NE_F_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_NE_F_EPS/rij*(0.35*SpExp(-0.3*rija)
			   
+0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}
//*************************************************************************
double config::HeCl(particle* I, particle *J){
  rija=rij/MO_HE_CL_a;
  Fij=MO_HE_CL_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_HE_CL_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_HE_CL_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_HE_CL_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_HE_CL_EPS/rij*(0.35*SpExp(-0.3*rija)
			   +0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}
//*************************************************************************
double config::ArCl(particle* I, particle *J){
  rija=rij/MO_AR_CL_a;
  Fij=MO_AR_CL_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_AR_CL_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_AR_CL_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_AR_CL_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_AR_CL_EPS/rij*(0.35*SpExp(-0.3*rija)
			   
+0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}
//************************************************************************
double config::NeCl(particle* I, particle *J){
  rija=rij/MO_NE_CL_a;
  Fij=MO_NE_CL_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_NE_CL_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_NE_CL_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_NE_CL_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_NE_CL_EPS/rij*(0.35*SpExp(-0.3*rija)
			   
+0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}
//*************************************************************************
double config::KrSi(particle* I, particle *J){
  rija=rij/MO_KR_SI_a;
  Fij=MO_KR_SI_EPS/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/MO_KR_SI_a)
			+0.55*SpExp(-1.2*rija)*(1/rij+1.2/MO_KR_SI_a)
			+0.1*SpExp(-6*rija)*(1/rij+6/MO_KR_SI_a))/rij*Rij;
  I->F += Fij; J->F -= Fij;
  return MO_KR_SI_EPS/rij*(0.35*SpExp(-0.3*rija)
			   
+0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
}

