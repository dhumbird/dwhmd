#include "config.h"
//**********************************************************************
double config::u_SiSi(){
  rija=1/(rij-SW_SI_a);
  if (rija < 0){
    return SW_SI_SI_EPSA*SpExp(SW_SI_SIGMA*rija)*
      (SW_SI_B_SIGMA4/pow(rij,4)-1);
  }
  else return 0;
}
//*********************************************************************
double config::u_SiSiSi(){
  rija=1/(rij-SW_SI_a); 
  rika=1/(rik-SW_SI_a); 
  rjka=1/(rjk-SW_SI_a);
  
  if(rija<0 && rika<0){
    hjik=SpExp(SW_SI_SI_SI_GAMMA*(rija+rika))
      *pow(cosOjik+ONE_THIRD,2);
  }
  else hjik=0;
  
  if(rija<0 && rjka<0){
    hijk=SpExp(SW_SI_SI_SI_GAMMA*(rija+rjka))
      *pow(cosOijk+ONE_THIRD,2);
  }
  else hijk=0;
  
  if(rika<0 && rjka<0){
    hikj2=SpExp(SW_SI_SI_SI_GAMMA*(rika+rjka))
      *pow(cosOikj+ONE_THIRD,2);
  }
  else hikj=0;
 
  return SW_SI_SI_SI_LAMBDA*SW_SI_SI_EPS*(hjik+hijk+hikj);
}


#define SW_F_SIGMA 1.2141
#define SW_F_EPS 1.66
#define SW_F_F_F_GAMMA1 1
#define SW_F_F_F_GAMMA2 3
const double SW_F_F_a1 = 3.6 * SW_F_SIGMA;
const double SW_F_F_a2 = 2.8 * SW_F_SIGMA;
const double SW_F_F_SIGMA8 = pow(SW_F_SIGMA, 8);
const double SW_F_F_SIGMA4 = pow(SW_F_SIGMA, 4);
const double SW_F_F_SIGMA2 = pow(SW_F_SIGMA, 2);
const double SW_F_F_A = 6.052463017 * SW_F_EPS;
const double SW_F_F_F_LAMBDA1 = SW_F_EPS * 8.4;
const double SW_F_F_F_LAMBDA2 = SW_F_EPS * 50;
const double SW_F_F_F_LAMBDA3 = SW_F_EPS * 25;
const double SW_F_GAM1SIG = SW_F_F_F_GAMMA1 * SW_F_SIGMA;
const double SW_F_GAM2SIG = SW_F_F_F_GAMMA2 * SW_F_SIGMA;
//****** WWC-SW F/Si *************************
const double SW_SI_F_a = 1.8 * SW_SI_SIGMA;

const double SW_SI_F_A = 21.199221 * SW_SI_SI_EPS;
#define SW_SI_F_B 0.546418
#define SW_SI_F_GAMMA 1.339450
const double SW_SI_F_GAMSIG = SW_SI_F_GAMMA * SW_SI_SIGMA;

const double SW_SI_F_AAB_A = 3.624533 * SW_SI_SI_EPS;
const double SW_SI_SI_F_LAMBDA = 0.218615 * SW_SI_F_AAB_A;
#define SW_SI_SI_F_GAMMA 0.463088

const double SW_SI_F_ABA_A = 50.874092 * SW_SI_SI_EPS;
#define SW_SI_F_SI_GAMMA 1.371580

const double SW_SI_F_ABB_A = 2.792073 * SW_SI_SI_EPS;

const double SW_SI_F_BAB_A = 22.406434 * SW_SI_SI_EPS;
const double SW_F_SI_F_LAMBDA = 2.068601 * SW_SI_SI_EPS;
#define SW_F_SI_F_GAMMA 0.890132


//************************************************************************
double config::u_FF(){
  rija=1/(rij-SW_F_F_a1);
  if (rija < 0){
    return SW_F_F_A*SpExp(SW_F_SIGMA*rija)*
      (SW_F_F_SIGMA8/pow(rij,8)-SW_F_F_SIGMA4/pow(rij,4));
  }
  else return 0;
}
//************************************************************************
double config::u_FFF(){
  rija=1/(rij-SW_F_F_a1); 
  rika=1/(rik-SW_F_F_a1); 
  rjka=1/(rjk-SW_F_F_a1);
  rija2=1/(rij-SW_F_F_a2); 
  rika2=1/(rik-SW_F_F_a2); 
  rjka2=1/(rjk-SW_F_F_a2);

  if (rija<0 && rika<0){ // jik
    expa2jik=(rija2<0 && rika2<0)?
      SpExp(SW_F_GAM2SIG*(rija2+rika2)):0;
    hjik=SW_F_F_F_LAMBDA1*pow(SW_F_F_SIGMA2/rij/rik,4)*SpExp(SW_F_GAM1SIG*(rija+rika))+(SW_F_F_F_LAMBDA2-cosOjik*cosOjik*SW_F_F_F_LAMBDA3)*expa2jik;
  }
  else hjik=0;
  
  if (rija<0 && rjka<0){ // ijk
    expa2ijk=(rija2<0 && rjka2<0)?
      SpExp(SW_F_GAM2SIG*(rija2+rjka2)):0;
    hijk=SW_F_F_F_LAMBDA1*pow(SW_F_F_SIGMA2/rij/rjk,4)*SpExp(SW_F_GAM1SIG*(rija+rjka))+(SW_F_F_F_LAMBDA2-cosOijk*cosOijk*SW_F_F_F_LAMBDA3)*expa2ijk;
  }
  else hijk=0;
  
  if (rika<0 && rjka<0){ // ikj
    expa2ikj=(rika2<0 && rjka2<0)?
      SpExp(SW_F_GAM2SIG*(rika2+rjka2)):0;
    hikj=SW_F_F_F_LAMBDA1*pow(SW_F_F_SIGMA2/rik/rjk,4)*SpExp(SW_F_GAM1SIG*(rika+rjka))+(SW_F_F_F_LAMBDA2-cosOikj*cosOikj*SW_F_F_F_LAMBDA3)*expa2ikj;
  }
  else hikj=0;
  
  return hjik+hijk+hikj;
}

//************************************************************************
double config::u_FSi(){
  rija=1/(rij-SW_SI_F_a);
  if (rija < 0){
    return SW_SI_F_A*(SW_SI_F_B*pow(SW_SI_SIGMA/rij,3)-pow(SW_SI_SIGMA/rij,2))
      *SpExp(SW_SI_F_GAMSIG*rija);
  }
  else return 0;
}
//************************************************************************
double config::u_FSiSi(particle* I, particle* J, particle* K){

  rija=1/(rij-SW_SI_F_a); 
  rika=1/(rik-SW_SI_F_a); 
  rjka=1/(rjk-SW_SI_F_a);

  if(rija<0 && rika<0)
    if (I->id==9)
      hjik=SW_SI_F_ABA_A*SpExp(SW_SI_F_SI_GAMMA*SW_SI_SIGMA*(rija+rika));
    else
      hjik=(SW_SI_F_AAB_A*pow(cosOjik+ONE_THIRD,2)-SW_SI_SI_F_LAMBDA)*
	SpExp(SW_SI_SI_F_GAMMA*SW_SI_SIGMA*(rija+rika));
  else hjik=0;

  if(rija<0 && rjka<0)
    if (J->id==9)
      hijk=SW_SI_F_ABA_A*SpExp(SW_SI_F_SI_GAMMA*SW_SI_SIGMA*(rija+rjka));
    else
      hijk=(SW_SI_F_AAB_A*pow(cosOijk+ONE_THIRD,2)-SW_SI_SI_F_LAMBDA)*
	SpExp(SW_SI_SI_F_GAMMA*SW_SI_SIGMA*(rija+rjka));
  else hijk=0;
  
  if(rika<0 && rjka<0)
    if (K->id==9)
      hikj=SW_SI_F_ABA_A*SpExp(SW_SI_F_SI_GAMMA*SW_SI_SIGMA*(rika+rjka));
    else
      hikj=(SW_SI_F_AAB_A*pow(cosOikj+ONE_THIRD,2)-SW_SI_SI_F_LAMBDA)*
	SpExp(SW_SI_SI_F_GAMMA*SW_SI_SIGMA*(rjka+rika));
  else hikj=0;
  
  return hjik+hijk+hikj;
}
//************************************************************************
double config::u_FFSi(particle* I, particle* J, particle* K){

  rija=1/(rij-SW_SI_F_a); 
  rika=1/(rik-SW_SI_F_a); 
  rjka=1/(rjk-SW_SI_F_a);

  if(rija<0 && rika<0)
    if (I->id==14)
      hjik=(SW_SI_F_BAB_A*pow(cosOjik-COS103,2)-SW_F_SI_F_LAMBDA)*
	SpExp(SW_F_SI_F_GAMMA*SW_SI_SIGMA*(rija+rika));
    else
      hjik=SW_SI_F_ABB_A*SpExp(SW_SI_SIGMA*(rija+rika));
  else hjik=0;
  
  if(rija<0 && rjka<0)
    if (J->id==14)
      hijk=(SW_SI_F_BAB_A*pow(cosOijk-COS103,2)-SW_F_SI_F_LAMBDA)*
	SpExp(SW_F_SI_F_GAMMA*SW_SI_SIGMA*(rija+rjka));
    else
      hijk=SW_SI_F_ABB_A*SpExp(SW_SI_SIGMA*(rija+rjka));
  else hijk=0;
  
  if(rika<0 && rjka<0)
    if (K->id==14)
      hikj=(SW_SI_F_BAB_A*pow(cosOikj-COS103,2)-SW_F_SI_F_LAMBDA)*
	SpExp(SW_F_SI_F_GAMMA*SW_SI_SIGMA*(rjka+rika));
    else
      hikj=SW_SI_F_ABB_A*SpExp(SW_SI_SIGMA*(rjka+rika));
  else hikj=0;
  
  return hjik+hijk+hikj;
}

