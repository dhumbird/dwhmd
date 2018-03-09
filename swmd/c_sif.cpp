#include "config.h"

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
double config::FSi(particle* I, particle* J){
  rija=1/(rij-SW_SI_F_a);
  if (rija < 0){
    expa1ijk=SpExp(SW_SI_F_GAMSIG*rija);
    hijk2=SW_SI_F_A*(SW_SI_F_B*pow(SW_SI_SIGMA/rij,3)-pow(SW_SI_SIGMA/rij,2))
      *expa1ijk;
    Fij=(-SW_SI_F_A*expa1ijk*
	 (2*SW_SI_SIGMA2/pow(rij,3)-3*SW_SI_F_B*SW_SI_SIGMA3/pow(rij,4))
	 +SW_SI_F_GAMSIG*rija*rija*hijk2)/rij*Rij;
    I->F+=Fij; J->F-=Fij;
    return hijk2;
  }
  else return 0;
}
//************************************************************************
double config::FSiSi(particle* I, particle* J, particle* K){

  rija=1/(rij-SW_SI_F_a); 
  rika=1/(rik-SW_SI_F_a); 
  rjka=1/(rjk-SW_SI_F_a);
  rija2=SW_SI_SIGMA*rija*rija;
  rika2=SW_SI_SIGMA*rika*rika;
  rjka2=SW_SI_SIGMA*rjka*rjka;

  if(rija<0 && rika<0){
    if (I->id==9){ //ABA
      hjik=SW_SI_F_ABA_A*SpExp(SW_SI_F_SI_GAMMA*SW_SI_SIGMA*(rija+rika));
      I->F+=(hjik*SW_SI_F_SI_GAMMA)*(rija2*Ruij+rika2*Ruik);
      J->F-=(hjik*SW_SI_F_SI_GAMMA*rija2)*Ruij;
      K->F-=(hjik*SW_SI_F_SI_GAMMA*rika2)*Ruik;
    }    
    else{ //AAB
      expa1jik=SpExp(SW_SI_SI_F_GAMMA*SW_SI_SIGMA*(rija+rika));
      cosOjik2=cosOjik+ONE_THIRD;
      hjik=(SW_SI_F_AAB_A*cosOjik2*cosOjik2-SW_SI_SI_F_LAMBDA)*expa1jik;
      expa1jik*=(2*SW_SI_F_AAB_A*cosOjik2);
      
      I->F-=expa1jik*((1/rij-cosOjik/rik)*Ruik+(1/rik-cosOjik/rij)*Ruij)
        -hjik*SW_SI_SI_F_GAMMA*(rija2*Ruij+rika2*Ruik);

      J->F-=expa1jik*((cosOjik/rij-1/rik)*Ruij-(rjk/rik/rij)*Rujk)
        +(hjik*SW_SI_SI_F_GAMMA*rija2)*Ruij;

      K->F-=expa1jik*((cosOjik/rik-1/rij)*Ruik+(rjk/rij/rik)*Rujk)
        +(hjik*SW_SI_SI_F_GAMMA*rika2)*Ruik;
    }
  }
  else hjik=0;

  if(rija<0 && rjka<0){
    if (J->id==9){ //ABA
      hijk=SW_SI_F_ABA_A*SpExp(SW_SI_F_SI_GAMMA*SW_SI_SIGMA*(rija+rjka));
      I->F+=(hijk*SW_SI_F_SI_GAMMA*rija2)*Ruij;
      J->F-=(hijk*SW_SI_F_SI_GAMMA)*(rija2*Ruij-rjka2*Rujk);
      K->F-=(hijk*SW_SI_F_SI_GAMMA*rjka2)*Rujk;
    }
    else{ //AAB
      expa1ijk=SpExp(SW_SI_SI_F_GAMMA*SW_SI_SIGMA*(rija+rjka));
      cosOijk2=cosOijk+ONE_THIRD;
      hijk=(SW_SI_F_AAB_A*cosOijk2*cosOijk2-SW_SI_SI_F_LAMBDA)*expa1ijk;
      expa1ijk*=(2*SW_SI_F_AAB_A*cosOijk2);
      
      I->F-=expa1ijk*((1/rjk-cosOijk/rij)*Ruij-(rik/rij/rjk)*Ruik)
        -(hijk*SW_SI_SI_F_GAMMA*rija2)*Ruij;
      
      J->F-=expa1ijk*((cosOijk/rij-1/rjk)*Ruij+(1/rij-cosOijk/rjk)*Rujk)
        +hijk*SW_SI_SI_F_GAMMA*(rija2*Ruij-rjka2*Rujk);
      
      K->F-=expa1ijk*((cosOijk/rjk-1/rij)*Rujk+(rik/rij/rjk)*Ruik)
        +(hijk*SW_SI_SI_F_GAMMA*rjka2)*Rujk;
    }
  }
  else hijk=0;
  
  if(rika<0 && rjka<0){
    if (K->id==9){ //ABA
      hikj=SW_SI_F_ABA_A*SpExp(SW_SI_F_SI_GAMMA*SW_SI_SIGMA*(rika+rjka));
      I->F+=(hikj*SW_SI_F_SI_GAMMA*rika2)*Ruik;
      J->F+=(hikj*SW_SI_F_SI_GAMMA*rjka2)*Rujk;
      K->F-=(hikj*SW_SI_F_SI_GAMMA)*(rjka2*Rujk+rika2*Ruik);
    }
    else{ //AAB
      expa1ikj=SpExp(SW_SI_SI_F_GAMMA*SW_SI_SIGMA*(rjka+rika));
      cosOikj2=cosOikj+ONE_THIRD;
      hikj=(SW_SI_F_AAB_A*cosOikj2*cosOikj2-SW_SI_SI_F_LAMBDA)*expa1ikj;
      expa1ikj*=(2*SW_SI_F_AAB_A*cosOikj2); 
      
      I->F-=expa1ikj*((1/rjk-cosOikj/rik)*Ruik-(rij/rik/rjk)*Ruij)
        -(hikj*SW_SI_SI_F_GAMMA*rika2)*Ruik;

      J->F-=expa1ikj*((1/rik-cosOikj/rjk)*Rujk+(rij/rik/rjk)*Ruij)
        -(hikj*SW_SI_SI_F_GAMMA*rjka2)*Rujk;
      
      K->F-=expa1ikj*((cosOikj/rik-1/rjk)*Ruik+(cosOikj/rjk-1/rik)*Rujk)
        +hikj*SW_SI_SI_F_GAMMA*(rika2*Ruik+rjka2*Rujk);
    }
  }
  else hikj=0;
  return hjik+hijk+hikj;
}
//************************************************************************
double config::FFSi(particle* I, particle* J, particle* K){

  rija=1/(rij-SW_SI_F_a); 
  rika=1/(rik-SW_SI_F_a); 
  rjka=1/(rjk-SW_SI_F_a);
  rija2=SW_SI_SIGMA*rija*rija;
  rika2=SW_SI_SIGMA*rika*rika;
  rjka2=SW_SI_SIGMA*rjka*rjka;

  if(rija<0 && rika<0){
    if (I->id==14){ //BAB
      expa1jik=SpExp(SW_F_SI_F_GAMMA*SW_SI_SIGMA*(rija+rika));
      cosOjik2=cosOjik-COS103;
      hjik=(SW_SI_F_BAB_A*cosOjik2*cosOjik2-SW_F_SI_F_LAMBDA)*expa1jik;
      expa1jik*=(2*SW_SI_F_BAB_A*cosOjik2);
      
      I->F-=expa1jik*((1/rij-cosOjik/rik)*Ruik+(1/rik-cosOjik/rij)*Ruij)
	-hjik*SW_F_SI_F_GAMMA*(rija2*Ruij+rika2*Ruik);

      J->F-=expa1jik*((cosOjik/rij-1/rik)*Ruij-(rjk/rik/rij)*Rujk)
	+(hjik*SW_F_SI_F_GAMMA*rija2)*Ruij;

      K->F-=expa1jik*((cosOjik/rik-1/rij)*Ruik+(rjk/rij/rik)*Rujk)
	+(hjik*SW_F_SI_F_GAMMA*rika2)*Ruik;
    }    
    else{ //ABB
      hjik=SW_SI_F_ABB_A*SpExp(SW_SI_SIGMA*(rija+rika));
      I->F+=hjik*(rija2*Ruij+rika2*Ruik);
      J->F-=(hjik*rija2)*Ruij;
      K->F-=(hjik*rika2)*Ruik;
    }
  }
  else hjik=0;
  
  if(rija<0 && rjka<0){
    if (J->id==14){ //BAB
      expa1ijk=SpExp(SW_F_SI_F_GAMMA*SW_SI_SIGMA*(rija+rjka));
      cosOijk2=cosOijk-COS103;
      hijk=(SW_SI_F_BAB_A*cosOijk2*cosOijk2-SW_F_SI_F_LAMBDA)*expa1ijk;
      expa1ijk*=(2*SW_SI_F_BAB_A*cosOijk2);
      
      I->F-=expa1ijk*((1/rjk-cosOijk/rij)*Ruij-(rik/rij/rjk)*Ruik)
	-(hijk*SW_F_SI_F_GAMMA*rija2)*Ruij;
      
      J->F-=expa1ijk*((cosOijk/rij-1/rjk)*Ruij+(1/rij-cosOijk/rjk)*Rujk)
	+hijk*SW_F_SI_F_GAMMA*(rija2*Ruij-rjka2*Rujk);
      
      K->F-=expa1ijk*((cosOijk/rjk-1/rij)*Rujk+(rik/rij/rjk)*Ruik)
	+(hijk*SW_F_SI_F_GAMMA*rjka2)*Rujk;
    }
    else{ //ABB
      hijk=SW_SI_F_ABB_A*SpExp(SW_SI_SIGMA*(rija+rjka));
      I->F+=(hijk*rija2)*Ruij;
      J->F-=hijk*(rija2*Ruij-rjka2*Rujk);
      K->F-=(hijk*rjka2)*Rujk;
    }
  }
  else hijk=0;
  
  if(rika<0 && rjka<0){
    if (K->id==14){ //BAB
      expa1ikj=SpExp(SW_F_SI_F_GAMMA*SW_SI_SIGMA*(rjka+rika));
      cosOikj2=cosOikj-COS103;
      hikj=(SW_SI_F_BAB_A*cosOikj2*cosOikj2-SW_F_SI_F_LAMBDA)*expa1ikj;
      expa1ikj*=(2*SW_SI_F_BAB_A*cosOikj2);
      
      I->F-=expa1ikj*((1/rjk-cosOikj/rik)*Ruik-(rij/rik/rjk)*Ruij)
	-(hikj*SW_F_SI_F_GAMMA*rika2)*Ruik;

      J->F-=expa1ikj*((1/rik-cosOikj/rjk)*Rujk+(rij/rik/rjk)*Ruij)
	-(hikj*SW_F_SI_F_GAMMA*rjka2)*Rujk;
      
      K->F-=expa1ikj*((cosOikj/rik-1/rjk)*Ruik+(cosOikj/rjk-1/rik)*Rujk)
	+hikj*SW_F_SI_F_GAMMA*(rika2*Ruik+rjka2*Rujk);
    }
    else{ //ABB
      hikj=SW_SI_F_ABB_A*SpExp(SW_SI_SIGMA*(rjka+rika));
      I->F+=(hikj*rika2)*Ruik;
      J->F+=(hikj*rjka2)*Rujk;
      K->F-=hikj*(rjka2*Rujk+rika2*Ruik);  
    }
  }
  else hikj=0;
  return hjik+hijk+hikj;
}
