#include "config.h"

#define SW_SI_CL_B 0.67
#define SW_SI_CL_GAMMA 1
const double SW_SI_CL_GAMSIG = SW_SI_CL_GAMMA * SW_SI_SIGMA;
const double SW_SI_CL_A = 28 * SW_SI_SI_EPS;
const double SW_SI_CL_a = 1.8 * SW_SI_SIGMA;
const double SW_SI_CL_AAB_A = 15 * SW_SI_SI_EPS;
const double SW_SI_CL_ABA_A = 50 * SW_SI_SI_EPS;
const double SW_SI_CL_ABB_A = 15 * SW_SI_SI_EPS;
const double SW_SI_CL_BAB_A = 30 * SW_SI_SI_EPS;
const double SW_SI_CL_LAMBDA = 0.5 * SW_SI_SI_EPS;

//************************************************************************
double config::ClSi(particle* I, particle* J){
  rija=1/(rij-SW_SI_CL_a);
  if (rija < 0){
    expa1ijk=SpExp(1.3*SW_SI_SIGMA*rija);
    hijk2=SW_SI_CL_A*(SW_SI_CL_B*pow(SW_SI_SIGMA/rij,2.2)
		      -pow(SW_SI_SIGMA/rij,0.9))*expa1ijk;
    Fij=(-SW_SI_CL_A*expa1ijk*
	 (0.9*pow(SW_SI_SIGMA/rij,0.9)/rij
	  -2.2*SW_SI_CL_B*pow(SW_SI_SIGMA/rij,2.2)/rij)
	 +1.3*SW_SI_SIGMA*rija*rija*hijk2)/rij*Rij;
    I->F+=Fij; J->F-=Fij;
    return hijk2;
  }
  else return 0;
}
//************************************************************************
double config::ClSiSi(particle* I, particle* J, particle* K){

  rija=1/(rij-SW_SI_CL_a); 
  rika=1/(rik-SW_SI_CL_a); 
  rjka=1/(rjk-SW_SI_CL_a);
  rija2=SW_SI_SIGMA*rija*rija;
  rika2=SW_SI_SIGMA*rika*rika;
  rjka2=SW_SI_SIGMA*rjka*rjka;

  if(rija<0 && rika<0){
    if (I->id==17){
      hjik=SW_SI_CL_ABA_A*SpExp(SW_SI_CL_GAMSIG*(rija+rika));
      I->F+=(hjik*SW_SI_CL_GAMMA)*(rija2*Ruij+rika2*Ruik);
      J->F-=(hjik*SW_SI_CL_GAMMA*rija2)*Ruij;
      K->F-=(hjik*SW_SI_CL_GAMMA*rika2)*Ruik;
    }    
    else{
      cosOjik2=cosOjik+ONE_THIRD;
      hjik=SW_SI_CL_AAB_A*cosOjik2*cosOjik2*SpExp(SW_SI_SIGMA*(rija+rika));
 
      I->F-=hjik*(2/cosOjik2*
		  ((1/rij-cosOjik/rik)*Ruik+(1/rik-cosOjik/rij)*Ruij)
		  -(rija2*Ruij+rika2*Ruik));
      
      J->F-=hjik*(2/cosOjik2*
		  ((cosOjik/rij-1/rik)*Ruij-(rjk/rik/rij)*Rujk)+(rija2*Ruij));
		  
      K->F-=hjik*(2/cosOjik2*
		  ((cosOjik/rik-1/rij)*Ruik+(rjk/rij/rik)*Rujk)+(rika2*Ruik));
    }
  }
  else hjik=0;

  if(rija<0 && rjka<0){
    if (J->id==17){
      hijk=SW_SI_CL_ABA_A*SpExp(SW_SI_CL_GAMSIG*(rija+rjka));
      I->F+=(hijk*SW_SI_CL_GAMMA*rija2)*Ruij;
      J->F-=(hijk*SW_SI_CL_GAMMA)*(rija2*Ruij-rjka2*Rujk);
      K->F-=(hijk*SW_SI_CL_GAMMA*rjka2)*Rujk;
    }
    else{
      cosOijk2=cosOijk+ONE_THIRD;
      hijk=SW_SI_CL_AAB_A*cosOijk2*cosOijk2*
	SpExp(SW_SI_SIGMA*(rija+rjka));
      
      I->F-=hijk*(2/cosOijk2*
		  ((1/rjk-cosOijk/rij)*Ruij-(rik/rij/rjk)*Ruik)-(rija2*Ruij));

      J->F-=hijk*(2/cosOijk2*
		  ((cosOijk/rij-1/rjk)*Ruij+(1/rij-cosOijk/rjk)*Rujk)
		  +(rija2*Ruij-rjka2*Rujk));
      
      K->F-=hijk*(2/cosOijk2*
		  ((cosOijk/rjk-1/rij)*Rujk+(rik/rij/rjk)*Ruik)+(rjka2*Rujk));
    }
  }
  else hijk=0;
  
  if(rika<0 && rjka<0){
    if (K->id==17){
      hikj=SW_SI_CL_ABA_A*SpExp(SW_SI_CL_GAMSIG*(rika+rjka));
      I->F+=(hikj*SW_SI_CL_GAMMA*rika2)*Ruik;
      J->F+=(hikj*SW_SI_CL_GAMMA*rjka2)*Rujk;
      K->F-=(hikj*SW_SI_CL_GAMMA)*(rjka2*Rujk+rika2*Ruik);
    }
    else{
      cosOikj2=cosOikj+ONE_THIRD;
      hikj=SW_SI_CL_AAB_A*cosOikj2*cosOikj2*
	SpExp(SW_SI_SIGMA*(rjka+rika));
      
      I->F-=hikj*(2/cosOikj2*
		  ((1/rjk-cosOikj/rik)*Ruik-(rij/rik/rjk)*Ruij)-(rika2*Ruik));

      J->F-=hikj*(2/cosOikj2*
		  ((1/rik-cosOikj/rjk)*Rujk+(rij/rik/rjk)*Ruij)-(rjka2*Rujk));

      K->F-=hikj*(2/cosOikj2*
		  ((cosOikj/rik-1/rjk)*Ruik+(cosOikj/rjk-1/rik)*Rujk)
		  +(rika2*Ruik+rjka2*Rujk));
    }
  }
  else hikj=0;
  return hjik+hijk+hikj;
}
//************************************************************************
double config::ClClSi(particle* I, particle* J, particle* K){

  rija=1/(rij-SW_SI_CL_a); 
  rika=1/(rik-SW_SI_CL_a); 
  rjka=1/(rjk-SW_SI_CL_a);
  rija2=SW_SI_SIGMA*rija*rija;
  rika2=SW_SI_SIGMA*rika*rika;
  rjka2=SW_SI_SIGMA*rjka*rjka;

  if(rija<0 && rika<0){
    
    expa1jik=SpExp(SW_SI_SIGMA*(rija+rika));
    if (I->id==14){
      cosOjik2=cosOjik-COS103;
      hjik=(SW_SI_CL_BAB_A*cosOjik2*cosOjik2-SW_SI_CL_LAMBDA)*expa1jik;
      expa2jik=expa1jik*2*SW_SI_CL_BAB_A*cosOjik2;
      I->F-=expa2jik*((1/rij-cosOjik/rik)*Ruik+(1/rik-cosOjik/rij)*Ruij)
	-hjik*(rija2*Ruij+rika2*Ruik);

      J->F-=expa2jik*((cosOjik/rij-1/rik)*Ruij-(rjk/rik/rij)*Rujk)
	+(hjik*rija2)*Ruij;

      K->F-=expa2jik*((cosOjik/rik-1/rij)*Ruik+(rjk/rij/rik)*Rujk)
	+(hjik*rika2)*Ruik;
    }    
    else{
      hjik=SW_SI_CL_ABB_A*expa1jik;
      I->F+=hjik*(rija2*Ruij+rika2*Ruik);
      J->F-=(hjik*rija2)*Ruij;
      K->F-=(hjik*rika2)*Ruik;
    }
  }
  else hjik=0;
  
  if(rija<0 && rjka<0){
    expa1ijk=SpExp(SW_SI_SIGMA*(rija+rjka));
    if (J->id==14){
      cosOijk2=cosOijk-COS103;
      hijk=(SW_SI_CL_BAB_A*cosOijk2*cosOijk2-SW_SI_CL_LAMBDA)*expa1ijk;
      expa2ijk=expa1ijk*2*SW_SI_CL_BAB_A*cosOijk2;
      
      I->F-=expa2ijk*((1/rjk-cosOijk/rij)*Ruij-(rik/rij/rjk)*Ruik)
	-(hijk*rija2)*Ruij;
      
      J->F-=expa2ijk*((cosOijk/rij-1/rjk)*Ruij+(1/rij-cosOijk/rjk)*Rujk)
	+hijk*(rija2*Ruij-rjka2*Rujk);
      
      K->F-=expa2ijk*((cosOijk/rjk-1/rij)*Rujk+(rik/rij/rjk)*Ruik)
	+(hijk*rjka2)*Rujk;
    }
    else{
      hijk=SW_SI_CL_ABB_A*expa1ijk;
      I->F+=(hijk*rija2)*Ruij;
      J->F-=hijk*(rija2*Ruij-rjka2*Rujk);
      K->F-=(hijk*rjka2)*Rujk;
    }
  }
  else hijk=0;
  
  if(rika<0 && rjka<0){
    expa1ikj=SpExp(SW_SI_SIGMA*(rjka+rika));
    if (K->id==14){
      cosOikj2=cosOikj-COS103;
      hikj=(SW_SI_CL_BAB_A*cosOikj2*cosOikj2-SW_SI_CL_LAMBDA)*expa1ikj;
      expa2ikj=expa1ikj*2*SW_SI_CL_BAB_A*cosOikj2;
      I->F-=expa2ikj*((1/rjk-cosOikj/rik)*Ruik-(rij/rik/rjk)*Ruij)
	-(hikj*rika2)*Ruik;

      J->F-=expa2ikj*((1/rik-cosOikj/rjk)*Rujk+(rij/rik/rjk)*Ruij)
	-(hikj*rjka2)*Rujk;
      
      K->F-=expa2ikj*((cosOikj/rik-1/rjk)*Ruik+(cosOikj/rjk-1/rik)*Rujk)
	+hikj*(rika2*Ruik+rjka2*Rujk);
    }
    else{
      hikj=SW_SI_CL_ABB_A*expa1ikj;
      I->F+=(hikj*rika2)*Ruik;
      J->F+=(hikj*rjka2)*Rujk;
      K->F-=hikj*(rjka2*Rujk+rika2*Ruik);  
    }
  }
  else hikj=0;
  return hjik+hijk+hikj;
}
