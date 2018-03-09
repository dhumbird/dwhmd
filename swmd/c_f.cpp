#include "config.h"

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

//************************************************************************
double config::FF(particle* I, particle* J){
  rija=1/(rij-SW_F_F_a1);
  if (rija < 0){
    expa1ijk=SW_F_F_A*SpExp(SW_F_SIGMA*rija);
    Fij=-expa1ijk*
      (SW_F_F_SIGMA4/pow(rij,4)*(4/rij+SW_F_SIGMA*rija*rija)-
       SW_F_F_SIGMA8/pow(rij,8)*(8/rij+SW_F_SIGMA*rija*rija))/rij*Rij;
    
    I->F+=Fij; J->F-=Fij;
    return expa1ijk*
      (SW_F_F_SIGMA8/pow(rij,8)-SW_F_F_SIGMA4/pow(rij,4));
  }
  else return 0;
}
//************************************************************************
double config::FFF(particle* I, particle* J, particle* K){
  rija=1/(rij-SW_F_F_a1); 
  rika=1/(rik-SW_F_F_a1); 
  rjka=1/(rjk-SW_F_F_a1);
  rija2=1/(rij-SW_F_F_a2); 
  rika2=1/(rik-SW_F_F_a2); 
  rjka2=1/(rjk-SW_F_F_a2);
  
  cosOjik2=SW_F_F_F_LAMBDA2-cosOjik*cosOjik*SW_F_F_F_LAMBDA3;
  cosOijk2=SW_F_F_F_LAMBDA2-cosOijk*cosOijk*SW_F_F_F_LAMBDA3;
  cosOikj2=SW_F_F_F_LAMBDA2-cosOikj*cosOikj*SW_F_F_F_LAMBDA3;

  if (rija<0 && rika<0){ // jik
    expa1jik=SpExp(SW_F_GAM1SIG*(rija+rika));
    expa2jik=(rija2<0 && rika2<0)?
      SpExp(SW_F_GAM2SIG*(rija2+rika2)):0;
    hjik2=SW_F_F_F_LAMBDA1*pow(SW_F_F_SIGMA2/rij/rik,4)*expa1jik;
    hjik=hjik2+cosOjik2*expa2jik;
  }
  else{expa1jik=expa2jik=hjik=hjik2=0;}
  
  if (rija<0 && rjka<0){ // ijk
    expa1ijk=SpExp(SW_F_GAM1SIG*(rija+rjka));
    expa2ijk=(rija2<0 && rjka2<0)?
      SpExp(SW_F_GAM2SIG*(rija2+rjka2)):0;
    hijk2=SW_F_F_F_LAMBDA1*pow(SW_F_F_SIGMA2/rij/rjk,4)*expa1ijk;
    hijk=hijk2+cosOijk2*expa2ijk;
  }
  else{expa1ijk=expa2ijk=hijk=hijk2=0;}
  
  if (rika<0 && rjka<0){ // ikj
    expa1ikj=SpExp(SW_F_GAM1SIG*(rika+rjka));
    expa2ikj=(rika2<0 && rjka2<0)?
      SpExp(SW_F_GAM2SIG*(rika2+rjka2)):0;
    hikj2=SW_F_F_F_LAMBDA1*pow(SW_F_F_SIGMA2/rik/rjk,4)*expa1ikj;
    hikj=hikj2+cosOikj2*expa2ikj;
  }
  else{expa1ikj=expa2ikj=hikj=hikj2=0;}
  
  rija*=rija*SW_F_GAM1SIG;
  rika*=rika*SW_F_GAM1SIG;
  rjka*=rjka*SW_F_GAM1SIG;
  rija2*=rija2*SW_F_GAM2SIG;
  rika2*=rika2*SW_F_GAM2SIG;
  rjka2*=rjka2*SW_F_GAM2SIG;

  I->F+=
    hjik2*((4/rij + rija)*Ruij+(4/rik + rika)*Ruik)
    +expa2jik*(2*SW_F_F_F_LAMBDA3*cosOjik*
	       ((1/rik-cosOjik/rij)*Ruij+(1/rij-cosOjik/rik)*Ruik)
	       +cosOjik2*(Ruij*rija2+Ruik*rika2))
    +hijk2*(4/rij+rija)*Ruij
    +expa2ijk*(2*SW_F_F_F_LAMBDA3*cosOijk*
	       ((1/rjk-cosOijk/rij)*Ruij-(rik/rjk/rij)*Ruik)
	       +(cosOijk2*rija2)*Ruij)
    +hikj2*(4/rik+rika)*Ruik
    +expa2ikj*(2*SW_F_F_F_LAMBDA3*cosOikj*
	       ((1/rjk-cosOikj/rik)*Ruik-(rij/rjk/rik)*Ruij)
	       +(cosOikj2*rika2)*Ruik);
  
  J->F-=
    hjik2*(4/rij+rija)*Ruij
    -expa2jik*(2*SW_F_F_F_LAMBDA3*cosOjik*
	       ((cosOjik/rij-1/rik)*Ruij-(rjk/rij/rik)*Rujk)
	       -(cosOjik2*rija2)*Ruij)
    +hijk2*((4/rij+rija)*Ruij-(4/rjk+rjka)*Rujk)
    -expa2ijk*(2*SW_F_F_F_LAMBDA3*cosOijk*
	       ((cosOijk/rij-1/rjk)*Ruij+(1/rij-cosOijk/rjk)*Rujk)
	       -cosOijk2*(Ruij*rija2-Rujk*rjka2))
    -hikj2*(4/rjk+rjka)*Rujk
    -expa2ikj*(2*SW_F_F_F_LAMBDA3*cosOikj*
	       ((1/rik-cosOikj/rjk)*Rujk+(rij/rjk/rik)*Ruij)
	       +(cosOikj2*rjka2)*Rujk);
  
  K->F-=
    hjik2*(4/rik+rika)*Ruik
    -expa2jik*(2*SW_F_F_F_LAMBDA3*cosOjik*
	       ((cosOjik/rik-1/rij)*Ruik+(rjk/rij/rik)*Rujk)
	       -(cosOjik2*rika2)*Ruik)
    +hijk2*(4/rjk+rjka)*Rujk
    -expa2ijk*(2*SW_F_F_F_LAMBDA3*cosOijk*
	       ((cosOijk/rjk-1/rij)*Rujk+(rik/rij/rjk)*Ruik)
	       -(cosOijk2*rjka2)*Rujk)
    +hikj2*((4/rik+rika)*Ruik+(4/rjk+rjka)*Rujk)
    -expa2ikj*(2*SW_F_F_F_LAMBDA3*cosOikj*
	       ((cosOikj/rik-1/rjk)*Ruik+(cosOikj/rjk-1/rik)*Rujk)
	       -cosOikj2*(Ruik*rika2+Rujk*rjka2));
  
  return hjik+hijk+hikj;
}
