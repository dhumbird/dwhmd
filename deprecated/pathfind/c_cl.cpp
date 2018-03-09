#include "config.h"

const double SW_CL_SIGMA = SW_SI_SIGMA;
const double SW_CL_EPS = SW_SI_SI_EPS;
#define SW_CL_CL_CL_GAMMA1 0.5795
#define SW_CL_CL_CL_GAMMA2 1.7386
const double SW_CL_CL_a1 = 2.0862 * SW_CL_SIGMA;
const double SW_CL_CL_a2 = 1.6226 * SW_CL_SIGMA;
const double SW_CL_CL_SIGMA6 = pow(SW_CL_SIGMA, 6);
const double SW_CL_CL_SIGMA5 = pow(SW_CL_SIGMA, 5);
const double SW_CL_CL_SIGMA2 = pow(SW_CL_SIGMA, 2);
const double SW_CL_CL_A = 8.611 * SW_CL_EPS;
#define SW_CL_CL_B 0.789
const double SW_CL_CL_CL_LAMBDA1 = SW_CL_EPS * 3;
const double SW_CL_CL_CL_LAMBDA2 = SW_CL_EPS * 2 * 23.778;
const double SW_CL_CL_CL_LAMBDA3 = SW_CL_EPS * 23.778;
const double SW_CL_GAM1SIG = SW_CL_CL_CL_GAMMA1 * SW_CL_SIGMA;
const double SW_CL_GAM2SIG = SW_CL_CL_CL_GAMMA2 * SW_CL_SIGMA;

//************************************************************************
double config::ClCl(particle* I, particle* J){
  rija=1/(rij-SW_CL_CL_a1);
  if (rija < 0){
    expa1ijk=SW_CL_CL_A*SpExp(SW_CL_GAM1SIG*rija);
    Fij=-expa1ijk*
      (SW_CL_CL_SIGMA5/pow(rij,5)*(5/rij+SW_CL_GAM1SIG*rija*rija)-
       SW_CL_CL_B*SW_CL_CL_SIGMA6/pow(rij,6)*
       (6/rij+SW_CL_GAM1SIG*rija*rija))/rij*Rij;
    
    I->F+=Fij; J->F-=Fij;
    return expa1ijk*
      (SW_CL_CL_B*SW_CL_CL_SIGMA6/pow(rij,6)-SW_CL_CL_SIGMA5/pow(rij,5));
  }
  else return 0;
}
//************************************************************************
double config::ClClCl(particle* I, particle* J, particle* K){
  rija=1/(rij-SW_CL_CL_a1); 
  rika=1/(rik-SW_CL_CL_a1); 
  rjka=1/(rjk-SW_CL_CL_a1);
  rija2=1/(rij-SW_CL_CL_a2); 
  rika2=1/(rik-SW_CL_CL_a2); 
  rjka2=1/(rjk-SW_CL_CL_a2);
  
  cosOjik2=SW_CL_CL_CL_LAMBDA2-cosOjik*cosOjik*SW_CL_CL_CL_LAMBDA3;
  cosOijk2=SW_CL_CL_CL_LAMBDA2-cosOijk*cosOijk*SW_CL_CL_CL_LAMBDA3;
  cosOikj2=SW_CL_CL_CL_LAMBDA2-cosOikj*cosOikj*SW_CL_CL_CL_LAMBDA3;

  if (rija<0 && rika<0){ // jik
    expa1jik=SpExp(SW_CL_GAM1SIG*(rija+rika));
    expa2jik=(rija2<0 && rika2<0)?
      SpExp(SW_CL_GAM2SIG*(rija2+rika2)):0;
    hjik2=SW_CL_CL_CL_LAMBDA1*pow(SW_CL_CL_SIGMA2/rij/rik,2.056)*expa1jik;
    hjik=hjik2+cosOjik2*expa2jik;
  }
  else{expa1jik=expa2jik=hjik=hjik2=0;}
  
  if (rija<0 && rjka<0){ // ijk
    expa1ijk=SpExp(SW_CL_GAM1SIG*(rija+rjka));
    expa2ijk=(rija2<0 && rjka2<0)?
      SpExp(SW_CL_GAM2SIG*(rija2+rjka2)):0;
    hijk2=SW_CL_CL_CL_LAMBDA1*pow(SW_CL_CL_SIGMA2/rij/rjk,2.056)*expa1ijk;
    hijk=hijk2+cosOijk2*expa2ijk;
  }
  else{expa1ijk=expa2ijk=hijk=hijk2=0;}
  
  if (rika<0 && rjka<0){ // ikj
    expa1ikj=SpExp(SW_CL_GAM1SIG*(rika+rjka));
    expa2ikj=(rika2<0 && rjka2<0)?
      SpExp(SW_CL_GAM2SIG*(rika2+rjka2)):0;
    hikj2=SW_CL_CL_CL_LAMBDA1*pow(SW_CL_CL_SIGMA2/rik/rjk,2.056)*expa1ikj;
    hikj=hikj2+cosOikj2*expa2ikj;
  }
  else{expa1ikj=expa2ikj=hikj=hikj2=0;}
  
  rija*=rija*SW_CL_GAM1SIG;
  rika*=rika*SW_CL_GAM1SIG;
  rjka*=rjka*SW_CL_GAM1SIG;
  rija2*=rija2*SW_CL_GAM2SIG;
  rika2*=rika2*SW_CL_GAM2SIG;
  rjka2*=rjka2*SW_CL_GAM2SIG;

  I->F+=
    hjik2*((2.056/rij + rija)*Ruij+(2.056/rik + rika)*Ruik)
    +expa2jik*(2*SW_CL_CL_CL_LAMBDA3*cosOjik*
	       ((1/rik-cosOjik/rij)*Ruij+(1/rij-cosOjik/rik)*Ruik)
	       +cosOjik2*(Ruij*rija2+Ruik*rika2))
    +hijk2*(2.056/rij+rija)*Ruij
    +expa2ijk*(2*SW_CL_CL_CL_LAMBDA3*cosOijk*
	       ((1/rjk-cosOijk/rij)*Ruij-(rik/rjk/rij)*Ruik)
	       +(cosOijk2*rija2)*Ruij)
    +hikj2*(2.056/rik+rika)*Ruik
    +expa2ikj*(2*SW_CL_CL_CL_LAMBDA3*cosOikj*
	       ((1/rjk-cosOikj/rik)*Ruik-(rij/rjk/rik)*Ruij)
	       +(cosOikj2*rika2)*Ruik);
  
  J->F-=
    hjik2*(2.056/rij+rija)*Ruij
    -expa2jik*(2*SW_CL_CL_CL_LAMBDA3*cosOjik*
	       ((cosOjik/rij-1/rik)*Ruij-(rjk/rij/rik)*Rujk)
	       -(cosOjik2*rija2)*Ruij)
    +hijk2*((2.056/rij+rija)*Ruij-(2.056/rjk+rjka)*Rujk)
    -expa2ijk*(2*SW_CL_CL_CL_LAMBDA3*cosOijk*
	       ((cosOijk/rij-1/rjk)*Ruij+(1/rij-cosOijk/rjk)*Rujk)
	       -cosOijk2*(Ruij*rija2-Rujk*rjka2))
    -hikj2*(2.056/rjk+rjka)*Rujk
    -expa2ikj*(2*SW_CL_CL_CL_LAMBDA3*cosOikj*
	       ((1/rik-cosOikj/rjk)*Rujk+(rij/rjk/rik)*Ruij)
	       +(cosOikj2*rjka2)*Rujk);
  
  K->F-=
    hjik2*(2.056/rik+rika)*Ruik
    -expa2jik*(2*SW_CL_CL_CL_LAMBDA3*cosOjik*
	       ((cosOjik/rik-1/rij)*Ruik+(rjk/rij/rik)*Rujk)
	       -(cosOjik2*rika2)*Ruik)
    +hijk2*(2.056/rjk+rjka)*Rujk
    -expa2ijk*(2*SW_CL_CL_CL_LAMBDA3*cosOijk*
	       ((cosOijk/rjk-1/rij)*Rujk+(rik/rij/rjk)*Ruik)
	       -(cosOijk2*rjka2)*Rujk)
    +hikj2*((2.056/rik+rika)*Ruik+(2.056/rjk+rjka)*Rujk)
    -expa2ikj*(2*SW_CL_CL_CL_LAMBDA3*cosOikj*
	       ((cosOikj/rik-1/rjk)*Ruik+(cosOikj/rjk-1/rik)*Rujk)
	       -cosOikj2*(Ruik*rika2+Rujk*rjka2));
  
  return hjik+hijk+hikj;
}
