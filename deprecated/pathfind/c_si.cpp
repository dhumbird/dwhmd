#include "config.h"
//**********************************************************************
double config::SiSi(particle* I, particle* J){
  rija=1/(rij-SW_SI_a);
  if (rija < 0){
    Fij=-SW_SI_SI_EPSA*SpExp(SW_SI_SIGMA*rija)
      *(SW_SI_SIGMA*rija*rija*(1-SW_SI_B_SIGMA4/pow(rij,4))
	-4*SW_SI_B_SIGMA4/pow(rij,5))/rij*Rij;
    I->F += Fij; J->F -= Fij;
    return SW_SI_SI_EPSA*SpExp(SW_SI_SIGMA*rija)*
      (SW_SI_B_SIGMA4/pow(rij,4)-1);
  }
  else return 0;
}
//*********************************************************************
double config::SiSiSi(particle* I, particle* J, particle* K){
  rija=1/(rij-SW_SI_a); 
  rika=1/(rik-SW_SI_a); 
  rjka=1/(rjk-SW_SI_a);
  
  if(rija<0 && rika<0){
    hjik2=SW_SI_SI_SI_LAMBDA*SpExp(SW_SI_SI_SI_GAMMA*(rija+rika))
      *(cosOjik+ONE_THIRD);
    hjik=hjik2*(cosOjik+ONE_THIRD);
  }
  else{hjik=hjik2=0;}
  
  if(rija<0 && rjka<0){
    hijk2=SW_SI_SI_SI_LAMBDA*SpExp(SW_SI_SI_SI_GAMMA*(rija+rjka))
      *(cosOijk+ONE_THIRD);
    hijk=hijk2*(cosOijk+ONE_THIRD);
  }
  else{hijk=hijk2=0;}
  
  if(rika<0 && rjka<0){
    hikj2=SW_SI_SI_SI_LAMBDA*SpExp(SW_SI_SI_SI_GAMMA*(rika+rjka))
      *(cosOikj+ONE_THIRD);
    hikj=hikj2*(cosOikj+ONE_THIRD);
  }
  else{hikj=hikj2=0;}
  
  rika=rika*rika*SW_SI_SI_SI_GEPS*(hjik+hikj);
  rija=rija*rija*SW_SI_SI_SI_GEPS*(hijk+hjik);
  rjka=rjka*rjka*SW_SI_SI_SI_GEPS*(hijk+hikj);
  
  I->F+=Ruik*(rika-SW_SI_SI_EPS2*(hjik2*(1/rij-cosOjik/rik)
				  +hikj2*(1/rjk-cosOikj/rik)
				  -hijk2*rik/rjk/rij))
    +Ruij*(rija-SW_SI_SI_EPS2*(hjik2*(1/rik-cosOjik/rij)
			       +hijk2*(1/rjk-cosOijk/rij)
			       -hikj2*rij/rjk/rik));
  
  J->F+=Rujk*(rjka-SW_SI_SI_EPS2*(hikj2*(1/rik-cosOikj/rjk)
				  +hijk2*(1/rij-cosOijk/rjk)
				  -hjik2*rjk/rik/rij))
    -Ruij*(rija+SW_SI_SI_EPS2*(hjik2*(cosOjik/rij-1/rik)
			       +hijk2*(cosOijk/rij-1/rjk)
			       +hikj2*rij/rjk/rik));
  
  K->F-=Ruik*(rika+SW_SI_SI_EPS2*(hjik2*(cosOjik/rik-1/rij)
				  +hijk2*rik/rij/rjk
				  +hikj2*(cosOikj/rik-1/rjk)))
    +Rujk*(rjka+SW_SI_SI_EPS2*(hjik2*rjk/rik/rij
			       +hijk2*(cosOijk/rjk-1/rij)
			       +hikj2*(cosOikj/rjk-1/rik)));
  
  return SW_SI_SI_EPS*(hjik+hijk+hikj);
}



