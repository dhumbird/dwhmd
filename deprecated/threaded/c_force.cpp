#include "config.h"

svector Rij, Rik, Rjk, Fij;	
svector Ruij, Ruik, Rujk;	
double rij, rik, rjk, cosOjik, cosOijk, cosOikj;
double rija, rjka, rika, d3;
double hjik, hjik2, hijk, hijk2, hikj, hikj2;

//**********************************************************************
void config::ForceEval(bool distrib){
  //evaluate forces, increment potential energy
  for (i=begin; i<end; i++)
    for(j0=i->nlist.begin(); j0!=i->nlist.end(); j0++)
      if(i<(*j0)){ //because the list contains all neighbors
	if (distrib){
	  d3=TwoBody(i, *j0);
	  u+=d3;
	  i->u+=d3;
	  (*j0)->u+=d3;
	}
	else u+=TwoBody(i, *j0);
	//only do what follows if a three-body pot. exists
	if (i->three_body){
	  trips_list.clear();

	  //trips_list contain a sorted, unique merge of i & j's neighbors
	  //add only j-neighbors that are > i
	  for (trip=(*j0)->nlist.begin(); trip!=(*j0)->nlist.end(); trip++)
	    if(*trip>i) trips_list.insert(*trip);

	  //add the rest of i's list to the merged list
	  trip=j0;
	  for (trip++; trip!=i->nlist.end(); trip++)
	    trips_list.insert(*trip);
		  
	  //this avoids repeating triplets by deleting anything that has
	  //already been tested as part of i's list.
	  for(trip=i->nlist.begin(); trip!=j0; trip++)
	    trips_list.erase(*trip);
	  
	  for(k0=trips_list.begin(); k0!=trips_list.end(); k0++)
	    if (distrib){
	      d3=ThreeBody(i, *j0, *k0);
	      u+=d3;
	      i->u+=d3;
	      (*j0)->u+=d3;
	      (*k0)->u+=d3;
	    }
	    else u+=ThreeBody(i, *j0, *k0);
	}
      }
}
//************************************************************************
double config::TwoBody(particle* I, particle* J){
  if((I->id + J->id)==28){  //Stillinger-Weber Si/Si Calculation
    Rij=(I->R-J->R); Rij.minimg(Lx, Ly, Lz);
    rij=Rij.mag();
    rija=1/(rij-a_SW);

    Fij=-EPSA_SW*SpExp(SIGMA_SW*rija)
      *(SIGMA_SW*rija*rija*(1-B_SIGMA4_SW/pow(rij,4))
	-4*B_SIGMA4_SW/pow(rij,5))/rij*Rij;
    
    I->F += Fij;
    J->F -= Fij;
        
    return EPSA_SW*SpExp(SIGMA_SW*rija)*(B_SIGMA4_SW/pow(rij,4)-1);
  }
  
  else if ((I->id + J->id)==32){ //Moliere-type Ar/Si calculation
    Rij=(I->R-J->R); Rij.minimg(Lx, Ly, Lz);
    rij=Rij.mag();
    rija=rij/a_MO_AR_SI;
    
    Fij=EPS_MO_AR_SI/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/a_MO_AR_SI)
			  +0.55*SpExp(-1.2*rija)*(1/rij+1.2/a_MO_AR_SI)
			  +0.1*SpExp(-6*rija)*(1/rij+6/a_MO_AR_SI))/rij*Rij;
    
    I->F += Fij;
    J->F -= Fij;
    
    return EPS_MO_AR_SI/rij*(0.35*SpExp(-0.3*rija)
			     +0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
  }
  else return 0;
}
//**********************************************************************
double config::ThreeBody (particle* I, particle* J, particle* K){
  if((I->id + J->id + K->id)==42){ //Still-Web Si/Si/Si calculation
    //Rij, rij have already been computed.
    Rik=(I->R - K->R); Rik.minimg(Lx, Ly, Lz);
    Rjk=(J->R - K->R); Rjk.minimg(Lx, Ly, Lz);
    rik=Rik.mag();
    rjk=Rjk.mag();
    Ruij=Rij*(1/rij);
    Ruik=Rik*(1/rik);
    Rujk=Rjk*(1/rjk);
    rija=1/(rij-a_SW); 
    rika=1/(rik-a_SW); 
    rjka=1/(rjk-a_SW);

    //remember that ^ means dot product
    cosOijk=(Rij^Rjk)*(-1/rij/rjk);
    cosOjik=(Rij^Rik)*(1/rij/rik);
    cosOikj=(Rik^Rjk)*(1/rik/rjk);

    if(rija<0 && rika<0){
      hjik2=LAMBDA_SW*SpExp(GAMMA_SW*(rija+rika))*(cosOjik+ONE_THIRD);
      hjik=hjik2*(cosOjik+ONE_THIRD);
    }
    else{hjik=0; hjik2=0;}
    
    if(rija<0 && rjka<0){
      hijk2=LAMBDA_SW*SpExp(GAMMA_SW*(rija+rjka))*(cosOijk+ONE_THIRD);
      hijk=hijk2*(cosOijk+ONE_THIRD);
    }
    else{hijk=0; hijk2=0;}
    
    if(rika<0 && rjka<0){
      hikj2=LAMBDA_SW*SpExp(GAMMA_SW*(rika+rjka))*(cosOikj+ONE_THIRD);
      hikj=hikj2*(cosOikj+ONE_THIRD);
    }
    else{hikj=0; hikj2=0;}

    rika=rika*rika*GEPS_SW*(hjik+hikj);
    rija=rija*rija*GEPS_SW*(hijk+hjik);
    rjka=rjka*rjka*GEPS_SW*(hijk+hikj);

    I->F+=Ruik*(rika-EPS2_SW*(hjik2*(1/rij-cosOjik/rik)
			      +hikj2*(1/rjk-cosOikj/rik)
			      -hijk2*rik/rjk/rij))
      +Ruij*(rija-EPS2_SW*(hjik2*(1/rik-cosOjik/rij)
			   +hijk2*(1/rjk-cosOijk/rij)
			   -hikj2*rij/rjk/rik));
    
    J->F+=Rujk*(rjka-EPS2_SW*(hikj2*(1/rik-cosOikj/rjk)
			      +hijk2*(1/rij-cosOijk/rjk)
			      -hjik2*rjk/rik/rij))
      -Ruij*(rija+EPS2_SW*(hjik2*(cosOjik/rij-1/rik)
			   +hijk2*(cosOijk/rij-1/rjk)
			   +hikj2*rij/rjk/rik));
    
    K->F-=Ruik*(rika+EPS2_SW*(hjik2*(cosOjik/rik-1/rij)
			      +hijk2*rik/rij/rjk
			      +hikj2*(cosOikj/rik-1/rjk)))
      +Rujk*(rjka+EPS2_SW*(hjik2*rjk/rik/rij
			   +hijk2*(cosOijk/rjk-1/rij)
			   +hikj2*(cosOikj/rjk-1/rik)));
    
    return EPS_SW*(hjik+hijk+hikj);
  }
  else return 0;
}
//**************************************************************************

