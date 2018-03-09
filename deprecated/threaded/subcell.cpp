#include "subcell.h"

//  //************************************************************************
//  void *subcell::th_s_ForceEval(void * sc){
//    set<subcell*> *This=(set<subcell*> *)sc;
//    for (set<subcell*>::iterator ssi=This->begin(); ssi!=This->end(); ssi++)
//      (*ssi)->s_ForceEval();
//    return sc;
//  }
//  //**************************************************************************
//  void subcell::s_ForceEval(){
//    //evaluate forces, increment potential energy
//    for (i=abegin(); i!=aend(); i++)
//      for(j=(*i)->nlist.begin(); j!=(*i)->nlist.end(); j++)
//        if(*i < *j){ //because the list contains all neighbors
//  	u+=TwoBody(*i, *j);
//   	//only do what follows if a three-body pot. exists
//  	if ((*i)->three_body){
//  	  trips_list.clear();
//    	  //trips_list contain a sorted, unique merge of i & j's neighbors
//    	  //add only j-neighbors that are > i
//  	  for (trip=(*j)->nlist.begin(); trip!=(*j)->nlist.end(); trip++)
//  	    if(*trip>*i) trips_list.insert(*trip);
//    	  //add the rest of i's list to the merged list
//    	  trip=j;
//    	  for (trip++; trip!=(*i)->nlist.end(); trip++)
//    	    trips_list.insert(*trip);
//    	  //this avoids repeating triplets by deleting anything that has
//    	  //already been tested as part of i's list.
//    	  for(trip=(*i)->nlist.begin(); trip!=j; trip++)
//    	    trips_list.erase(*trip);
//    	  for(k=trips_list.begin(); k!=trips_list.end(); k++)
//    	    u+=ThreeBody(*i, *j, *k);
//  	}
//        }
//  }
//  //*************************************************************************
//  double subcell::TwoBody(particle* I, particle* J){
//    if((I->id + J->id)==28){ //Stillinger-Weber Si/Si Calculation
//      Rij=(I->R-J->R); Rij.minimg(Lx, Ly, Lz);
//      rij=Rij.mag();
    
//      Fij=-1*EPS_SW*A_SW*SpExp(SIGMA_SW/(rij-a_SW))
//        *(SIGMA_SW/pow(rij-a_SW,2)*(1-B_SIGMA4_SW/pow(rij,4))
//  	-4*B_SIGMA4_SW/pow(rij,5))/rij*Rij;
    
//      I->F += Fij;
//      J->F -= Fij;
        
//      return EPS_SW*(A_SW*SpExp(SIGMA_SW/(rij-a_SW))*(B_SIGMA4_SW/pow(rij,4)-1));
//    }
  
//    else if ((I->id + J->id)==32){ //Moliere-type Ar/Si calculation
//      Rij=(I->R-J->R); Rij.minimg(Lx, Ly, Lz);
//      rij=Rij.mag();
//      rija=rij/a_MO_AR_SI;
    
//      Fij=EPS_MO_AR_SI/rij*(0.35*SpExp(-0.3*rija)*(1/rij+0.3/a_MO_AR_SI)
//  			  +0.55*SpExp(-1.2*rija)*(1/rij+1.2/a_MO_AR_SI)
//  			  +0.1*SpExp(-6*rija)*(1/rij+6/a_MO_AR_SI))/rij*Rij;
    
//      I->F += Fij;
//      J->F -= Fij;
    
//      return EPS_MO_AR_SI/rij*(0.35*SpExp(-0.3*rija)
//  			     +0.55*SpExp(-1.2*rija)+0.1*SpExp(-6*rija));
//    }
//    else return 0;
//  }
//  //****************************************************************************
//  double subcell::ThreeBody (particle* I, particle* J, particle* K){
  
//    if((I->id + J->id + K->id) == 42){ //Still-Web Si/Si/Si calculation
//      //Rij, rij have already been computed.
//      Rik=(I->R - K->R); Rik.minimg(Lx, Ly, Lz);
//      Rjk=(J->R - K->R); Rjk.minimg(Lx, Ly, Lz);
//      rik=Rik.mag();
//      rjk=Rjk.mag();
//      Ruij=Rij*(1/rij);
//      Ruik=Rik*(1/rik);
//      Rujk=Rjk*(1/rjk);
//      rija=1/(rij-a_SW); 
//      rika=1/(rik-a_SW); 
//      rjka=1/(rjk-a_SW);

//      //remember that ^ means dot product
//      cosOijk=(Rij^Rjk)*(-1/rij/rjk);
//      cosOjik=(Rij^Rik)*(1/rij/rik);
//      cosOikj=(Rik^Rjk)*(1/rik/rjk);

//      if(rija<0 && rika<0){
//        hjik2=LAMBDA_SW*SpExp(GAMMA_SW*(rija+rika))*(cosOjik+ONE_THIRD);
//        hjik=hjik2*(cosOjik+ONE_THIRD);
//      }
//      else{hjik=0; hjik2=0;}
    
//      if(rija<0 && rjka<0){
//        hijk2=LAMBDA_SW*SpExp(GAMMA_SW*(rija+rjka))*(cosOijk+ONE_THIRD);
//        hijk=hijk2*(cosOijk+ONE_THIRD);
//      }
//      else{hijk=0; hijk2=0;}
    
//      if(rika<0 && rjka<0){
//        hikj2=LAMBDA_SW*SpExp(GAMMA_SW*(rika+rjka))*(cosOikj+ONE_THIRD);
//        hikj=hikj2*(cosOikj+ONE_THIRD);
//      }
//      else{hikj=0; hikj2=0;}

//      rika=rika*rika*geps*(hjik+hikj);
//      rija=rija*rija*geps*(hijk+hjik);
//      rjka=rjka*rjka*geps*(hijk+hikj);

//      I->F+=Ruik*(rika-eps2*(hjik2*(1/rij-cosOjik/rik)+hikj2*(1/rjk-cosOikj/rik)
//      		   -hijk2*rik/rjk/rij))
//      +Ruij*(rija-eps2*(hjik2*(1/rik-cosOjik/rij)+hijk2*(1/rjk-cosOijk/rij)
//      		-hikj2*rij/rjk/rik));
    
//      J->F+=Rujk*(rjka-eps2*(hikj2*(1/rik-cosOikj/rjk)+hijk2*(1/rij-cosOijk/rjk)
//      			   -hjik2*rjk/rik/rij))
//       -Ruij*(rija+eps2*(hjik2*(cosOjik/rij-1/rik)+hijk2*(cosOijk/rij-1/rjk)
//      		+hikj2*rij/rjk/rik));
    
//      K->F-=Ruik*(rika+eps2*(hjik2*(cosOjik/rik-1/rij)+hijk2*rik/rij/rjk
//  			   +hikj2*(cosOikj/rik-1/rjk)))
//        +Rujk*(rjka+eps2*(hjik2*rjk/rik/rij+hijk2*(cosOijk/rjk-1/rij)
//  			+hikj2*(cosOikj/rjk-1/rik)));
    
//      return EPS_SW*(hjik+hijk+hikj);
//    }
//    else return 0;
//  }
//**************************************************************************
void *subcell::th_s_ReNeighbor(void * sc){
  set<subcell*> *This=(set<subcell*> *)sc;
  for (set<subcell*>::iterator ssi=This->begin(); ssi!=This->end(); ssi++)
    (*ssi)->s_ReNeighbor();
  return sc;
}
//**************************************************************************
void subcell::s_ReNeighbor(){
  for (i=abegin(); i!=aend(); i++)  
      for (c=nbegin(); c!=nend(); c++)
	for (j=(*c)->abegin(); j!=(*c)->aend(); j++)
	  if (*i < *j)
	    if (((*i)->R - (*j)->R).minsqmag(Lx, Ly, Lz) < 
		((*i)->rc >? (*j)->rc)){
	      (*i)->nlist.insert(*j);
	      (*j)->nlist.insert(*i);
	    }
}
