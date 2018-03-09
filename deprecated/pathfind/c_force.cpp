#include "config.h"

//**********************************************************************
void config::ForceEval(){
  //evaluate forces, increment potential energy
  for (i=begin; i<end; i++){
    for(j0=i->nlist.begin(); j0!=i->nlist.end(); j0++)
      if(i<(*j0)){ //because the list contains all neighbors
	u+=TwoBody(i, *j0);
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
	    u+=ThreeBody(i, *j0, *k0);
	}
      }
  }
}
//************************************************************************
double config::TwoBody(particle* I, particle* J){
  Rij=(I->R-J->R); Rij.minimg(Lx, Ly);
  rij=Rij.mag();
  
  switch(I->id + J->id){
  case 11: 
    return HeF(I,J); break;
  case 16: 
    return HeSi(I,J); break;
  case 18: 
    return FF(I,J); break;
  case 19: 
    if (I->id==10 || J->id==10) return NeF(I,J); 
    else return HeCl(I,J);
    break;
  case 23: 
    return FSi(I,J); break;
  case 24: 
    return NeSi(I,J); break;
  case 27:
    if (I->id==10 || J->id==10) return NeCl(I,J); 
    else return ArF(I,J); 
    break;
  case 28: 
    return SiSi(I,J); break;
  case 31: 
    return ClSi(I,J); break;
  case 32: 
    return ArSi(I,J); break;
  case 34: 
    return ClCl(I,J); break;
  case 35: 
    return ArCl(I,J); break;
  case 50:
    return KrSi(I,J); break;
  default:
    return 0; break;
  }
}
//**********************************************************************
double config::ThreeBody (particle* I, particle* J, particle* K){
  //Rij, rij have already been computed.
  //  Rij=(I->R - J->R); Rij.minimg(Lx,Ly);
  Rik=(I->R - K->R); Rik.minimg(Lx, Ly);
  Rjk=(J->R - K->R); Rjk.minimg(Lx, Ly);
  //rij=Rij.mag();
  rik=Rik.mag();
  rjk=Rjk.mag();
  Ruij=Rij*(1/rij);
  Ruik=Rik*(1/rik);
  Rujk=Rjk*(1/rjk);
  //remember that ^ means dot product
  cosOijk=(Rij^Rjk)*(-1/rij/rjk);
  cosOjik=(Rij^Rik)*(1/rij/rik);
  cosOikj=(Rik^Rjk)*(1/rik/rjk);
  
  switch(I->id + J->id + K->id){ 
  case 27:
      return FFF(I,J,K); break;
  case 32:
    return FFSi(I,J,K); break;
  case 37:
    return FSiSi(I,J,K); break;
  case 42:
    return SiSiSi(I,J,K); break;
  case 45:
    return ClSiSi(I,J,K); break;
  case 48:
    return ClClSi(I,J,K); break;
  case 51:
    return ClClCl(I,J,K); break;
  default:
    return 0; break;
  }
}
//************************************************************************
double config::u_TwoBody(particle* I, particle* J){
  Rij=(I->R-J->R); Rij.minimg(Lx, Ly);
  rij=Rij.mag();
  
  switch(I->id + J->id){
  case 18: 
    return u_FF(); break;
  case 23: 
    return u_FSi(); break;
  case 28: 
    return u_SiSi(); break;
  default:
    return 0; break;
  }
}
//**********************************************************************
double config::u_ThreeBody (particle* I, particle* J, particle* K){
  //Rij, rij have already been computed.
  //Rij=(I->R - J->R); Rij.minimg(Lx,Ly);
  Rik=(I->R - K->R); Rik.minimg(Lx, Ly);
  Rjk=(J->R - K->R); Rjk.minimg(Lx, Ly);
  //rij=Rij.mag();
  rik=Rik.mag();
  rjk=Rjk.mag();
  //remember that ^ means dot product
  cosOijk=(Rij^Rjk)*(-1/rij/rjk);
  cosOjik=(Rij^Rik)*(1/rij/rik);
  cosOikj=(Rik^Rjk)*(1/rik/rjk);
  
  switch(I->id + J->id + K->id){ 
  case 27:
      return u_FFF(); break;
  case 32:
    return u_FFSi(I,J,K); break;
  case 37:
    return u_FSiSi(I,J,K); break;
  case 42:
    return u_SiSiSi(); break;
  default:
    return 0; break;
  }
}
