//atom.cpp

#include "atom.h"

extern map<short, double> MASS;
// class nbr_less : public binary_function<nbr,nbr,bool>{
// public:
//   bool operator()(nbr n1, nbr n2){return n1.a2 < n2.a2;}
// };

//***constructor**********************************************
atom::atom(){
  F.clear();
  nlist.clear();
  m=0;
  R.clear();
  V.clear();
  P.clear();
  id=0; ix=0;
  is_fixed=0;
  my_cell=p_clust=NULL;
}
//***************************************************************
// void atom::Sort_nlist(){
//   sort(nlist.begin(), nlist.end(), nbr_less());
//   for (stop=nlist.begin(); stop!=nlist.end(); stop++)
//     if (stop->a2 > this) break;
// }
//***************************************************************
VNI atom::FindNbr(atom* j){
  VNI s=nlist.begin(); 
  while(s!=nlist.end() && (atom*)s->a2 != j) s++;
  return s;
}
VNJI atom::LJ_FindNbr(atom* j){
  VNJI s=ljnlist.begin(); 
  while(s!=ljnlist.end() && (atom*)s->a2 != j) s++;
  return s;
}
//*************************************************************
void atom::EraseNbr(atom* j){
  VNI s=nlist.begin(); 
  for (s=nlist.begin(); s!=nlist.end(); s++)
    if ((atom*)s->a2==j) break;
  if (s!=nlist.end())
    nlist.erase(s);
}
//*************************************************************
void atom::LJ_EraseNbr(atom* j){
  VNJI s=ljnlist.begin(); 
  for (s=ljnlist.begin(); s!=ljnlist.end(); s++)
    if ((atom*)s->a2==j) break;
  if (s!=ljnlist.end())
    ljnlist.erase(s);
}
//**************************************************************
void atom::PostComp(){
  Nt=NC=NH=0;
  for (VNI ij=nlist.begin(); ij!=nlist.end(); ij++){
    switch (((atom*)ij->a2)->id){
    case 1:
      NH+=ij->f; break;
    case 6:
      NC+=ij->f; break;
    }
  }
  Nt=NH+NC;
}
//**********************************************************
void atom::LJ_PreComp(){
  for (VNI ij=nlist.begin(); ij!=nlist.end(); ij++){
    atom* j=(atom*)ij->a2;
    VNJI l_ij=LJ_FindNbr(j);
    if (l_ij != ljnlist.end()){
      if (ij->f==1) LJ_EraseNbr(j);
      else if (ij->f > l_ij->C){
	l_ij->C = ij->f;
	l_ij->n1 = &(*ij);
	l_ij->n2=l_ij->n3=NULL;
      }
    }
    for (VNI jk=j->nlist.begin(); jk!=j->nlist.end(); jk++){
      atom* k=(atom*)jk->a2;
      VNJI l_ik=LJ_FindNbr(k);
      if (l_ik != ljnlist.end()){
	double ff=jk->f * ij->f;
	if (ff == 1) LJ_EraseNbr(k);
	else if (ff > l_ik->C){
	  l_ik->C = ff;
	  l_ik->n1 = &(*ij);
	  l_ik->n2 = &(*jk);
	  l_ik->n3=NULL;
	}
      }
      for (VNI kl=k->nlist.begin(); kl!=k->nlist.end(); kl++){
	atom* l=(atom*)kl->a2;
	VNJI l_il=LJ_FindNbr(l);
	if (l_il != ljnlist.end()){
	  double fff = ij->f * jk->f * kl->f;
	  if (fff == 1) LJ_EraseNbr(l);
	  else if (fff > l_il->C){
	    l_il->C = fff;
	    l_il->n1 = &(*ij);
	    l_il->n2 = &(*jk);
	    l_il->n3 = &(*kl);
	  }
	}
      }
    }
  }
}
//**************************************************************
void atom::SetProps(){
  m=MASS[id]/9648.531;
}
//******************************************************************
void atom::PosUpdate(double& dt, double& dtsq2, double& Lx, double& Ly){
  R+=dt*V+dtsq2*F/m;
  while (abs(R.x/Lx)>0.5)
    if (R.x < 0) R.x+=Lx;
    else R.x-=Lx;
      
  while (abs(R.y/Ly)>0.5)
    if (R.y < 0) R.y+=Ly;
    else R.y-=Ly;
}
//**************************************************************
ostream& operator << (ostream& out, atom& p){
  out<<p.ix<<" "<<p.id<<" "<<p.R<<" "<<p.V<<" "<<p.is_fixed;
  return out;
}
//**************************************************************
istream& operator >> (istream& in, atom& p){
  in>>p.ix>>p.id>>p.R>>p.V>>p.is_fixed;
  return in;
}






