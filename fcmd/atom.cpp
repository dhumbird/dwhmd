//atom.cpp

#include "atom.h"

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
  is_fixed=0; Nmap[6]=Nmap[9]=Nmap[14]=Nmap[17]=0;
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
//*************************************************************
void atom::EraseNbr(atom* j){
  VNI s=nlist.begin(); 
  for (s=nlist.begin(); s!=nlist.end(); s++)
    if ((atom*)s->a2==j) break;
  if (s!=nlist.end())
    nlist.erase(s);
}
//**************************************************************
void atom::AddNbrs(){
  Nt=0;
  for (map<short,double>::iterator m=Nmap.begin(); m!=Nmap.end(); m++){
    Nt+=m->second;
    //cerr<<ix<<" "<<Nt<<endl;
  }
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






