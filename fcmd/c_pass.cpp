#include "config.h"
#define VDW_C 1.7
#define VDW_H 1.2
#define C_C 1.54
#define C_H 1.1
extern map<short,double> RC_SQ;
//***********************************************************************
double config::U_on_i(atom* s){
  int type; svector R; double r;
  for (VNI J0=s->nlist.begin(); J0!=s->nlist.end(); J0++){
    ((atom*)J0->a2)->EraseNbr(s);
  }
  s->nlist.clear();
  subcell* subc = (subcell*)s->my_cell;
  for (c=subc->nbegin(); c!=subc->nend(); c++)
    for (sai_j=(*c)->abegin(); sai_j!=(*c)->aend(); sai_j++){
      atom* J=*sai_j;
      if (s!=J){
	type = s->id + J->id;
	R = s->R-J->R; R.minimg(Lx,Ly);
	if ((r=R.sqmag()) < RC_SQ[type]){
	  nbr n;
	  n.type=type;
	  r=sqrt(r);
	  n.a1=(void*)s;
	  n.a2=(void*)J;
	  n.r=r;
	  n.Rhat=R/r;
	  n.PreComp();
	  s->nlist.push_back(n);
	  n.Invert();
	  J->nlist.push_back(n);
	}
      }
    }
  s->AddNbrs();
  for (VNI J0=s->nlist.begin(); J0!=s->nlist.end(); J0++){
    ((atom*)J0->a2)->AddNbrs();
  }
  return Uonly()-u;
}
//**************************************************************************
void config::SurfAcc(short id){
  surf_list.clear();
//   u=Uonly();
//   float VDW=VDW_C;
//   if (id==1) VDW+=VDW_H;
//   else VDW+=VDW_C;
//   float VDW_SQ=pow(VDW,2);
//   svector P;
//   float bl=0;
//   if (id==1) bl=C_H;
//   else if (id==6) bl=C_C;
//   bool flag; int count;
//   float sin_phi_rad,phi_rad,increm, uu;
//   VAI newatom=append(id);
  for (VAI s=begin; s<end-1; s++){
    if (!s->is_fixed && s->id!=9){
      if (s->Nt < 3.3)
	
//       s->u=1000;
//       for (float phi=0; phi<=180; phi+=3){
// 	phi_rad=phi/180*PI;
// 	sin_phi_rad=sin(phi_rad);
// 	increm=PI/sin_phi_rad/90.0;
// 	for (float theta=0; phi!=0 && phi!=180 && theta<2*PI; theta+=increm){
// 	  P.set(VDW*sin_phi_rad*cos(theta), VDW*sin_phi_rad*sin(theta),
// 		VDW*cos(phi_rad));
// 	  P+=s->R;
// 	  P.minimg(Lx,Ly);
// 	  flag=1;
// 	  for (VAI J=begin; flag && J<end; J++)
// 	    if (J!=s && J->id==6)
// 	      if ((J->R - P).minsqmag(Lx,Ly) < VDW_SQ)
// 		flag=0;
// 	  if (flag){
	surf_list.insert(s->ix);
      // 	    P.set(bl*sin_phi_rad*cos(theta), bl*sin_phi_rad*sin(theta),
// 		  bl*cos(phi_rad));
// 	    P+=s->R;
// 	    P.minimg(Lx,Ly);
// 	    newatom->R=P;
// 	    uu=U_on_i(&(*newatom));
// 	    if (uu < s->u){
// 	      s->u=uu;
// 	      s->P=P;
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
//   erase(newatom);
    }
  }
}
//*************************************************************************
void config::Passivate(string id){
//   SurfAcc(id);
//   set<int>::iterator sii;
//   set<int> pass_list;
//   for (sii=surf_list.begin(); sii!=surf_list.end(); sii++){
//     VAI s=atomix(*sii);
//     if (s->id==6 && s->Nt<4){
//       pass_list.insert(*sii);
//     }
//   }
  
//   if (!pass_list.empty()){
//     int kept=pass_list.size();
//     cerr<<"* Passivating. Found "<<kept<<" sites...";
//     float bl=0;
//     svector P;
//     if (id==1) bl=C_H;
//     else if (id==6) bl=C_C;
//     else {cerr<<"\nI shouldn't passivate with this species!\n"; exit(1);}
    
//     for (sii=pass_list.begin(); sii!=pass_list.end(); sii++){
//       VAI newatom=append(id);
//       VAI target=atomix(*sii);
//       cerr<<"Target atom "<<*sii<<" at "<<target->P<<endl;
//       newatom->R=target->P;
//       RunQuenched(0.01,0,0,0,0,-1,1,1,0,300);
      
//       //if it is less than 110% of the bond length away, keep it.
//       if ((newatom->R - target->R).minmag(Lx,Ly)<bl*1.1){
//  	newatom->V.clear();
//       }
//       else{
//  	kept--;
//  	erase(newatom);
//  	Nmax--;
//       }
//     }
//     TimeStepInit(); ReNeighbor();
//     cerr<<" Kept "<<kept<<". Total atoms: "<<N<<endl;
//   }
}
