#include "rcfg.h"

class ix_less : public binary_function<atom,atom,bool>{
public:
  bool operator()(atom x, atom y){return x.ix < y.ix;}
};

void rcfg::Load(string cfgfile, float radscale){
  if (FileExists(cfgfile)){
    int critical=0;
    int commands=0;
    ifstream fin(cfgfile.c_str(), ios::in);
    string s="#";
    string arg;
    while (!fin.eof()){
      getline(fin,s);
      if (sfind(s,"#")){
        istringstream ist(s);
        ist>>arg;
        if (arg=="#N") ist>>N;
	else if (sfind(arg, "ax")) ist>>Nmax;
	else if (sfind(arg, "imensions")) ist>>Lx>>Ly>>Lz;
	else if (sfind(arg, "#-")) break;
	else if (sfind(arg, "ime"))ist>>t;
      }
    }
    resize(N);
    Lz=fabs(Lz);
    the_top=-1000;
    the_bottom=1000;
    for (i=begin; i<end; i++){ 
      fin>>*i;
      i->SetProps(radscale);
      if (i->R.z > the_top) the_top=i->R.z;
      if (i->R.z < the_bottom) the_bottom=i->R.z;
    }
    fin.close();
    name=cfgfile.replace(cfgfile.find("_"),100,"");
  }
  else{
    cerr<<"Can't open specified cfg file: "<<cfgfile<<endl;
    exit(1);
  } 
}
//*********************************************************************
void rcfg::init(){
  atom_list=new vector<atom>;
  name="temp";
  N=0; 
  t=Lx=Ly=Lz=0; 
}
//*****************************************************************************
VAI rcfg::atomix(int index){
  VAI p=begin;
  while (p<end && p->ix!=index) p++;
  return p;
}
//**************************************************************************
void rcfg::resize(int s){
  N=s;
  atom_list->resize(s);
  begin=atom_list->begin();
  end=atom_list->end();
}
//**********************************************************************
void rcfg::ix_sort_asc(){
  sort(atom_list->begin(), atom_list->end(), ix_less());
  begin=atom_list->begin(); end=atom_list->end();
}

void rcfg::ix_sort_desc(){
  ix_sort_asc();
  reverse(atom_list->begin(), atom_list->end());
  begin=atom_list->begin(); end=atom_list->end();
}
//**********************************************************************
void rcfg::ReNeighbor(){
  double cut=0;
  for (i=begin; i<end; i++) i->nlist.clear();
  for (i=begin; i<end; i++)
    for (j=begin; j<end; j++)  
      if (i < j){
	switch (i->id + j->id){
	case 12: 
	  //cut=2.092; break;
	  cut=3; break;
	case 7: 
	  cut=2.25; break;
	case 2: 
	  cut=1.21; break;
	case 20: 
	  cut=4.141225; break;
	case 28: 
	  cut=6.682; break;
	  //cut=8; break;
	case 23: 
	  cut=3.1006; break;
	case 15: 
	  cut=1.941; break;
	}
	if ((i->R - j->R).minsqmag(Lx, Ly) < cut){
	  i->nlist.insert(&(*j));
	  j->nlist.insert(&(*i));
	}
      }
}
// //**********************************************************************
// double rcfg::BondOrder(nbr* bond_ij, double Pre){
//   svector Rik, Rjk, R;
//   double rik, rjk;
//   iNconj=jNconj=0;
  
//   atom* i=(atom*)bond_ij->a1;
//   atom* j=(atom*)bond_ij->a2;
 
//   double NF_ij, NC_ij, NSi_ij, Nt_ij, NF_ji, NC_ji, NSi_ji, Nt_ji;
//   double NCl_ij, NCl_ji;
//   double P_ij, P_ji, dFP_ij, dFP_ji, dCP_ij, dCP_ji, dlam;
//   double cosO, el, elf, g, g1, xik, Yik, Y1ik, F_ij, dFi_ij, dFj_ij, dFc_ij;
//   double eta_i, eta_j, delta_i, delta_j;

//   NF_ij  = i->Nmap[9]  - (j->id==9)  * bond_ij->f;
//   NC_ij  = i->Nmap[6]  - (j->id==6)  * bond_ij->f;
//   NSi_ij = i->Nmap[14] - (j->id==14) * bond_ij->f;
//   NCl_ij = i->Nmap[17] - (j->id==17) * bond_ij->f;
//   NF_ji  = j->Nmap[9]  - (i->id==9)  * bond_ij->f;
//   NC_ji  = j->Nmap[6]  - (i->id==6)  * bond_ij->f;
//   NSi_ji = j->Nmap[14] - (i->id==14) * bond_ij->f;
//   NCl_ji = j->Nmap[17] - (i->id==17) * bond_ij->f;

//   Nt_ij = NF_ij + NC_ij + NSi_ij + NCl_ij;
//   Nt_ji = NF_ji + NC_ji + NSi_ji + NCl_ji;

//   P_ij=P_ji=dFP_ij=dFP_ji=dCP_ij=dCP_ji=0;
  
//   if (bond_ij->type==12){
//     Pcc_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
//     Pcc_bicubicint(NF_ji, NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
//   }
//   else if (bond_ij->type==15){
//     if (i->id==6) Pcf_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
//     else Pcf_bicubicint(NF_ji, NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
//   }
//   else if (bond_ij->type==23){
//     if (i->id==14) Psif_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
//     else Psif_bicubicint(NF_ji,  NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
//   }
//   else if (bond_ij->type==31){
//     if (i->id==14) Psicl_bicubicint(NCl_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
//     else Psicl_bicubicint(NCl_ji,  NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
//   }
//   bsp_ij=P_ij; bsp_ji=P_ji;

//   //***************for k neighbor of i*************************
//   for (VNI bond_ik=i->nlist.begin(); bond_ik!=i->nlist.end(); bond_ik++){
//     atom* k=(atom*)bond_ik->a2;
//     if (k!=j){
//       VNI bond_jk=j->FindNbr(k);
//       if (bond_jk!=j->nlist.end()){
// 	rjk=bond_jk->r;
// 	Rjk=bond_jk->Rhat;
//       }
//       else{
// 	Rjk=j->R-k->R;
// 	Rjk.minimg(Lx,Ly);
// 	rjk=Rjk.mag();
// 	Rjk/=rjk;
//       }
//       //*******b-sigma-pi calculations************
//       elambda(i->id, j->id, k->id, bond_ij->r, bond_ik->r, &el, &dlam);
//       elf=el*bond_ik->f;
      
//       cosO = bond_ij->Rhat ^ bond_ik->Rhat;
//       G_angle(i->id, j->id, k->id, cosO, &g, &g1);
//       bsp_ij += g*elf;
//       //*********conjugation*****************
//       if (bond_ij->type==12 && k->id==6){
// 	PolySwitch(k->Nt - bond_ik->f - 2, &Yik, &Y1ik);
// 	iNconj += bond_ik->f * Yik;
//       }
//     }
//   }

//   //***************for k neighbor of j*************************
//   for (VNI bond_jk=j->nlist.begin(); bond_jk!=j->nlist.end(); bond_jk++){
//     atom* k=(atom*)bond_jk->a2;
//     if (k!=i){
//       VNI bond_ik=i->FindNbr(k);
//       if (bond_ik!=i->nlist.end()){
// 	rik=bond_ik->r;
// 	Rik=bond_ik->Rhat;
//       }
//       else{
// 	Rik=i->R-k->R;
// 	Rik.minimg(Lx,Ly);
// 	rik=Rik.mag();
// 	Rik/=rik;
//       }
//       //*********b-sigma-pi calculation*************
//       elambda(j->id, i->id, k->id, bond_ij->r, bond_jk->r, &el, &dlam);
//       elf=el*bond_jk->f;
      
//       cosO = -(bond_ij->Rhat ^ bond_jk->Rhat);
//       G_angle(j->id, i->id, k->id, cosO, &g, &g1);
//       bsp_ji += g * elf;
//       //*********conjugation************
//       if ( bond_ij->type==12 && k->id==6){
// 	PolySwitch(k->Nt - bond_jk->f - 2, &Yik, &Y1ik);
// 	jNconj += bond_jk->f * Yik;
//       }
//     }
//   }
  
//   //***************conjugation************************
//   if (bond_ij->type==12){
//     Nconj_ij=1+iNconj+jNconj;
//     Fcc_tricubicint(Nt_ij, Nt_ji, Nconj_ij, &F_ij, &dFi_ij, &dFj_ij, &dFc_ij);
//   }
//   else 
//     F_ij=dFi_ij=dFj_ij=dFc_ij=0;
//   //***********************************************
//   double dtemp=bsp_ji;
//   set_eta(i->id, j->id, &eta_i, &delta_i);
//   set_eta(j->id, i->id, &eta_j, &delta_j);
//   bbar_ij = 0.5*(pow(1 + pow(bsp_ij, eta_i), -delta_i)
// 		 + pow(1 + pow(bsp_ji, eta_j), -delta_j)
// 		 +F_ij);
//   return bbar_ij;
// }
