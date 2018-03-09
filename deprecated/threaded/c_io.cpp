#include "config.h"

//**********************************************************************
void config::SetProps(){
  if (!chem_is_def){
    ChemDef();
    chem_is_def=1;
  }
  for (i=begin; i<end; i++) i->SetProps();
}
//*******************************************************************
void config::Load(string cfgfile){
  ifstream fin(cfgfile.c_str(), ios::in);
  if(!fin){
    cerr<<"Can't open specified cfg file: "<<cfgfile<<endl;
    exit(1);
  }
  fin>>N>>NFixed>>mat>>T>>Lx>>Ly>>Lz>>t;
  resize(N);

  for (i=begin; i<end; i++) fin>>*i;
   
  fin.close();
  SetProps();
  name=cfgfile.replace(cfgfile.find("_"),11,"");
  Lz2=fabs(Lz/2);
}
//************************************************************************
void config::Dump(bool ascii){
  
  string out_cfg=name;
  out_cfg+="_";
  out_cfg+=time2string(t);
  if(!ascii){
    out_cfg+=".cfg";
    fstream fout(out_cfg.c_str(), ios::out);
    fout.fill(' ');
    fout.precision(12);
    fout.setf(ios::fixed);
    fout<<N<<" "<<NFixed<<" "<<mat<<" "<<T<<" "<<Lx<<" "<<Ly<<" "<<Lz
	<<" "<<t<<endl;
    for (i=begin; i<end; i++){
      fout<<*i<<endl;
    }
  }
  else{
    out_cfg+=".txt";
    fstream fout(out_cfg.c_str(), ios::out);
    fout<<N<<" "<<NFixed<<" "<<mat<<" "<<T<<" "<<Lx<<" "<<Ly<<" "<<Lz
	<<" "<<t<<endl;
    for (i=begin; i<end; i++){
      fout<<i->ix<<" "<<i->id<<" "<<i->R.x<<" "<<i->R.y<<" "<<i->R.z<<" "
	  <<i->V.x<<" "<<i->V.y<<" "<<i->V.z<<" "<<i->is_fixed<<endl;
    }
  }
}



