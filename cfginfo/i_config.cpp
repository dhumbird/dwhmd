#include "i_config.h"
extern map<short,double> RC_SQ;

//**********************************************************************
void config::ReNeighbor(){
  int type; double r; svector R;
  for (VAI i=begin; i<end; i++) i->nlist.clear();
  for (VAI i=begin; i<end; i++){
    for (VAI j=begin; j<end; j++){
      if (i<j){
	R = i->R-j->R; R.minimg(Lx,Ly);
	if ((r=R.sqmag()) < RC_SQ[i->id+j->id]){
	  i->nlist.insert(&(*j));
	  j->nlist.insert(&(*i));
	}
      }
    }
  }
}
//**********************************************************************
void config::Load(string cfgfile){
  if (FileExists(cfgfile)){
    int critical=0;
    int commands=0;
    ifstream fin(cfgfile.c_str(), ios::in);
    string s="#";
    string arg;
    bool kill_me=0;
    while (!kill_me && !fin.eof()){
      getline(fin,s);
      if (sfind(s,"#")){
	istringstream ist(s);
	ist>>arg;
	if (arg=="#N"){ist>>N; critical++;}
	else if (sfind(arg, "ixed")){ist>>NFixed; critical++;}
	else if (sfind(arg, "max")){ist>>Nmax; critical++;}
	else if (sfind(arg, "emperature")){ist>>T; critical++;}
	else if (sfind(arg, "imensions")){ist>>Lx>>Ly>>Lz; critical++;}
	else if (sfind(arg, "ime")){ist>>t; critical++;}
	else if (sfind(arg, "#-")) kill_me=1;
	else{
          cerr<<"I don't understand the argument "<<arg<<" in "<<cfgfile<<endl;
          exit(1);
        }
      }
    }
    if (critical < 6){
      cerr<<"Hmm, something seems to be missing from "<<cfgfile<<endl;
      exit(1);
    }
    resize(N);
    for (VAI i=begin; i<end; i++) fin>>*i;
    fin.close();
    SetProps();
    name=cfgfile.replace(cfgfile.find("_"),14,"");
  }
  else{
    cerr<<"Can't open specified cfg file: "<<cfgfile<<endl;
    exit(1);
  } 
}
//*********************************************************************
void config::init(){
  ChemDef();
  atom_list=new vector<atom>;
  name="temp";
  t=0; //doubles
  N=NFixed=Nmax=0; //ints
  Lx=Ly=Lz=T=0; //doubles
}
//*****************************************************************************
VAI config::atomix(int index){
  VAI p=begin;
  while (p<end && p->ix!=index) p++;
  return p;
}
//**************************************************************************
void config::resize(int s){
  N=s;
  atom_list->resize(s);
  begin=atom_list->begin();
  end=atom_list->end();
}
