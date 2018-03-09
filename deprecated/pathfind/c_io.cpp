#include "config.h"

//*******************************************************************
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
	else if (sfind(arg, "mat")){ist>>mat; critical++;}
	else if (sfind(arg, "emperature")){ist>>T; critical++;}
	else if (sfind(arg, "imensions")){ist>>Lx>>Ly>>Lz; critical++;}
	else if (sfind(arg, "ime")){ist>>t; critical++;}
	else if (sfind(arg, "ask")){
	  b_masked=1;
	  my_mask=new mask;
	  ist>>*my_mask;
	}
	else if (sfind(arg, "#-")) kill_me=1;
	else{
          cerr<<"I don't understand the argument "<<arg<<" in "<<cfgfile<<endl;
          exit(1);
        }
      }
    }
    if (critical < 7){
      cerr<<"Hmm, something seems to be missing from "<<cfgfile<<endl;
      exit(1);
    }
    resize(N);
    Lz=fabs(Lz);
    for (i=begin; i<end; i++) fin>>*i;
    fin.close();
    SetProps();
    name=cfgfile.replace(cfgfile.find("_"),13,"");
  }
  else{
    cerr<<"Can't open specified cfg file: "<<cfgfile<<endl;
    exit(1);
  } 
}
//************************************************************************
void config::Dump(){
  string out_cfg=name+"_"+time2string(t)+".cfg";
  fstream fout(out_cfg.c_str(), ios::out);
  fout.fill(' '); fout.precision(12); fout.setf(ios::fixed);
  fout<<"%"<<DateTime()<<" "<<out_cfg<<endl;
  fout<<"%file format revision 04152001"<<endl;
  fout<<"#N "<<N<<endl;
  fout<<"#Nfixed "<<NFixed<<endl;
  fout<<"#Nmax "<<Nmax<<endl;
  fout<<"#material "<<mat<<endl;
  fout<<"#temperature "<<T<<endl;
  fout<<"#dimensions "<<Lx<<" "<<Ly<<" "<<Lz<<endl;
  fout<<"#time "<<t<<endl;
  if (b_masked) fout<<"#mask "<<*my_mask<<endl;
  fout<<"#----------"<<endl;
  for (i=begin; i<end; i++) fout<<*i<<endl;
}



