#include "chemtab.h"

map<short, cheminfo> chemtab;

void ChemDef(){  

  //silicon (14)
  chemtab[14].mass=M_SI;
  chemtab[14].cutoff=RC_SI;
  chemtab[14].three_body=TB_SI;
  
  //argon (18)
  chemtab[18].mass=M_AR;
  chemtab[18].cutoff=RC_AR;
  chemtab[18].three_body=TB_AR;

  if (FileExists("chem.def")){
    ifstream fin("chem.def", ios::in);
    string s;
    short id=0;
    float m=0; float rc=0;
    bool tb=0;
    while(!fin.eof()){
      getline(fin,s);
      if (s.find("#")!=0){
	istringstream ist(s);
	ist>>id>>m>>rc>>tb;
	chemtab[id].mass=m;
	chemtab[id].cutoff=rc;
	chemtab[id].three_body=tb;
      }
    }
  }
}
