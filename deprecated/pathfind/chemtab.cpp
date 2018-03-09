#include "chemtab.h"

map<short, cheminfo> chemtab;

void ChemDef(){  

  //silicon (14)
  chemtab[14].mass=M_14;
  chemtab[14].cutoff=RC_14;
  chemtab[14].three_body=TB_14;
  
  //argon (18)
  chemtab[18].mass=M_18;
  chemtab[18].cutoff=RC_18;
  chemtab[18].three_body=TB_18;

  //helium (2)
  chemtab[2].mass=M_2;
  chemtab[2].cutoff=RC_2;
  chemtab[2].three_body=TB_2;

  //fluorine (9)
  chemtab[9].mass=M_9;
  chemtab[9].cutoff=RC_9;
  chemtab[9].three_body=TB_9;

  //neon (10)
  chemtab[10].mass=M_10;
  chemtab[10].cutoff=RC_10;
  chemtab[10].three_body=TB_10;

  //chlorine (17)
  chemtab[17].mass=M_17;
  chemtab[17].cutoff=RC_17;
  chemtab[17].three_body=TB_17;

  //krypton (36)
  chemtab[36].mass=M_36;
  chemtab[36].cutoff=RC_36;
  chemtab[36].three_body=TB_36;

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
