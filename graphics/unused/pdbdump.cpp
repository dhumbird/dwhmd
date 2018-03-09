#include "../config.h"
#include <cstdio>

int main(int argc, char* argv[]){
  string arg, pdbfile;
  vector<string> filelist;
  map<int, string> name;
  name[14]="Si";
  name[9]="F";
  name[18]="Ar";
  int Nmax=0;
  
  for (int i=1; i<argc; i++){
    arg=argv[i];
    if (sfind(arg,"cfg")) filelist.push_back(arg);
  }
  if (filelist.size()>1){
    config cfg(filelist.back());
    Nmax=cfg.Nmax;
    pdbfile=cfg.name+".pdb";
  }
  else{
    config cfg(filelist.front());
    pdbfile=cfg.name+"_"+time2string(cfg.t)+".pdb";
  }
  
  int ione=1;
  float fone=1; float fzero=0;

  fstream pdb(pdbfile.c_str(), ios::out);
  pdb<<"HEADER  PROTEIN\n";
  pdb<<"COMPND  "<<pdbfile<<endl;
  pdb<<"AUTHOR  JOE MAMMA\n";
  char outbuf[256];

  for (int i=0; i<filelist.size(); i++){
    if (filelist.size()>1){
      sprintf(outbuf, "MODEL %4d\n", i+1);
      pdb<<outbuf;
    }
    int serial=-1;
    config cfg(filelist[i]);
    for (VPI s=cfg.begin; s<cfg.end; s++){
      if (filelist.size()>1)
	while(serial+1 < s->ix){
	  sprintf(outbuf, "ATOM  %5d\n", ++serial);
	  pdb<<outbuf;
	}
      
      //record name (string 6)
      sprintf(outbuf, "ATOM  "); pdb<<outbuf;
      
      //serial number (integer 5)
      sprintf(outbuf, "%5d", ++serial); pdb<<outbuf;
      
      //skip column 12
      pdb<<" ";
      
      //atom name (string 4)
      sprintf (outbuf, "%4s", name[s->id].c_str()); pdb<<outbuf;
      
      //alternate location indicator (char);
      //residue name (string 3); skip column 21; chain identifier (char)
      pdb<<" UNK  ";
      
      //residue sequence number (integer 4)
      sprintf (outbuf, "%4d", s->ix); pdb<<outbuf;

      //code for insertion of residues (char); skip columns 28-30
      pdb<<"    ";

      //x,y,z coord (angstroms) (8.3, 8.3, 8.3)
      sprintf (outbuf, "%8.3f%8.3f%8.3f",s->R.x,s->R.y,s->R.z); pdb<<outbuf;
      
      //occupancy factor (6.2) --i use this for color defs
      //if ($colordef){
      //  printf PDB "%6.2f",$color{$ix};
      //}
      //    else{
      sprintf (outbuf, "%6.2f", fone); pdb<<outbuf;
      //}
      
      //temperature factor (6.2)
      if (!s->is_fixed)
	if (s->Ek() < 9){
	  sprintf(outbuf, "%6.2f", (9-s->Ek())/9.0); 
	  pdb<<outbuf;
	}
	else{
	  sprintf(outbuf, "%6.2f", fzero); 
	  pdb<<outbuf;
	}
      else{
	sprintf (outbuf, "%6.2f",fone); 
	pdb<<outbuf;
      }
      
      //end of record
      pdb<<endl;
    }
    if (filelist.size()>1){
      while (serial < Nmax-1){
	sprintf(outbuf, "ATOM  %5d\n", ++serial);
	pdb<<outbuf;
      }
      pdb<<"ENDMDL\n";
    }
  }
  pdb<<"END\n";
  cout<<pdbfile;
}
