#include "../config.h"
#include <list>
#include <map>

void SiF3Statistics (int argc, char * argv[]){
  list<string> filelist;
  map<int, list<int> > partners;
  string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) filelist.push_back(arg);
  }
  string outfile;
  {
    config cfg(filelist.back());
    outfile=cfg.name+".fstat";
  }
  fstream fout(outfile.c_str(), ios::out);
  list<string>::iterator f;
  int impact=0;
  for (f=filelist.begin(); f!=filelist.end(); f++){
    impact++;
    config cfg(*f);
    int totalF=0;
    int boundF=0;
    for (VPI i=cfg.begin; i<cfg.end; i++){
      if (i->id==9){
	totalF++;
	int partner=-1;
	double d1=1000; double d2=1000;
	for (SPI s=i->nlist.begin(); s!=i->nlist.end(); s++){
	  d2=d1;
	  d1 = d1 <? ((*s)->R - i->R).minmag(cfg.Lx,cfg.Ly);
	  if (d2>d1) partner=(*s)->ix;
	}
	if (partner - i->ix >=3) boundF++;
      }
    }
    if (totalF>0){
      fout<<impact<<"\t"<<(float)boundF/totalF<<endl;
    }
    else fout<<impact<<"\t0"<<endl;
  }
}



