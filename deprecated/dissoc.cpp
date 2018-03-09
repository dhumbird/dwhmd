#include "../config.h"
#include <list>

void InfoSif3Diss(int argc, char * argv[]){
  list<string> filelist;
  int diss_zero=0;
  int diss_one=0;
  int diss_two=0;
  int diss_three=0;
  string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) filelist.push_back(arg);
  }
  list<string>::iterator f;
  int files=filelist.size()-1;
  f=filelist.begin();
  for (f++; f!=filelist.end(); f++){
    config cfg(*f);
    VPI Si=cfg.atomix(cfg.Nmax-4);
    if (Si==cfg.end) diss_three++;
    else if (Si->id==14){
      int num_F=0;
      for (SPI j=Si->nlist.begin(); j!=Si->nlist.end(); j++)
	if ((*j)->ix > Si->ix)
	  if ((*j)->id==9)
	    num_F++;
      switch (num_F){
      case 0:
	diss_three++; break;
      case 1:
	diss_two++; break;
      case 2:
	diss_one++; break;
      case 3:
	diss_zero++; break;
      }
    }
  }
  cout<<"Si-F Bond Dissociation Probability:\n";
  cout<<"0:\t"<<(float)diss_zero/files<<endl;
  cout<<"1:\t"<<(float)diss_one/files<<endl;
  cout<<"2:\t"<<(float)diss_two/files<<endl;
  cout<<"3:\t"<<(float)diss_three/files<<endl;
}
