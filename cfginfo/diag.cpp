#include "i_config.h"

void InfoTellme(int argc, char * argv[]){
  string scfg="temp_0000-000.cfg";
  string arg;
  int ix=0;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else ix=atoi(arg.c_str());
  }
  config cfg(scfg);
  VAI i=cfg.atomix(ix);
  if (i!=cfg.end){
    int bondix=-1;
    float bondr=1000;
    cout<<*i<<endl;
    for (SAI j=i->nlist.begin(); j!=i->nlist.end(); j++){
      if ((i->R - (*j)->R).minsqmag(cfg.Lx,cfg.Ly) < bondr){
	bondr=(i->R - (*j)->R).minsqmag(cfg.Lx,cfg.Ly);
	bondix=(*j)->ix;
      }
    }
    cout<<"Nbrs ("<<i->nlist.size()<<"): ";
    for (SAI j=i->nlist.begin(); j!=i->nlist.end(); j++){
      cout<<(*j)->ix;
      if ((*j)->ix == bondix) cout<<"*";
      cout<<" ";
    }
    cout<<endl;
    cout<<" Ek: "<<i->Ek()<<" eV."<<endl;
  }
  else cout<<"Atom "<<ix<<" not found!"<<endl;
}
//*********************************************************************
// void InfoNeighborDist(int argc, char * argv[]){
//   string scfg="temp_000000.cfg";
//   string arg;
//   int ix=0;
//   for (int i=2; i<argc; i++){
//     arg=argv[i];
//     if (sfind(arg, ".cfg")) scfg=arg;
//     else ix=atoi(arg.c_str());
//   }
//   config cfg(scfg);
//   VPI sweet=cfg.atomix(ix);
//   if (sweet!=cfg.end){
//     cout<<"Neighbors of "<<ix<<" ("<<sweet->nlist.size()<<"):"<<endl;
//     for (SPI sw=sweet->nlist.begin(); sw!=sweet->nlist.end(); sw++){
//       cout<<(*sw)->ix<<"\t"<<(*sw)->id<<"\t"<<
// 	((*sw)->R-sweet->R).minmag(cfg.Lx, cfg.Ly)<<endl;
//     }
//   }
//   else cout<<"Atom "<<ix<<" not found!"<<endl;
// }
//************************************************************************
// void InfoTrajectory(int argc, char * argv[]){
//   string outfile;
//   int ix=0;
//   vector<string>filelist;
//   string arg;
//   for (int j=2; j<argc; j++){
//     arg=argv[j];
//     if (sfind(arg, ".cfg")) filelist.push_back(arg);
//     else if (arg=="-ix") ix=atoi(argv[++j]);
//   }
//   int runs=filelist.size();
//   {
//     config cfg(filelist.front());
//     outfile=cfg.name+"."+itoa(ix,4)+".trj";
//   }
//   fstream fout(outfile.c_str(), ios::out);
//   for (int run=0; run<runs; run++){
//     config cfg(filelist[run]);
//     VPI I=cfg.atomix(ix);
//     if (I!=cfg.end){
//       char outbuf[256];
//       sprintf(outbuf,"%i\t%6.2f\t%6.2f\t%6.2f\n",run, 
// 	      I->R.x+I->xpass*cfg.Lx, I->R.y+I->ypass*cfg.Ly, I->R.z);
//       fout<<outbuf;
//     }
//   }
// }
//*********************************************************************

