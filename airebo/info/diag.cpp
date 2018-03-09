#include "../config.h"

void InfoTellme(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  string arg;
  int ix=0;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else ix=atoi(arg.c_str());
  }
  config cfg(scfg);
  VAI sweet=cfg.atomix(ix);
  if (sweet!=cfg.end){
    int bondix=-1;
    float bondr=1000;
    cout<<*sweet<<endl;
    for (VNI sw=sweet->nlist.begin(); sw!=sweet->nlist.end(); sw++){
      atom* jj = (atom*)sw->a2;
      if ((sweet->R - jj->R).minsqmag(cfg.Lx,cfg.Ly) < bondr){
	bondr=(sweet->R - jj->R).minsqmag(cfg.Lx,cfg.Ly);
	bondix=jj->ix;
      }
    }
    cout<<"Nbrs ("<<sweet->nlist.size()<<"): ";
    for (VNI sw=sweet->nlist.begin(); sw!=sweet->nlist.end(); sw++){
      atom* jj = (atom*)sw->a2;
      cout<<jj->ix;
      if (jj->ix == bondix) cout<<"*";
      cout<<" ";
    }
    cout<<endl;
    //cout<<"U: "<<cfg.U_on_i(sweet)<<" eV.";
    cout<<" Ek: "<<sweet->Ek()<<" eV."<<endl;
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

