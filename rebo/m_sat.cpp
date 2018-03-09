#include "config.h"
#include "errors.h"

void MainSat(int argc, char * argv[]){
  string arg;
  string scfg="temp_00000000.cfg";
  short id=1;
  int rounds=10;
  int begin=1;
  float T=300;
  float runtime=0.5;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-id") id=atoi(argv[++i]);
    else if (arg=="-rounds") rounds=atoi(argv[++i]);
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-T") T=atof(argv[++i]);
  }
  config cfg(scfg);
  VAI i, target;
  float u_min;
  for (int run=begin; run<=begin+rounds-1; run++){
    cerr<<"------------"<<DateTime()<<"--------round "<<run
	<<"-----------------------"<<endl;
    cfg.SurfAcc(id);
    u_min=1000;
    for (i=cfg.begin; i!=cfg.end; i++){
      if (i->id==6 && !i->is_fixed){
	if (i->u < u_min){
	  target=i;
	  u_min=i->u;
	}
      }
    }
    VAI test=cfg.append(id);
    test->R=target->P;
    cerr<<"* Adding atom "<<test->ix<<" to "<<target->ix<<" at "<<test->R<<endl;
    cerr<<"* Well depth was "<<u_min<<" eV."<<endl;
    cfg.RunQuenched(runtime,0,0,0,0,-1,0,0,0,T);
    cfg.Dump();
  }
  cerr<<"---------------------------------------------------------------------------"<<endl;
  cerr<<DateTime()<<" Requested rounds ("<<rounds<<") completed."<<endl;
}
  
void MainPath(int argc, char* argv[]){
  string arg;
  string scfg="temp_00000000.cfg";
  int ix=-1;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-ix") ix=atoi(argv[++i]);
  }
  config cfg(scfg);
  VAI s=cfg.atomix(ix);
  if (s==cfg.end){
    cerr<<"That atom not found!\n";
    exit(1);
  }
  else{
    cfg.UPath(s);
  }
}
