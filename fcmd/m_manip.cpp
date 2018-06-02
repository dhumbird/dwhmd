#include "config.h"
#include "errors.h"

void MainReset(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
  }
  config cfg(scfg);
  cfg.reset();
  cfg.Dump();
  exit(0);
}

void MainDelete(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  int ix=0; string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else ix=atoi(arg.c_str());
  }
  config cfg(scfg);
  if (cfg.erase(cfg.atomix(ix))) cfg.Dump();
  else cerr<<"Atom "<<ix<<" not found!\n";
  exit(0);
}

void MainShiftOrigin(int argc, char * argv[])
{
  short ix=-1; string arg;
  float x=0; float y=0;
  vector<string> filelist;
  vector<string>::iterator f;
  for (int i=2; i<argc; i++)
  {
    arg=argv[i];
    if (sfind(arg, ".cfg")) filelist.push_back(arg);
    else if (arg=="-x") x=atof(argv[++i]);
    else if (arg=="-y") y=atof(argv[++i]);
    else if (arg=="-ix") ix=atoi(argv[++i]);
    else CmdError(arg.c_str());
  }
  for (f=filelist.begin(); f!=filelist.end(); f++)
  {
    config cfg(*f);
    svector R;
    if (ix!=-1)
    {
      if (cfg.atomix(ix)!=cfg.end)
      {
        R=cfg.atomix(ix)->R;
      }
      else
      {
	      cerr<<"That atom not found!"<<endl;
	      exit(0);
      }
    }
    else
    {
      R.x=x;
      R.y=y;
    }
    R.z=0;
    for (VAI i=cfg.begin; i!=cfg.end; i++){
      i->R-=R; i->R.minimg(cfg.Lx, cfg.Ly);
    }
    cfg.Dump();
  }
  exit(0);
}

void MainRandomize(int argc, char * argv[]){
  short ix=-1; string arg;
  float T=300; string scfg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-T") T=atof(argv[++i]);
    else CmdError(arg.c_str());
  }
  config cfg(scfg);
  cfg.Maxwell(T,0);
  cfg.Dump();
  exit(0);
}

void MainMisc(int argc, char * argv[]){
  string scfg="temp_0000-000.cfg";
  string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else CmdError(arg.c_str());
  }
  config cfg(scfg);
  cfg.Sweep();
  //cfg.Dump();
  exit(0);
}

// void MainPass(int argc, char * argv[]){
//   string scfg="temp_0000-000.cfg";
//   string id="H"; string arg;
//   for (int i=2; i<argc; i++){
//     arg=argv[i];
//     if (sfind(arg, ".cfg")) scfg=arg;
//     else if (arg=="-id") id=argv[++i];
//     else CmdError(arg.c_str());
//   }
//   config cfg(scfg);
//   cfg.Passivate(id);
//   cfg.Dump();
//   exit(0);
// }
