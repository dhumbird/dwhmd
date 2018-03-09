#include "bombard.h"

int IonBombard(int argc, char * argv[]){
  
  int i;
  string arg;
  
  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_000000.cfg";
  float runtime=0.2; float e=200; float Tion=300;
  float inc=0; float azi=0; float T=300; bool confine=0;
  bool round=0; float conf=0;
  double dt=0; double bhb_dt=0;
  int cfgevery=0; int runs=1; int id=18; int part=0; int begin=1;
  bool tso=0; bool printe=0; bool dump_end=0; bool traj=0;
  bool choose_inc=0; bool choose_azi=0; bool no_noise=0; bool choose_pos=1;
  short cfgthermo=1;
  bool refresh=0; bool dump_each=0;

  for (i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-cfgevery") cfgevery=atoi(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-tso") tso=1;
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-printe") printe=1;
    else if (arg=="-runs") runs=atoi(argv[++i]);
    else if (arg=="-e") e=atof(argv[++i]);
    else if (arg=="-T") T=atof(argv[++i]);
    else if (arg=="-Tion") Tion=atof(argv[++i]);
    else if (arg=="-id") id=atoi(argv[++i]);
    else if (arg=="-bhb_dt") bhb_dt=atof(argv[++i]);
    else if (arg=="-norand") choose_pos=0;
    else if (arg=="-no_noise") no_noise=1;
    else if (arg=="-part") part=atoi(argv[++i]);
    else if (arg=="-traj") traj=1;
    else if (arg=="-confine"){
      confine=1;
      arg=argv[++i];
      if (arg=="round") round=1;
      else if (arg=="square") round=0;
      else CmdError(arg.c_str());
      conf=atof(argv[++i]);
    }
    else if (arg=="-nocfgthermo") cfgthermo=0;
    else if (arg=="-refresh") refresh=1;
    else if (arg=="-dump"){
      arg=argv[++i];
      if (arg=="end") dump_end=1;
      else if (arg=="each") dump_each=1;
      else cfgevery=atoi(arg.c_str());
    }
    else if (arg=="-inc"){
      arg=argv[++i];
      if (arg=="rand") choose_inc=1;
      else inc=atof(arg.c_str());
    }
    else if (arg=="-azi"){
      arg=argv[++i];
      if (arg=="rand") choose_azi=1;
      else azi=atof(arg.c_str());
    }
    else CmdError(arg.c_str());
  }
  //**************************************************************************
  
  if (!refresh){
    config cfg(scfg);
    for (int run=begin; run<=begin+runs-1; run++){
      cerr<<"------------"<<DateTime()<<"--------run "<<run
	  <<"-----------------------"<<endl;
      cfg.AddIon(e, id, Tion, choose_inc, inc, choose_azi, azi, 0, 0, no_noise,
		 choose_pos, confine, round, conf);
      VPI I;
      if (traj) for (I=cfg.begin; I<cfg.end; I++) I->P=I->R;
      if (part>0) AddMultiIons(part, &cfg);
      cfg.Run(runtime, cfgevery, dt, printe, tso, 0);
      if (traj){
	string tr_file=cfg.name+"_"+time2string(cfg.t)+".trj";
	ofstream tr_out(tr_file.c_str(), ios::out);
	for (I=cfg.begin; I<cfg.end; I++)
	  if (!I->is_fixed)
	    if (I->id==cfg.mat)
	      if ((I->P - I->R).minsqmag(cfg.Lx, cfg.Ly, cfg.Lz) > 0.5)
		tr_out<<I->ix<<"\t"<<I->P<<"\t "<<I->R<<endl;
      }
      cfg.Thermo(T, bhb_dt, cfgthermo*cfgevery, dt, printe, tso, 0);
      if(part>0) DelMultiIons(part, &cfg);
      cfg.DelIon(0);
      if (dump_each) cfg.Dump(0);
    }
    if (dump_end) cfg.Dump(0);
    cerr<<"---------------------------------------------------------------------------"<<endl;
    cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
    cerr<<cfg.Nsput<<" atoms were sputtered."<<endl;
    return(0);
  }

  else
    for (int run=begin; run<=begin+runs-1; run++){
      cerr<<"------------"<<DateTime()<<"--------run "<<run
	  <<"-----------------------"<<endl;
      config cfg(scfg);
      cfg.AddIon(e, id, Tion, choose_inc, inc, choose_azi, azi, 0, 0, no_noise,
		 choose_pos, confine, round, conf);
      VPI I;
      if (traj) for (I=cfg.begin; I<cfg.end; I++) I->P=I->R;
      cfg.Run(runtime, 0, dt, printe, tso, 0);
      if (traj){
	string tr_file=cfg.name+"_"+time2string(cfg.t)+".trj";
	ofstream tr_out(tr_file.c_str(), ios::out);
	for (I=cfg.begin; I<cfg.end; I++)
	  if (!I->is_fixed)
	    if (I->id==cfg.mat)
	      if ((I->P - I->R).minsqmag(cfg.Lx, cfg.Ly, cfg.Lz) > 0.5)
		tr_out<<I->ix<<"\t"<<I->P<<"\t "<<I->R<<endl;
      }
      if (dump_each){
	cfg.DelIon(0);
	cfg.Dump(0);
      }
    }
}

//****************************************************************************
void AddMultiIons(int part, config * cfg){

  /*this code divides the surface into even partitions and copies the
  end atom (presumably an ion) identically above each of them. I don't
  think it will work for part!=3 (ie 9 partitions).*/

  VPI I, ion; 
  ion=cfg->end-1;
  for (int p=1; p<part; p++){
    cfg->CopyAtom(ion);
    I=cfg->end-1;
    I->R.x -= ((float)p)*cfg->Lx/part;
    I->R.minimg(cfg->Lx, cfg->Ly, cfg->Lz);
    cfg->CopyAtom(ion);
    I=cfg->end-1;
    I->R.y -= ((float)p)*cfg->Ly/part;
    I->R.minimg(cfg->Lx, cfg->Ly, cfg->Lz);
    cfg->CopyAtom(ion);
    I=cfg->end-1;
    I->R.x -= ((float)p)*cfg->Lx/part;
    I->R.y -= ((float)p)*cfg->Ly/part;
    I->R.minimg(cfg->Lx, cfg->Ly, cfg->Lz);
    if(p>1){
      cfg->CopyAtom(ion);
      I=cfg->end-1;
      I->R.x -= ((float)(p-1))*cfg->Lx/part;
      I->R.y -= ((float)p)*cfg->Ly/part;
      I->R.minimg(cfg->Lx, cfg->Ly, cfg->Lz);
      cfg->CopyAtom(ion);
      I=cfg->end-1;
      I->R.x -= ((float)p)*cfg->Lx/part;
      I->R.y -= ((float)(p-1))*cfg->Ly/part;
      I->R.minimg(cfg->Lx, cfg->Ly, cfg->Lz);
    }
  }
}
//************************************************************************
void DelMultiIons(int part, config * cfg){
  for (int p=0; p<part*part-1; p++) cfg->DelIon(1);
}
