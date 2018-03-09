#include "config.h"
#include "errors.h"

void MainSat(int argc, char * argv[]){
  string arg;
  string scfg="temp_00000000.cfg";
  short id=9;
  int rounds=10;
  int begin=1;
  float T=300;
  float aftertime=0.5;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-id") id=atoi(argv[++i]);
    else if (arg=="-rounds") rounds=atoi(argv[++i]);
    else if (arg=="-time") aftertime=atof(argv[++i]);
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-T") T=atof(argv[++i]);
  }
    
  for (int run=begin; run<=begin+rounds-1; run++){
    int seed=TimeSeed();
    srand(seed);
    cerr<<"------"<<DateTime()<<"-----round "<<run
	<<"---seed "<<seed<<"------------"<<endl;
    bool kept=0;
    while(!kept){
      set<subcell*>::iterator c;
      config cfg(scfg);
      float zmin=cfg.MinZ();
      //pick a random trajectory
      VPI newatom=cfg.append(id);
      newatom->R.set(rand01()*fabs(cfg.Lx), rand01()*fabs(cfg.Ly), 0);
      newatom->R.minimg(cfg.Lx, cfg.Ly);
      //assign z to be top plane
      newatom->R.z=cfg.MaxZ()+sqrt(newatom->rc)+0.001;
      cfg.Partition();
      //pick angles
      svector unitlen=-1*svector(1, rand01()*2*PI, rand01()*0.5*PI);
      unitlen.Sph2Rec();
      map<float,int> si_nbrs;
      float rij=0;
      bool goodpath=1;
      do{
	si_nbrs.clear();
	newatom->R+=unitlen/10;
	newatom->R.minimg(cfg.Lx,cfg.Ly);
	if (newatom->R.z < zmin) goodpath=0;
	subcell* sc=cfg.WhichCell(newatom->R);
	for (c=sc->nbegin(); c!=sc->nend(); c++)
	  for (SPI j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	    if (&(*newatom)!=(*j0)){
	      rij=(newatom->R - (*j0)->R).minsqmag(cfg.Lx, cfg.Ly);
	      if (rij < 1) goodpath=0;
	      if (rij < 6.25)
		if ((*j0)->id==14) si_nbrs[rij]=(*j0)->ix;
	    }
      }while (goodpath && si_nbrs.size()<2);
      
      double u=1000; int a1,a2;
      if (goodpath){
	map<float,int>::iterator mfi=si_nbrs.begin();
	a1=mfi->second; mfi++;
	a2=mfi->second;
	u=cfg.U_on_i(newatom);
      }
      if (u < 6){
	cerr<<"* Adding atom "<<newatom->ix<<" at "<<newatom->R<<endl;
	cerr<<"* Well depth was "<<u<<" eV btwn atoms "
	    <<a1<<" and "<<a2<<endl;
	newatom->is_fixed=1;
	cfg.RunQuenched(0.2,0,0,0,0,-1,0,0,1,T);
	newatom->is_fixed=0;
	cfg.RunQuenched(0.3,0,0,0,0,-1,0,0,1,T);
	cfg.CheckSput();
	cfg.CheckDepth();
	cfg.Dump();
	scfg=cfg.name+"_"+time2string(cfg.t)+".cfg";
	kept=1;
	break;
      }
      else{
	cfg.erase(newatom);
	cfg.Nmax--;
	cfg.TimeStepInit();
      }
    }
  }
  cerr<<"---------------------------------------------------------------------------"<<endl;
  cerr<<DateTime()<<" Requested rounds ("<<rounds<<") completed."<<endl;
}
  
