#include "config.h"
#include "errors.h"
#include "consts.h"
#include "sysfun.h"
#include "strfun.h"
#include <list>

void MainGrowFilm(int argc, char * argv[]){
  
  int i;
  string arg;
  
  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_000000.cfg";
  float runtime=0.2; float e_ion=0; float Tion=300;
  float inc=0; float azi=0; float T=300; 
  double dt=0; double bhb_dt=0;
  int cfgevery=0; int runs=1; int ion_id=18; int begin=1;
  bool printe=0; bool dump_end=0; bool tso=0;
  bool choose_inc=0; bool choose_azi=0; bool dump_only=0;
  short cfgthermo=1; bool pickup=0;
  bool refresh=0; bool dump_each=0;
  list<int> seeds;
  string species="unknown";
  string mix=""; float ionevery=0;
  map<string,float> ratio;
  for (i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-cfgevery") cfgevery=atoi(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-printe") printe=1;
    else if (arg=="-runs"||arg=="-rounds") runs=atoi(argv[++i]);
    else if (arg=="-e"||arg=="-eion") e_ion=atof(argv[++i]);
    else if (arg=="-T") T=atof(argv[++i]);
    else if (arg=="-Tion") Tion=atof(argv[++i]);
    else if (arg=="-type") species=argv[++i];
    else if (arg=="-bhb_dt") bhb_dt=atof(argv[++i]);
    else if (arg=="-pickup") pickup=1;
    else if (arg=="-tso") tso=1;
    else if (arg=="-nocfgthermo") cfgthermo=0;
    else if (arg=="-refresh") refresh=1;
    else if (arg=="-ionevery") ionevery=atof(argv[++i]);
    else if (arg=="-dump"){
      arg=argv[++i];
      if (arg=="end") dump_end=1;
      else if (arg=="each") dump_each=1;
      else if (arg=="only") dump_only=1;
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
    else if (arg=="-seeds")
      for (int j=0; j<runs; j++)
	seeds.push_back(atoi(argv[++i]));
    else CmdError(arg.c_str());
  }

  if (seeds.empty()) 
    for (int j=0; j<runs; j++)
      seeds.push_back(0);
  //*******************************************************************
  
  if (ionevery != 0) ionevery=1/ionevery;
  config cfg(scfg);
  cfg.u=cfg.Uonly();
  list<int>::iterator seed=seeds.begin();
  int nsave=-1;
  for (int run=begin; run<=begin+runs-1; run++){
    cerr<<"------------"<<DateTime()<<"--------run "<<run
 	<<"-----------------------"<<endl;
    svector P;
    vector<int> dbC;
    vector<int> Hs;
    for (VAI s=cfg.begin; s<cfg.end; s++){
      if (s->id==1) Hs.push_back(s->ix);
      if (s->id == 6 && !s->is_fixed && s->Nt < 3.3){
	dbC.push_back(s->ix);
	if (s->Nt < 2.3) dbC.push_back(s->ix);
	if (s->Nt < 1.3) dbC.push_back(s->ix);
      }
    }
    //pick one to resaturate
    srand(TimeSeed()-1);
    if (dbC.size() > 0){
      bool res=1;
      do{
	int resat=(int)floor(rand01()*dbC.size());
	atom* newatom=&(*cfg.append(1));
	VAI s=cfg.atomix(dbC[resat]);
	s->u=1000;
	float bl=1.1;
	float uu;
	for (int c=0; c<500; c++){
	  float phi=rand01()*PI;
	  float theta=rand01()*2*PI;
	  P.set(bl*sin(phi)*cos(theta), bl*sin(phi)*sin(theta), bl*cos(phi));
	  P+=s->R;
	  P.minimg(cfg.Lx,cfg.Ly);
	  newatom->R=P;
	  uu=cfg.U_on_i(newatom);
	  if (uu < s->u){
	    s->u=uu;
	    s->P=P;
	  }
	}
	//attempt to resat.
	newatom->R=s->P;
	cfg.RunQuenched(0.01,0,0,0,0,-1,1,1,1,T);
	if (newatom->NC > 0.5 ){
	  if (nsave!=-1){
	    newatom->ix=nsave;
	    cfg.Nmax--;
	  }	  
	  cerr<<"* Put an H (# "<<newatom->ix<<") on Carbon "<<s->ix<<endl;
	  dbC.erase(dbC.begin()+resat);
	  res=1;
	} 
	else{
	  cfg.erase(cfg.atomix(newatom->ix));
	  res=0;
	}
      }while (!res);
    }
    //remove random H's until there are 5% db
    int dbcsize=dbC.size();
    float target=0.05*Hs.size();
    vector<int> Herase;
    while (dbcsize < target){
      int take=(int)floor(rand01()*Hs.size());
      atom* anH = &(*cfg.atomix(Hs[take]));
      //do a line-of-sight test (2 Ang. cylinder)
      svector dir;
      float zmax=cfg.the_top;
      for (int c=0; c<500; c++){
	float phi=rand01()*PI/2;
	float theta=rand01()*2*PI;
	dir.set(1, theta, phi);
	dir.Sph2Rec();
	P=anH->R;
	bool goodpath=1;
	do{
	  P+=dir/10;
	  P.minimg(cfg.Lx,cfg.Ly);
	  if (P.z > zmax) break;
	  subcell* sc=cfg.WhichCell(P);
	  set<subcell*>::iterator c;
	  for (c=sc->nbegin(); c!=sc->nend(); c++)
	    for (SAI j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	      if (&(*anH)!=(*j0)){
		if ((P - (*j0)->R).minsqmag(cfg.Lx, cfg.Ly) < 1){
		  goodpath=0;
		}
	      }
	}while (goodpath);
	if (goodpath){
	  Herase.push_back(anH->ix);
	  Hs.erase(Hs.begin()+take);
	  dbcsize++;
	  break;
	}
      }
    }
    if (Herase.size() > 0){
      for (vector<int>::iterator vi=Herase.begin(); vi!=Herase.end(); vi++){
	cfg.erase(cfg.atomix(*vi));
	cerr<<"* Removed Hydrogen "<<*vi<<endl;
	cfg.RunQuenched(0.01,0,0,0,0,-1,1,1,1,T);
	nsave=*vi;
      }
    }
    else nsave=-1;
    dbC.clear();
    for (VAI s=cfg.begin; s<cfg.end; s++){
      if (s->id == 6 && !s->is_fixed && s->Nt < 3.3){
	dbC.push_back(s->ix);
	if (s->Nt < 2.3) dbC.push_back(s->ix);
	if (s->Nt < 1.3) dbC.push_back(s->ix);
      }
    }
    double ipmin=1000;
    VAI I;
    species="CH3";
    int No=cfg.N;

    I=cfg.AddIon(-1, species, Tion, choose_inc, inc, choose_azi, 
		 azi, *seed, 0);

    svector Vu;
    for (VAI s=I; s<cfg.end; s++) Vu+=s->V*s->m;
    Vu/=Vu.mag();
    int coun=-1;
    while (ipmin > 16 && coun < 1000){
      //I'm going to translate the ion around laterally until it's good.
      svector Rrand(rand01(), rand01(), 0);
      for (VAI s=I; s<cfg.end; s++){
	s->R+=Rrand;
	s->R.minimg(cfg.Lx, cfg.Ly);
      }
      for (vector<int>::iterator vii=dbC.begin(); vii!=dbC.end(); vii++){
	VAI targ=cfg.atomix(*vii);
	double impt=(targ->R.z-I->R.z)/Vu.z;
	svector Rproj(I->R.x+Vu.x*impt, I->R.y+Vu.y*impt, targ->V.z);
	Rproj.minimg(cfg.Lx, cfg.Ly);
	double ipar=(Rproj-targ->R).minsqmag(cfg.Lx,cfg.Ly);
	if (ipar<ipmin) ipmin=ipar;
      }
      coun++;
    }

    cfg.Partition();
    //minimum z-position
    bool got_neighbors=0;
    SAI sai_j;
    set<subcell*>::iterator c;
    while (!got_neighbors){
      for (VAI s=I; s<cfg.end; s++){
        subcell* sc=cfg.WhichCell(s->R);
        for (c=sc->nbegin(); c!=sc->nend(); c++)
          for (sai_j=(*c)->abegin(); sai_j!=(*c)->aend(); sai_j++)
            if (!got_neighbors)
              if((*sai_j)->ix < I->ix){
                if ((s->R - (*sai_j)->R).minsqmag(cfg.Lx, cfg.Ly) < 4)
                  got_neighbors=1;
	      }
      }
      if (!got_neighbors) 
	for (VAI s=I; s<cfg.end; s++){ 
	  s->R+=0.01*Vu; s->R.minimg(cfg.Lx,cfg.Ly);
	  cfg.Partition();
	}
    }
    for (VAI s=I; s<cfg.end; s++){
      s->R-=0.01*Vu;
      s->R.minimg(cfg.Lx,cfg.Ly);
    }
    cfg.Partition(); cfg.ReNeighbor();

    cerr<<"* Translation settled on "<<I->R;
    if (coun==1000) cerr<<" (Gave up.)"<<endl;
    else cerr<<" ("<<coun<<" skipped)"<<endl;

    for (VAI s=cfg.begin; s<cfg.end; s++) s->P=s->R;
    int Hix=I->ix;
    double nbegin=I->NC;

    bool quit=0;
    while (!quit){
      //step 10 fs at a time.
      cfg.Run(0.01, 0, dt, printe, 0, -1, 0, tso,0);
      I=cfg.atomix(Hix);
      if (I==cfg.end){
	if (species=="CH3") cfg.Nmax-=4;
	quit=1;
	if (cfg.N==No) cerr<<"* Species was scattered.\n";
	else{ 
	  cerr<<"* Species abstracted something.\n";
	  if (dump_only) cfg.Dump();
	}
	if (fabs(1-T/cfg.T) > 0.05)
	  cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0);
      }
      else if (I->NC > 0.5){
	cerr<<"* Species appears to be bound. "<<I->R<<endl;
	cerr<<"  Begin: "<<nbegin<<" End: "<<I->NC<<endl;
	quit=1;
	cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0.5);
	if (dump_only) cfg.Dump();
      }
    }
    cerr<<"* Impact complete. Time: "<<cfg.t<<endl;
    if (dump_each) cfg.Dump();
    seed++;
  }
  if (dump_only || dump_end) cfg.Dump();
  cerr<<"---------------------------------------------------------------------------"<<endl;
  cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
}

