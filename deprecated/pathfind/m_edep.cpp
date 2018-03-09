//edep
//Perform energy deposition runs

#include "errors.h"
#include "config.h"
#include "consts.h"
#include "sysfun.h"
#include "strfun.h"
#include <fstream>

void MainEnergyDep(int argc, char * argv[]){
  
  //******************local variables************************************
  string arg; string outfile;
  string energy; string angle="0deg";
  char out_buffer [128];
  fstream fout;
  int i, k;

  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_000000.cfg";
  float runtime=0.4; float e=200; float Tion=300;
  float inc=0; float azi=0; float res=0.5;
  double dt=0;
  int runs=1; int id=18; int seed=0; int begin=1;
  bool addion=1; bool tso=0;
  bool choose_inc=0; bool choose_azi=0;
  int num=1;
  
  for (i=2; i<argc; i++){
    arg=argv[i];
    if (arg.find(".cfg")!=string::npos){
      scfg=arg;
    }
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-runs") runs=atoi(argv[++i]);
    else if (arg=="-e")  e=atof(argv[++i]);
    else if (arg=="-T") Tion=atof(argv[++i]);
    else if (arg=="-tso") tso=1;
    else if (arg=="-id") id=atoi(argv[++i]);
    else if (arg=="-num") num=atoi(argv[++i]);
    else if (arg=="-seed") seed=atoi(argv[++i]);
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-res") res=atof(argv[++i]);
    else if (arg=="-noion") addion=0;
    else if (arg=="-inc"){
      arg=argv[++i];
      if (arg=="rand"){
	choose_inc=1;
	angle="randeg";
      }
      else{
	inc=atof(arg.c_str());
	angle=ftoa(inc, 3)+"deg";
      }
    }
    else if (arg=="-azi"){
      arg=argv[++i];
      if (arg=="rand") choose_azi=1;
      else azi=atof(arg.c_str());
    }
    else CmdError(arg.c_str());
  }
  string Num=itoa(num,2);
  //**********************begin code****************************************
  if (addion){
    energy=ftoa(e, 3);  
    energy+="eV";
  }  

  //set up the spatial array.
  double zmin=-5;
  int Nz;
  float my_top;
  {
    config cfg(scfg);
    my_top=cfg.AvgTop();
    Nz=(int)(cfg.Lz_full/res)+1;  
    outfile=cfg.name+"."+Num+"."+energy+"."+angle+".dep";
  }
  float box[Nz];
  float megabox[Nz];
  for (i=0; i<Nz; i++){
    box[i]=0;
    megabox[i]=0;
  }
  
  //***************integration*********************************  
  for (int run=begin; run<=begin+runs-1; run++){
    cerr<<"----"<<DateTime()<<"--EDEP Run "<<run<<"-------------"<<endl;
    config cfg(scfg);
    VPI I,ion;
    if (dt!=0) cfg.Setdt(dt);
    if (addion) ion=cfg.AddIon(e, id, Tion, choose_inc, inc, choose_azi, azi, 
			       seed, 0);
    else ion=cfg.end-1;

    ion->P.z=ion->P.y=ion->P.x=ion->Ek();

    double finish=cfg.t+runtime;
    double z_cutoff=ion->R.z+0.1;
    
    //here's the integration.
    while (cfg.t<finish && ion->R.z <= z_cutoff){
      if (tso) cfg.dtOptimize();
      if (cfg.u!=0) cfg.FirstVV();
      cfg.TimeStepInit();    
      cfg.ReNeighbor();    
      cfg.ForceEval();
      cfg.SecondVV();
      cfg.t+=cfg.dt;
      ion->P.x=ion->P.y-ion->Ek();
      ion->P.y=ion->Ek();
      i=(int)((my_top-ion->R.z)/res);
      if (i<Nz) box[i]+=ion->P.x;
    }
    float ek_dep=ion->P.z-ion->P.y;
    float ek_sum=0;
    for (i=0; i<Nz; i++)
      ek_sum+=box[i];
    float ek_rat=ek_dep/ek_sum;
    for (i=0; i<Nz; i++){
      megabox[i]+=box[i]*ek_rat;
      box[i]=0;
    }
  }
  fout.open(outfile.c_str(), ios::out);
  for (i=0; i<Nz; i++){
    sprintf(out_buffer, "%f\t%f\n", i*res, megabox[i]/runs);
    fout<<out_buffer;
  }
  fout.close();
}




