//edep
//Perform energy deposition runs

#include "errors.h"
#include "config.h"
#include "consts.h"
#include "sysfun.h"
#include "strfun.h"
#include <fstream>

double AvgTop(config& cfg){
  int nx=(int)(cfg.Lx)+1;
  int ny=(int)(cfg.Ly)+1;
  double box[nx][ny];
  double my_top=0;
  for (int m=0; m<nx; m++) for (int n=0; n<ny; n++) box[m][n]=cfg.the_top+5;
  svector R;
  bool flag;
  for (int m=0; m<nx; m++){
    for (int n=0; n<ny; n++){
      flag=1;
      R.x=m-cfg.Lx/2;
      R.y=cfg.Ly/2-n;
      while(flag && box[m][n]>cfg.the_bottom){
	box[m][n]-=1;
	R.z=box[m][n];
	int x=(int)((R.x + cfg.Lx/2)/cfg.Lcx)%cfg.Nx;
	int y=(int)(fabs((R.y - cfg.Ly/2)/cfg.Lcy))%cfg.Ny;
	int z=(int)(fabs((R.z - cfg.the_top)/cfg.Lcz))%cfg.Nz;
	int a=x+cfg.Nx*y+cfg.Nx*cfg.Ny*z;
	set<subcell*>::iterator c;
	for (c=cfg.cells[a].nbegin(); c!=cfg.cells[a].nend(); c++)
	  for (SAI j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	    if (flag)
	      if ((R - (*j0)->R).minsqmag(cfg.Lx, cfg.Ly) > 6.25)
		flag=1;
	      else{
		flag=0;
		box[m][n]=(*j0)->R.z;
		break;
	      }
      }
    }
  }
  int total=0;
  for (int m=0; m<nx; m++) 
    for (int n=0; n<ny; n++) 
      if (box[m][n]>cfg.the_bottom){
	my_top+=box[m][n];
	total++;
      }
  return my_top/(total);
}

void MainEnergyDep(int argc, char * argv[]){
    //******************local variables************************************
  string arg; string outfile;
  string energy; string angle="0deg";
  char out_buffer [128];
  fstream fout;
  int i, k;

  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_0000-000.cfg";
  float runtime=0.4; float e=200; float Tion=300;
  float inc=0; float azi=0; float res=0.5;
  double dt=0;
  int runs=1; string id="Ar"; int seed=0; int begin=1;
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
    else if (arg=="-id") id=argv[++i];
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
  //double Llong=cfg.the_top-cfg.the_bottom;
  double zmin=-5;
  int Nz;
  float my_top;
  {
    config cfg(scfg);
    //my_top=AvgTop(cfg);
    my_top=cfg.the_top;
    Nz=(int)((cfg.the_top-cfg.the_bottom)/res)+1;  
    outfile=cfg.name+"."+Num+"."+energy+"."+angle+".dep";
  }
  vector<double> depth(Nz);
  vector<double> megabox(Nz);
  for (i=0; i<Nz; i++){
    depth[i]=0;
    megabox[i]=0;
  }
  
  //***************integration*********************************  
  for (int run=begin; run<=begin+runs-1; run++){
    cerr<<"----"<<DateTime()<<"--EDEP Run "<<run<<"-------------"<<endl;
    config cfg(scfg);
    VAI ion;
    if (dt!=0) cfg.Setdt(dt);
    if (addion) ion=cfg.AddIon(e, id, Tion, choose_inc, inc, choose_azi, azi, 
			       seed, 0, 0);
    else ion=cfg.end-1;

    double ek=ion->Ek();

    double finish=cfg.t+runtime;
    double z_cutoff=ion->R.z+0.1;
    int total_time=0;
    while (cfg.t<finish && ion->R.z <= z_cutoff && ion->R.z > cfg.the_bottom){
      if (tso) cfg.dtOptimize();
      if (cfg.u!=0) cfg.FirstVV();
      cfg.TimeStepInit();    
      cfg.ReNeighbor();    
      cfg.ForceEval();
      cfg.SecondVV();
      cfg.t+=cfg.dt;
      total_time++;
      i=(int)((my_top-ion->R.z)/res); //cerr<<i<<endl;
      if (i>0 && i<Nz) depth[i]+=ion->Ek();
    }
    for (i=0; i<Nz; i++){
      //depth[i]/=total_time;
      megabox[i]+=depth[i];
      depth[i]=0;
    }
    fout.open(outfile.c_str(), ios::out);
    for (i=0; i<Nz; i++){
      sprintf(out_buffer, "%f\t%f\n", i*res, megabox[i]/(float)run);
      fout<<out_buffer;
    }
    fout.close();
  }
}
