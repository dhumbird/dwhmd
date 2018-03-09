#include "edep.h"

class fbox{
public:
  float f;
  int count;
  void clear(){f=0; count=0;}
  float avg(){return (count >0) ? f/count : 0;}
};


int EnergyDep(int argc, char * argv[]){
  
  //******************local variables************************************
  string arg; string outfile;
  string host=Hostname();
  host=host[0];
  string energy; string angle="0deg";
  char out_buffer [128];
  fstream fout;
  int i, k;

  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_000000.cfg";
  float runtime=0.2; float e=200; float Tion=300;
  float inc=0; float azi=0; float noisetime=0.1; float res=0.5;
  double dt=0;
  int runs=1; int id=18; int seed=0; int begin=1;
  bool tso=0; bool addion=1;
  bool choose_inc=0; bool choose_azi=0;
  
  for (i=2; i<argc; i++){
    arg=argv[i];
    if (arg.find(".cfg")!=string::npos){
      scfg=arg;
    }
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-tso") tso=1;
    else if (arg=="-runs") runs=atoi(argv[++i]);
    else if (arg=="-e")  e=atof(argv[++i]);
    else if (arg=="-T") Tion=atof(argv[++i]);
    else if (arg=="-id") id=atoi(argv[++i]);
    else if (arg=="-seed") seed=atoi(argv[++i]);
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-noisetime") noisetime=atof(argv[++i]);
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
  
  //**********************begin code****************************************
  if (addion){
    energy=ftoa(e, 3);  
    energy+="eV";
  }  

  //set up the spatial array.
  double zmin=-5;
  int Nr, Nz;
  {
    config cfg(scfg);
    Nr=(int)((sqrt(pow(cfg.Lx/2,2)+pow(cfg.Ly/2,2)))/res)+1;
    Nz=(int)((fabs(cfg.Lz)-zmin)/res)+1;  
    outfile=cfg.name+"."+host+"."+energy+"."+angle+".dep";
  }
  fbox box[Nr][Nz];
  for (i=0; i<Nr; i++) for (k=0; k<Nz; k++) box[i][k].clear();

  //***************integration*********************************  
  for (int run=begin; run<=begin+runs-1; run++){
   
    if(addion) if (seed==0) seed=(unsigned)clock()+(unsigned)time(NULL);      

    config cfg(scfg);
    VPI I;
    if (dt!=0) cfg.Setdt(dt);

    //take note of origins
    for (I=cfg.begin; I<cfg.end; I++) I->P=I->R;

    if (addion) cfg.AddIon(e, id, Tion, choose_inc, inc, choose_azi,
			   azi, seed, 1, 0, 1, 0, 0, 0);
      
    VPI ion=cfg.end-1;
    svector O=ion->R;
    double start=cfg.t;
    double finish=start+runtime;
    double ek_cutoff=ion->Ek()*0.05 <? 5;
    double z_cutoff=ion->R.z+0.1;
    double zmax=0;
    
    for (I=cfg.begin; I<ion; I++) zmax = zmax >? I->R.z;
    
    //find point of impact in zmax plane
    double ximp=O.x+(zmax-O.z)*ion->V.x/ion->V.z;
    double yimp=O.y+(zmax-O.z)*ion->V.y/ion->V.z;
    svector Impact(ximp, yimp, zmax);
    O-=Impact;

   
    //construct Euler matrix
    double A[9];
    {
      double phi=-1*atan(O.x/O.y);
      double theta=(ion->V.y >0)?
	PI+acos(O.z/O.mag()):
	PI-acos(O.z/O.mag());
      double q0=cos(0.5*theta)*cos(0.5*phi);
      double q1=sin(0.5*theta)*cos(0.5*phi);
      double q2=sin(0.5*theta)*sin(0.5*phi);
      double q3=cos(0.5*theta)*sin(0.5*phi);
      A[0]=q0*q0+q1*q1-q2*q2-q3*q3;
      A[1]=2*(q1*q2+q0*q3);
      A[2]=2*(q1*q3-q0*q2);
      A[3]=2*(q1*q2-q0*q3);
      A[4]=q0*q0-q1*q1+q2*q2-q3*q3;
      A[5]=2*(q2*q3+q0*q1);
      A[6]=2*(q1*q3+q0*q2);
      A[7]=2*(q2*q3-q1*q0);
      A[8]=q0*q0-q1*q1-q2*q2+q3*q3;  
    }

    //here's the integration.
    while (ion->Ek() > ek_cutoff && cfg.t<finish && ion->R.z <= z_cutoff){
      if (tso) cfg.dtOptimize();
      if (cfg.u!=0) cfg.OFirstVV();
      cfg.TimeStepInit();    
      cfg.ReNeighbor();    
      cfg.ForceEval(0);
      cfg.SecondVV();
      cfg.t+=cfg.dt;
    }
    cfg.DelIon(1);

    svector Trans;
    for (I=cfg.begin; I<cfg.end; I++)
      if(!I->is_fixed)
	if(I->R.z > -1*cfg.Lz2){
	  Trans=I->P-Impact;
	  Trans.set(A[0]*Trans.x + A[1]*Trans.y + A[2]*Trans.z,
		    A[3]*Trans.x + A[4]*Trans.y + A[5]*Trans.z,
		    A[6]*Trans.x + A[7]*Trans.y + A[8]*Trans.z);
	  Trans.Rec2Cyl(); //x is now r
	  i=(int)(Trans.x/res);
	  k=(int)((Trans.z-zmin)/res);
	  if (i<Nr && k<Nz){
	    box[i][k].f+=(I->P - I->R).minmag(cfg.Lx, cfg.Ly, cfg.Lz);
	    box[i][k].count++;
	  }
	}
  }
     
  fout.open(outfile.c_str(), ios::out);
  for (i=0; i<Nr; i++)
    for(k=0; k<Nz; k++)
      if (box[i][k].count>0){
	sprintf(out_buffer, "%f\t%f\t%f\n", i*res, k*res+zmin, 
		box[i][k].f);
	fout<<out_buffer;
      }
  fout.close();
  cerr<<endl;
  return(0);
}
