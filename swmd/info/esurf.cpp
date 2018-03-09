#include "ndens.h"

void NDensity(int argc, char * argv[]){
  int i, j, k;
  float res=0.5;
  int mc_tries=1000;
  string arg; string scfg;
  bool rec=0;
  for (i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-rec") rec=1;
    else if (arg=="-res") res=atof(argv[++i]);
    else if (arg=="-mc") mc_tries=atoi(argv[++i]);
    else CmdError(arg.c_str());
  }

  config cfg(scfg);
  double Lx=abs(cfg.Lx);
  double Ly=abs(cfg.Ly);
  double Lz=abs(cfg.Lz);
  double Lx2=Lx/2;
  double Ly2=Ly/2;
  double Lz2=Lz/2;

  float aSi=5.434;
  float aSi2=aSi/2;
  float aSisq=aSi2*aSi2;

  int seed=(unsigned)clock()+(unsigned)time(NULL);
  srand(seed);

  float adder=1/(float)mc_tries;
  
  if (!rec){
    int Nr=(int)((Lx2 >? Ly2)/res)+1;
    int Nz=(int)(Lz/res)+1;
    float box[Nr][Nz];
    double vol=0;

    for (i=0; i<Nr; i++) for (k=0; k<Nz; k++) box[i][k]=0;
    
    svector R;
    int atom=cfg.N;
    int a=0;
    vector<particle>::iterator I;
    double x=0; double y=0; double z=0;
    for (I=cfg.begin; I<cfg.end; I++){
      cerr<<"\r"<<atom--;
      if(!I->is_fixed)
	if(I->id==cfg.mat)
	  for (a=0; a < mc_tries; a++){
	    x=aSi2; y=aSi2; z=aSi2;
	    while(x*x+y*y+z*z > aSisq){
	      x=rand01()*aSi-aSi2;
	      y=rand01()*aSi-aSi2;
	      z=rand01()*aSi-aSi2;
	    }
	    R.set(x,y,z);
	    R+=I->R;
	    R.minimg(cfg.Lx, cfg.Ly, cfg.Lz);
	    R.Rec2Cyl();
	    i=(int)(R.x/res);
	    k=(int)((fabs(cfg.Lz/2)-R.z)/res);
	    if (i<Nr && k<Nz) box[i][k]+=adder;
	  }
    }
    cerr<<endl;
  
    string file=cfg.name+"_"+time2string(cfg.t)+".Nrz";
    ofstream fout(file.c_str(), ios::out);
    char outbuf[256];
    for (i=0; i<Nr; i++){
      vol=PI*res*res*res*(pow(i+1,2)-pow(i,2));
      for (k=0; k<Nz; k++)
	if(box[i][k]>0){
	  sprintf(outbuf, "%f\t%f\t%f\n", i*res, Lz2-k*res, box[i][k]/vol);
	  fout<<outbuf;
	}
    }
  }
  
  else{
    int Nx=(int)(Lx/res)+1;
    int Ny=(int)(Ly/res)+1;
    int Nz=(int)(Lz/res)+1;
    float box[Nx][Ny][Nz];
    double vol=pow(res,3);

    for (i=0; i<Nx; i++) for(j=0; j<Ny; j++) for (k=0; k<Nz; k++)
      box[i][j][k]=0;
    
    svector R;
    int atom=cfg.N;
    int a=0;
    vector<particle>::iterator I;
    double x=0; double y=0; double z=0;
    for (I=cfg.begin; I<cfg.end; I++){
      cerr<<"\r"<<atom--;
      if(!I->is_fixed)
	if(I->id==cfg.mat)
	  for (a=0; a < mc_tries; a++){
	    x=aSi2; y=aSi2; z=aSi2;
	    while(x*x+y*y+z*z > aSi2){
	      x=rand01()*aSi-aSi2;
	      y=rand01()*aSi-aSi2;
	      z=rand01()*aSi-aSi2;
	    }
	    R.set(x,y,z);
	    R+=I->R;
	    R.minimg(cfg.Lx, cfg.Ly, cfg.Lz);
	    i=(int)((Lx2-R.x)/res);
	    j=(int)((Ly2-R.y)/res);
	    k=(int)((Lz2-R.z)/res);
	    if (i<Nx && j<Ny && k<Nz) box[i][j][k]+=adder;
	  }
    }
    cerr<<endl;
  
    string file=cfg.name+"_"+time2string(cfg.t)+".Nxyz";
    ofstream fout(file.c_str(), ios::out);
    char outbuf[256];
    for (i=0; i<Nx; i++)
      for (j=0; j<Ny; j++)
	for (k=0; k<Nz; k++)
	  if(box[i][j][k]>0){
	    sprintf(outbuf, "%f\t%f\t%f\t%f\n", Lx2-i*res, Ly2-j*res,
		    Lz2-k*res, box[i][j][k]/vol);
	    fout<<outbuf;
	  }
  }
}

