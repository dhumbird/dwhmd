#include "../config.h"
#include "../errors.h"

//***********************************************************************
void ZDistribution(int argc, char * argv[]){
  int bin=0; int j;
  double r; double rMin=1000; double rMax=-1000;
  vector<particle>::iterator i;
  string arg; string scfg;
  float res=0.05;
  for (j=2; j<argc; j++){
    arg=argv[j];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-res") res=atof(argv[++j]);
    else CmdError(arg.c_str());
  }
  config cfg(scfg);

  //find atom with highest z 
  for (i=cfg.begin; i<cfg.end; i++){
    if(i->R.z > neg(cfg.Lz/2)){
      r=i->R.z;
      rMin = r <? rMin;
      rMax = r >? rMax;
    }
  }
  int num_bins=(int)((rMax-rMin)/res)+1;
  int hist[num_bins];
  for (j=0; j<num_bins; j++) hist[j]=0;
  
  for (i=cfg.begin; i<cfg.end; i++){
    bin=(int)((rMax-i->R.z)/res);
    if (bin<num_bins) hist[bin]+=1;
  }
  for (j=0; j<num_bins; j++){
    cout<<(float)j*res<<"\t"<<hist[j]<<endl;
  }
}
//**************************************************************************
void RDistribution(int argc, char * argv[]){
  int bin=0; int j;
  double r; double rMin=1000; double rMax=-1000;
  vector<particle>::iterator i;
  string arg; string scfg;
  float zg=-1000;
  float res=0.05;
  for (j=2; j<argc; j++){
    arg=argv[j];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-res") res=atof(argv[++j]);
    else if (arg=="-zg") zg=atof(argv[++j]);
    else CmdError(arg.c_str());
  }
  config cfg(scfg);

  int num_bins=(int)((fabs(cfg.Lx/2) >? fabs(cfg.Ly/2))/res)+1;
  if (zg==-1000) zg=neg(cfg.Lz/2);
  int hist[num_bins];
  for (j=0; j<num_bins; j++) hist[j]=0;
  svector R;
  for (i=cfg.begin; i<cfg.end; i++)
    if (i->R.z > zg){
      R=i->R;
      R.Rec2Cyl();
      bin=(int)(R.x/res);
      if (bin<num_bins) hist[bin]+=1;
    }
  for (j=0; j<num_bins; j++){
    cout<<(float)j*res<<"\t"<<hist[j]<<endl;
  }
}
//*************************************************************************

