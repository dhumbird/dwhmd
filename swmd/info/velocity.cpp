#include "../config.h"

void SpeedDistribution(int argc, char * argv[]){
  int bin=0; int j;
  vector<particle>::iterator i;
  string infile=argv[2];
  config cfg(infile);
  double v; double vMin=1000; double vMax=-1000;

  for (i=cfg.begin; i<cfg.end; i++){
    if (i->id==14){
      if(!i->is_fixed){
	v=i->V.mag();
	if(v<vMin) vMin=v;
	if(v>vMax) vMax=v;
      }
    }
  }
  
  double dr=0.1;
  int num_bins=(int)(vMax/dr)+1;
  int hist[num_bins];
  for (j=0; j<num_bins; j++){
    hist[j]=0;
  }
  for (i=cfg.begin; i<cfg.end; i++){
    if (!i->is_fixed){
      bin=(int)(i->V.mag()/dr);
      if (bin<num_bins) hist[bin]+=1;
    }
  }
  for (j=0; j<num_bins; j++){
    cout<<(float)j*dr<<"\t"<<hist[j]<<endl;
  }
}   
//***********************************************************************
void VelocityDistribution(int argc, char * argv[]){
  string arg=argv[2];
  int bin=0; int j;
  double v; double vMin=1000; double vMax=-1000;
  vector<particle>::iterator i;
  string infile=argv[3];
  config cfg(infile);
  
  for (i=cfg.begin; i<cfg.end; i++){
    if (i->id==14){
      if(!i->is_fixed){
	if (arg=="x") v=i->V.x;
	else if (arg=="y") v=i->V.y;
	else if (arg=="z") v=i->V.z;
	if(v<vMin) vMin=v;
	if(v>vMax) vMax=v;
      }
    }
  }
  double dr=0.1;
  int num_bins=(int)(vMax/dr)+1;
  int hist[num_bins];
  for (j=0; j<num_bins; j++){
    hist[j]=0;
  }
  
  for (i=cfg.begin; i<cfg.end; i++){
    if (!i->is_fixed){
      if (arg=="x") bin=(int)(i->V.x/dr);
      else if (arg=="y") bin=(int)(i->V.y/dr);
      else if (arg=="z") bin=(int)(i->V.z/dr);
      if (bin<num_bins) hist[bin]+=1;
    }
  }
  for (j=0; j<num_bins; j++){
    cout<<(float)j*dr<<"\t"<<hist[j]<<endl;
  }
} 

