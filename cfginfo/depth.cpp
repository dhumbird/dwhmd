#include "i_config.h"

extern map<short,double> F_MIN, F_MAX;

//***************************************************************************
void InfoDepth(int argc, char * argv[]){
  string outfile;
  int mono=1;
  float every=-1;
  vector<string>filelist;
  string arg, orig, scfg;
  bool _stdout=0;
  svector R; double r,f,f1;
  float res=2.5;
  for (int q=2; q<argc; q++){
    arg=argv[q];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-mono") mono=atoi(argv[++q]);
    else if (arg=="-every") every=atof(argv[++q]);
    else if (arg=="-stdout") _stdout=1;
    else if (arg=="-res") res=atof(argv[++q]);
  }
  config cfg(scfg);
  
  float top=-1000;
  float bottom=1000;
 
  for (VAI i=cfg.begin; i<cfg.end; i++){
    if (i->R.z > top) top=i->R.z;
    if (i->R.z < bottom) bottom=i->R.z;
  }
  
  int Nz=(int)((top-bottom)/res)+1;
  int zdensC[Nz];
  int zdensF[Nz];
  int zdensSi[Nz];
  int zdensC_C[Nz];
  int zdensC_F[Nz];
  int zdensC_Si[Nz];
  int zdensSi_F[Nz];

  res=(top-bottom)/(float)Nz;

  for (int x=0; x<Nz; x++){
    zdensC[x]=zdensF[x]=zdensSi[x]=0;
    zdensC_C[x]=zdensC_F[x]=zdensC_Si[x]=zdensSi_F[x]=0;
  }
  int z, type;
  for (VAI i=cfg.begin; i<cfg.end; i++){
    z=(int)((i->R.z-bottom)/res);
    if (z>=0 && z<Nz){
      if (i->id==6) zdensC[z]++;
      if (i->id==9) zdensF[z]++; 
      if (i->id==14) zdensSi[z]++; 
    }
    for (SAI j=i->nlist.begin(); j!=i->nlist.end(); j++){
      if (i->ix > (*j)->ix){
	type=i->id+(*j)->id;
	r=(i->R.z + (*j)->R.z)/2;
	z=(int)((r-bottom)/res);
	if (z>=0 && z<Nz){
	  if (type==12) zdensC_C[z]++;
	  if (type==15) zdensC_F[z]++;
	  if (type==20) zdensC_Si[z]++;
	  if (type==23) zdensSi_F[z]++;
	}
      }
    }
  }
  char outbuf[256];
  cout<<"     z     C     F    Si   C-C   C-F  C-Si  Si-F\n";
  printf("%6.2f  %4i  %4i  %4i  %4i  %4i  %4i  %4i\n",0.0,0,0,0,0,0,0,0);
  for (z=Nz-1; z>=0; z--){
    sprintf(outbuf,
	    "%6.2f  %4i  %4i  %4i  %4i  %4i  %4i  %4i\n",
	    (top-bottom)-z*res, zdensC[z], zdensF[z], zdensSi[z],
	    zdensC_C[z], zdensC_F[z], zdensC_Si[z], zdensSi_F[z]);
    cout<<outbuf;
  }
}
//***************************************************************
// void InfoSurf(int argc, char * argv[]){
//   string outfile;
//   string arg, orig, scfg;
//   bool _stdout=0;
//   svector R; double r,f,f1;
//   float res=.25;
//   for (int q=2; q<argc; q++){
//     arg=argv[q];
//     if (sfind(arg, ".cfg")) scfg=arg;
//     else if (arg=="-stdout") _stdout=1;
//     else if (arg=="-res") res=atof(argv[++q]);
//   }
//   config cfg(scfg);
  
//   float top=-1000;
//   float bottom=1000;
 
//   for (VAI i=cfg.begin; i<cfg.end; i++){
//     if (i->R.z > top) top=i->R.z;
//     if (i->R.z < bottom) bottom=i->R.z;
//   }
  
//   int Nz=(int)((top-bottom)/res)+1;
//   int zdensC[Nz];
//   int zdensF[Nz];
//   int zdensSi[Nz];
//   int zdensC_C[Nz];
//   int zdensC_F[Nz];
//   int zdensC_Si[Nz];
//   int zdensSi_F[Nz];

//   res=(top-bottom)/(float)Nz;

//   for (int x=0; x<Nz; x++){
//     zdensC[x]=zdensF[x]=zdensSi[x]=0;
//     zdensC_C[x]=zdensC_F[x]=zdensC_Si[x]=zdensSi_F[x]=0;
//   }
//   int z, type;
//   for (VAI i=cfg.begin; i<cfg.end; i++){
//     z=(int)((i->R.z-bottom)/res);
//     if (z>=0 && z<Nz){
//       if (i->id==6) zdensC[z]++;
//       if (i->id==9) zdensF[z]++; 
//       if (i->id==14) zdensSi[z]++; 
//     }
//     for (SAI j=i->nlist.begin(); j!=i->nlist.end(); j++){
//       if (i->ix > (*j)->ix){
// 	type=i->id+(*j)->id;
// 	r=(i->R.z + (*j)->R.z)/2;
// 	z=(int)((r-bottom)/res);
// 	if (z>=0 && z<Nz){
// 	  if (type==12) zdensC_C[z]++;
// 	  if (type==15) zdensC_F[z]++;
// 	  if (type==20) zdensC_Si[z]++;
// 	  if (type==23) zdensSi_F[z]++;
// 	}
//       }
//     }
//   }
//   char outbuf[256];
//   cout<<"     z     C     F    Si   C-C   C-F  C-Si  Si-F\n";
//   printf("%6.2f  %4i  %4i  %4i  %4i  %4i  %4i  %4i\n",0.0,0,0,0,0,0,0,0);
//   for (z=Nz-1; z>=0; z--){
//     sprintf(outbuf,
// 	    "%6.2f  %4i  %4i  %4i  %4i  %4i  %4i  %4i\n",
// 	    (top-bottom)-z*res, zdensC[z], zdensF[z], zdensSi[z],
// 	    zdensC_C[z], zdensC_F[z], zdensC_Si[z], zdensSi_F[z]);
//     //if (!_stdout) fout<<outbuf;
//     cout<<outbuf;
//   }
// }
