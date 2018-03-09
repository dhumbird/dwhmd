#include "i_config.h"

//*********************************************************************
// void EkMaximum(int argc, char * argv[]){
//   int j; int TotalSiAtoms=0;
//   string infile;
//   double ek=0; double TotalEk=0; double AvgEk=0;
//   double EkMin=1000; double EkMax=0;
//   vector<particle>::iterator i;

//   for (j=2; j<argc; j++){
//     infile=argv[j];
//     {
//       config cfg(infile);
//       for (i=cfg.begin; i<cfg.end; i++)
// 	if (i->id==14)
// 	  if(!i->is_fixed){
// 	    TotalSiAtoms++;
// 	    ek=i->Ek();
// 	    TotalEk+=ek;
// 	    if(ek<EkMin) EkMin=ek;
// 	    if(ek>EkMax) EkMax=ek;
// 	  }
//     }
//   }
//   AvgEk=TotalEk/TotalSiAtoms;
//   cout<<EkMax<<" max, "<<EkMin<<" min, "<<AvgEk<<" average"<<endl;
// }
//***********************************************************************
void EkDistribution(int argc, char * argv[]){
  int bin=0; int j;
  string infile=argv[2];
  config cfg(infile);
  double ek; double EkMin=1000; double EkMax=0;
  VAI i;

  for (i=cfg.begin; i<cfg.end; i++)
    if(!i->is_fixed){
      ek=i->Ek();
      if(ek<EkMin) EkMin=ek;
      if(ek>EkMax) EkMax=ek;
    }
  
  double dr=EkMax/100;
  int num_bins=(int)(EkMax/dr)+1;
  int hist[num_bins];
  for (j=0; j<num_bins; j++)
    hist[j]=0;
  for (i=cfg.begin; i<cfg.end; i++)
    if(!i->is_fixed){
      bin=(int)(i->Ek()/dr);
      if (bin<num_bins) hist[bin]+=1;
    }
  for (j=0; j<num_bins; j++)
    cout<<(float)j*dr<<"\t"<<(float)j*dr*11600<<"\t"<<hist[j]<<endl;
}
//***********************************************************************
// void EkZ(int argc, char * argv[]){
//   string infile=argv[2];
//   config cfg(infile);
//   cfg.Rz_sort_asc();
//   for (VPI i=cfg.begin; i<cfg.end; i++)
//     cout<<i->ix<<"\t"<<i->Ek()<<endl;
// }
//**********************************************************************
  

