//binlist.cpp
//function listing for binlist class

#include "binlist.h"

binlist::binlist(double MIN, double MAX, double INC){
  int i;
  min=MIN;
  max=MAX;
  inc=INC;
  bins=(int)((max-min)/inc)+1;
  histo=new int[bins];
  for (i=0; i<bins; i++) histo[i]=0;
}

binlist::~binlist(){
  delete list;
}

void binlist::toss(double r){
  int i=0;
  double bottom=0;
  while(r>bottom){
    bottom+=inc;
    i++;
  }
  histo[i-1]++;
}

void binlist::display(double c){
  int i;
  for (i=0; i<bins; i++){
    cout<<i*inc<<"\t"<<histo[i]*c<<endl;
  }
}
