#include "strfun.h"
#include <cstdio>
#include <stdlib.h>
//***************************************************************************
bool sfind (string s1, string s2){
  return (s1.find(s2)!=string::npos);
}
//**************************************************************************
string time2string (double t){
  char strt[12];
  sprintf(strt, "%10.3f", t);
  string b=strt;
  string front=b.substr(0,b.find("."));
  for (int s=0; s<6; s++)
    if (sfind(front," "))
      front.replace(front.find(" "),1,"0");
  return front+"-"+b.substr(7,11);
}
//**************************************************************************
string ftoa (float f, int digits){
  return itoa((int)f, digits);
}
//**************************************************************************
string itoa(int n, int digits){
  n=abs(n);
  string s1="";
  if (n==0) s1="0";
  while (n!=0){
    s1.insert(s1.begin(),n%10+'0');
    n/=10;
  }
  while(s1.size()<digits) s1.insert(s1.begin(),'0');
  return s1;
}
//***************************************************************
string itoa(int n){
  string s1="";
  if (n==0) s1="0";
  while (n!=0){
    s1.insert(s1.begin(),n%10+'0');
    n/=10;
  }
  return s1;
}
//********************************************
string RectifySpecies(string s){
  if (s=="Ar" || s=="AR" || s=="ar") return "Ar";
  else if (s=="H" || s=="h") return "H";
  else if (s=="F" || s=="f") return "F";
  else if (s=="F2" || s=="f2") return "F2";
  else if (s=="CH3" || s=="ch3" || s=="Ch3") return "CH3"; 
  else if (s=="CH2" || s=="ch2" || s=="Ch2")  return "CH2";
  else if (s=="CH" || s=="ch" || s=="Ch") return "CH";
  else if (s=="CF3" || s=="cf3" || s=="Cf3") return "CF3"; 
  else if (s=="CF2" || s=="cf2" || s=="Cf2")  return "CF2";
  else if (s=="CF" || s=="cf" || s=="Cf") return "CF";
  else if (s=="SiF3" || s=="sif3") return "SiF3"; 
  else if (s=="SiF2" || s=="sif2")  return "SiF2";
  else if (s=="SiF" || s=="sif") return "SiF";
  else if (s=="C" || s=="c") return "C";
  else if (s=="Cl" || s=="cl" || s=="CL") return "Cl";
  else if (s=="Cl2" || s=="cl2" || s=="CL2") return "Cl2";
  else{
    cerr<<"Unknown species: "<<s<<endl;
    exit(1);
  }
}


