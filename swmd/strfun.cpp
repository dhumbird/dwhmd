#include "strfun.h"
#include <cstdio>
//***************************************************************************
bool sfind (string s1, string s2){
  return (s1.find(s2)!=string::npos);
}
//**************************************************************************
string time2string (double t){
  char strt[12];
  sprintf(strt, "%8.3f", t);
  string bitch=strt;
  string front=bitch.substr(0,bitch.find("."));
  for (int s=0; s<4; s++)
    if (sfind(front," "))
      front.replace(front.find(" "),1,"0");
  return front+"-"+bitch.substr(5,9);
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

