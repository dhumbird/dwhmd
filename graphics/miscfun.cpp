#include "miscfun.h"
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
//**********************************************************************
bool FileExists(string file){
  ifstream ifs(file.c_str(), ios::in);
  if (ifs) return 1;
  else return 0;
}
//************************************************************************
void CmdError(const char * arg){
  cerr<<"error. "<<arg;
  cerr<<" is not recognized."<<endl;
  exit(0);
}
