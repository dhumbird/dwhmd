#include "strfun.h"

//***************************************************************************
bool sfind (string s1, string s2){
  return (s1.find(s2)!=string::npos);
}
//**************************************************************************
string time2string (double t){
  return itoa((int)(t*1000), 6);
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
//**************************************************************************
pair<int, int> CheckIntRange(string s){
  string::size_type div_loc;
  int begin; int end;
  pair<int, int> p;

  div_loc=s.find(":");
  if (div_loc==string::npos){ 
    cerr<<"Illegal range separator.\n";
    exit(1);
  }
  else{
    s.erase(s.find("["), 1);
    s.erase(s.find("]"), 1);
    begin=atoi((s.substr(0, div_loc-1)).c_str());
    end=atoi((s.substr(div_loc, string::npos)).c_str());
    if (begin < end){
      p=make_pair(begin, end);
      return p;
    }
    else{
      cerr<<"Illegal range.\n";
      exit(1);
    }
  } 
}
//**************************************************************************
pair<float, float> CheckFloatRange(string s){
  string::size_type div_loc;
  float begin; float end;
  pair<float, float> p;

  div_loc=s.find(":");
  if (div_loc==string::npos){ 
    cerr<<"Illegal range separator.\n";
    exit(1);
  }
  else{
    s.erase(s.find("["), 1);
    s.erase(s.find("]"), 1);
    begin=atof((s.substr(0, div_loc-1)).c_str());
    end=atof((s.substr(div_loc, string::npos)).c_str());
    if (begin < end){
      p=make_pair(begin, end);
      return p;
    }
    else{
      cerr<<"Illegal range.\n";
      exit(1);
    }
  } 
}
