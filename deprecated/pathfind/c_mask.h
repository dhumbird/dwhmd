#ifndef MASK_H
#define MASK_H

#include "svector.h"

class mask{
 public:
  svector start, finish;
  short dim;
  mask(){start.clear(); finish.clear();}
  bool InsideMask(const svector& v){
    switch (dim){
    case 0:
      return (v.x > start.x && v.x < finish.x); break;
    case 1:
      return (v.y > start.y && v.y < finish.y); break;
    }
  } 
  friend ostream& operator << (ostream& out, const mask& m){
    switch (m.dim){
    case 0:
      out<<"x"<<" "<<m.start.x<<" "<<m.finish.x-m.start.x; break;
    case 1: 
      out<<"y"<<" "<<m.start.y<<" "<<m.finish.y-m.start.y; break;
    }
    return out;
  }
  friend istream& operator >> (istream& in, mask& m){
    char c1; double d1, d2;
    in>>c1>>d1>>d2;
    if (c1=='x'){
      m.dim=0;
      m.start.x=d1;
      m.finish.x=m.start.x+d2;
    }
    else if (c1=='y'){
      m.dim=1;
      m.start.y=d1;
      m.finish.y=m.start.y+d2;
    }
    return in;
  }
};


#endif
