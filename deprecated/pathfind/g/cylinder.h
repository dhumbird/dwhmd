#ifndef CYL_H
#define CYL_H

#include "../svector.h"
#include "rgb.h"

class cylinder{
 public:
  svector R1, R2;
  float radius;
  rgb color;
  float alpha;
  cylinder(){radius=0; alpha=0;}
  cylinder(svector f, svector g){R1=f; R2=g;}
  void Color(float r, float g, float b){color.set(r,g,b);}
  void Transform(float phi, float theta, float psi,
		 svector Offset, float scale){
    R1.EulerTrans(phi, theta, psi);
    R2.EulerTrans(phi, theta, psi);
    R1-=Offset;
    R2-=Offset;
    R1/=scale;
    R2/=scale;
    radius/=scale;
  }
  void cout_r3d(){
    if (alpha>0){
      cout<<"8"<<endl;
      cout<<"20.0 -1.0   1. 1. 1. "<<alpha<<"  0 0 0 0"<<endl;
    }
    cout<<"5\n";
    Print(R1,cout);
    cout<<" "<<radius<<" ";
    Print(R2,cout);
    cout<<" 0 "<<color<<endl;
    if (alpha>0) cout<<"9\n";
  }
};

#endif
