#ifndef SPHERE_H
#define SPHERE_H

#include "../svector.h"
#include "rgb.h"

class sphere{
 public:
  svector R;
  float radius;
  rgb color;
  float alpha;
  sphere(){radius=0; alpha=0;}
  sphere(svector R1){R=R1;}
  void Color(float r, float g, float b){color.set(r,g,b);}
  void Transform(float phi, float theta, float psi,
		 svector Offset, float scale){
    R.EulerTrans(phi, theta, psi);
    R-=Offset;
    R/=scale;
    radius/=scale;
  }
  void cout_r3d(){
    if (alpha>0){
      cout<<"8"<<endl;
      cout<<"20.0 -1.0   1. 1. 1. "<<alpha<<"  0 0 0 0"<<endl;
    }
    cout<<"2\n";
    Print(R,cout);
    cout<<" "<<radius<<" "<<color<<endl;
    if (alpha>0) cout<<"9\n";
  }
};

#endif
