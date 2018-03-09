#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "svector.h"
#include "rgb.h"

class triangle{
 public:
  svector v1, v2, v3;
  rgb c;
  float alpha;
  void Color(float r, float g, float b){c.set(r,g,b);}
  triangle(){alpha=0;}
  triangle(svector r1, svector r2, svector r3, rgb r)
    {v1=r1; v2=r2; v3=r3; c=r; alpha=0;}
  void Transform(float phi, float theta, float psi,
		 svector Offset, float scale){
    v1.EulerTrans(phi, theta, psi);
    v2.EulerTrans(phi, theta, psi);
    v3.EulerTrans(phi, theta, psi);
    v1-=Offset;
    v2-=Offset;
    v3-=Offset;
    v1/=scale;
    v2/=scale;
    v3/=scale;
  }
  void cout_r3d(){
    if (alpha > 0){
      cout<<"8"<<endl;
      cout<<"20.0 -1.0   1. 1. 1. "<<alpha<<"  0 0 0 0"<<endl;
    }
    cout<<"1\n";
    v1.Print();
    v2.Print();
    v3.Print(); cout<<c<<endl;
    if (alpha > 0) cout<<"9\n";
  }
};

#endif
