#ifndef PLANE_H
#define PLANE_H

#include "../svector.h"
#include "rgb.h"
#include <map>
#include "r3d_util.h"

class plane{
 public:
  svector v1, v2, v3, v4;
  rgb c;
  float alpha;
  void Color(float r, float g, float b){c.set(r,g,b);}
  plane(){alpha=0;}
  plane(svector r1, svector r2, svector r3, svector r4, rgb r)
    {v1=r1; v2=r2; v3=r3; v4=r4; c=r;}
  void Transform(float phi, float theta, float psi,
		 svector Offset, float scale){
    v1.EulerTrans(phi, theta, psi);
    v2.EulerTrans(phi, theta, psi);
    v3.EulerTrans(phi, theta, psi);
    v4.EulerTrans(phi, theta, psi);
    v1-=Offset;
    v2-=Offset;
    v3-=Offset;
    v4-=Offset;
    v1/=scale;
    v2/=scale;
    v3/=scale;
    v4/=scale;
  }
  void cout_r3d(){
    if (alpha > 0){
      cout<<"8"<<endl;
      cout<<"20.0 -1.0   1. 1. 1. "<<alpha<<"  0 0 0 0"<<endl;
    }
    map<float, short> rmap;
    rmap[(v1-v2).sqmag()]=12;
    rmap[(v1-v3).sqmag()]=13;
    rmap[(v1-v4).sqmag()]=14;
    rmap[(v2-v3).sqmag()]=23;
    rmap[(v2-v4).sqmag()]=24;
    rmap[(v3-v4).sqmag()]=34;
    map<float, short>::iterator m=rmap.end(); m--;
    switch (m->second){
    case 12:
      cout<<"1"<<endl; PrintV(v1); PrintV(v2); PrintV(v3); cout<<c<<endl;
      cout<<"1"<<endl; PrintV(v1); PrintV(v2); PrintV(v4); cout<<c<<endl;
      break;
    case 13:
      cout<<"1"<<endl; PrintV(v1); PrintV(v3); PrintV(v4); cout<<c<<endl;
      cout<<"1"<<endl; PrintV(v1); PrintV(v3); PrintV(v2); cout<<c<<endl;
      break;
    case 14:
      cout<<"1"<<endl; PrintV(v1); PrintV(v4); PrintV(v3); cout<<c<<endl;
      cout<<"1"<<endl; PrintV(v1); PrintV(v4); PrintV(v2); cout<<c<<endl;
      break;
    case 23:
      cout<<"1"<<endl; PrintV(v2); PrintV(v3); PrintV(v1); cout<<c<<endl;
      cout<<"1"<<endl; PrintV(v2); PrintV(v3); PrintV(v4); cout<<c<<endl;
      break;
    case 24:
      cout<<"1"<<endl; PrintV(v2); PrintV(v4); PrintV(v1); cout<<c<<endl;
      cout<<"1"<<endl; PrintV(v2); PrintV(v4); PrintV(v3); cout<<c<<endl;
      break;
    case 34:
      cout<<"1"<<endl; PrintV(v3); PrintV(v4); PrintV(v1); cout<<c<<endl;
      cout<<"1"<<endl; PrintV(v3); PrintV(v4); PrintV(v2); cout<<c<<endl;
      break;
    }
    if (alpha > 0) cout<<"9\n";
  }
};

#endif
