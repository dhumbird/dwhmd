//svector.h
//class listing for svector class

#ifndef SVECTOR_H
#define SVECTOR_H 1

#include <iostream>
#include <cmath>

class svector{
 public:
  double x,y,z;

  svector() {x=y=z=0;}
  svector(double X, double Y, double Z) {x=X; y=Y; z=Z;}
  void clear() {x=y=z=0;}
  bool is_zero() {return (x==0.0 && y==0.0 && z==0.0);}
  double sqmag() {return x*x+y*y+z*z;}
  double mag() {return sqrt(x*x+y*y+z*z);}
  double minsqmag(double Lx, double Ly, double Lz)
    {minimg(Lx, Ly, Lz); return sqmag();}
  double minmag(double Lx, double Ly, double Lz)
    {return sqrt(minsqmag(Lx, Ly, Lz));}
  void minimg(double Lx, double Ly, double Lz){
    if (abs(x/Lx)>0.5) x = (x<0) ? x+Lx : x-Lx;
    if (abs(y/Ly)>0.5) y = (y<0) ? y+Ly : y-Ly;
    if (abs(z/Lz)>0.5) z = (z<0) ? z+Lz : z-Lz;
  }
  void set(double X, double Y, double Z) {x=X; y=Y; z=Z;}
  void operator = (const svector& v) {x=v.x; y=v.y; z=v.z;}
  void operator *= (double f) {x*=f; y*=f; z*=f;}
  void operator /= (double f) {x/=f; y/=f; z/=f;}
  void operator += (double f) {x+=f; y+=f; z+=f;}
  void operator -= (double f) {x-=f; y-=f; z-=f;}
  void operator += (const svector& v) {x+=v.x; y+=v.y; z+=v.z;}
  void operator -= (const svector& v) {x-=v.x; y-=v.y; z-=v.z;}
  friend double operator ^ (const svector& V1, const svector& V2) //dot
    {return V1.x*V2.x+V1.y*V2.y+V1.z*V2.z;}
  friend svector operator - (const svector& V1, const svector& V2)
    {return svector(V1.x-V2.x, V1.y-V2.y, V1.z-V2.z);}
  friend svector operator * (const svector& V1, const svector& V2)  //cross
    {return svector(V1.y*V2.z-V1.z*V2.y, 
		    V1.z*V2.x-V1.x*V2.z,
		    V1.x*V2.y-V1.y*V2.x);}
};

#endif
