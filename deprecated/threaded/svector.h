//svector.h
//class listing for svector class

#ifndef SVECTOR_H
#define SVECTOR_H 1

#include <iostream>
#include <cmath>

class svector{
 public:
  double x;
  double y;
  double z;

  svector(){x=0; y=0; z=0;} //generic constructor
  svector(double X, double Y, double Z) {x=X; y=Y; z=Z;}
  void clear(){x=0; y=0; z=0;}
  bool is_zero(){return (x==0.0 && y==0.0 && z==0.0);}
  double sqmag() {return x*x+y*y+z*z;}
  double mag() {return sqrt(x*x+y*y+z*z);}
  double minsqmag(double&, double&, double&);
  double minmag(double Lx, double Ly, double Lz)
    {return sqrt(minsqmag(Lx, Ly, Lz));}
  void minimg(double&, double&, double&);
  void set(double X, double Y, double Z) {x=X; y=Y; z=Z;}
  void operator = (const svector& v) {x=v.x; y=v.y; z=v.z;}
  void operator *= (double f) {x*=f; y*=f; z*=f;}
  void operator /= (double f) {x/=f; y/=f; z/=f;}
  void operator += (const svector& v) {x+=v.x; y+=v.y; z+=v.z;}
  void operator -= (const svector& v) {x-=v.x; y-=v.y; z-=v.z;}
  void operator += (double f) {x+=f; y+=f; z+=f;}
  void operator -= (double f) {x-=f; y-=f; z-=f;}
  friend double operator ^ (const svector& V1, const svector& V2) //dot
    {return V1.x*V2.x+V1.y*V2.y+V1.z*V2.z;}
  friend svector operator * (const svector& V1, const svector& V2)  //cross
    {return svector(V1.y*V2.z-V1.z*V2.y,
		    V1.z*V2.x-V1.x*V2.z,
		    V1.x*V2.y-V1.y*V2.x);}
  friend svector operator * (const svector& V1, const double f)
    {return svector(V1.x*f, V1.y*f, V1.z*f);}
  friend svector operator * (const double f, const svector& V1)
    {return svector(V1.x*f, V1.y*f, V1.z*f);}
  friend svector operator / (const svector& V1, const double f)
    {return svector(V1.x/f, V1.y/f, V1.z/f);}
  friend svector operator + (const svector& V1, const svector& V2)
    {return svector(V1.x+V2.x, V1.y+V2.y, V1.z+V2.z);}
  friend svector operator - (const svector& V1, const svector& V2)
    {return svector(V1.x-V2.x, V1.y-V2.y, V1.z-V2.z);}
  void EulerTrans(double, double, double);
  void Rec2Cyl();
  void Cyl2Rec();
  void Rec2Sph();
  void Sph2Rec();
  friend void Print(svector, ostream&);
  friend ostream& operator << (ostream&, const svector&); //output
  friend istream& operator >> (istream&, svector&); //input
};

#endif
