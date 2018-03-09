//vector.cpp
//function listing for vector class

#include "vector.h"

//***********constructors, destructor*******************
vector::vector(){
  x=0.0; y=0.0; z=0.0;
}

vector::vector(double X, double Y, double Z){
  x=X; y=Y; z=Z;
}

vector::~vector(){}

//********************member functions*****************
void vector::clear(){  //reinit elements
  x=0.0; y=0.0; z=0.0;
}

double vector::sqmag()const{  //find square of magnitude
  double d=x*x+y*y+z*z;
  return d;
}

double vector::mag()const{ //find magnitude
  double d=sqrt(x*x+y*y+z*z);
  return d;
}

void vector::minimg(const double L){ //min image xform
  if (fabs(x/L)>0.5){
    if (x<0.0) x+=L;
    else x-=L;
  }
  if (fabs(y/L)>0.5){
    if (y<0.0) y+=L;
    else y-=L;
  }
  if (fabs(z/L)>0.5){
    if (z<0.0) z+=L;
    else z-=L;
  }
}

void vector::minimg(double Lx, double Ly, double Lz){ //min image xform in 3-D
  if (fabs(x/Lx)>0.5){
    if (x<0) x+=Lx;
    else x-=Lx;
  }
  if (fabs(y/Ly)>0.5){
    if (y<0) y+=Ly;
    else y-=Ly;
  }
  if (fabs(z/Lz)>0.5){
    if (z<0) z+=Lz;
    else z-=Lz;
  }
}

void vector::set(double X, double Y, double Z){
  x=X; y=Y; z=Z;
}

double vector::X()const{
  return x;
}


double vector::Y()const{
  return y;
}

double vector::Z()const{
  return z;
}

//************arithmetic operators********************
void vector::operator = (const vector& v){
  x=v.x; y=v.y; z=v.z;
}

void vector::operator *= (const double d){
  x*=d; y*=d; z*=d;
}

void vector::operator /= (const double d){
  x/=d; y/=d; z/=d;
}

void vector::operator += (const vector& v){
  x+=v.x; y+=v.y; z+=v.z;
}

void vector::operator -= (const vector& v){
  x-=v.x; y-=v.y; z-=v.z;
}

void vector::operator += (const double d){
  x+=d; y+=d; z+=d;
}

void vector::operator -= (const double d){
  x-=d; y-=d; z-=d;
}

vector operator * (const vector& v, const double d){
  vector nv(v.x*d, v.y*d, v.z*d);
  return nv;
}

vector operator * (const double d, const vector& v){
  vector nv(v.x*d, v.y*d, v.z*d);
  return nv;
}

vector operator / (const vector& v, const double d){
  vector nv(v.x/d, v.y/d, v.z/d);
  return nv;
}

double operator ^ (const vector& v1, const vector& v2){
  double d=v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
  return d;
}

vector operator * (const vector& v1, const vector& v2){
  double X=v1.y*v2.z-v1.z*v2.y;
  double Y=v1.z*v2.x-v1.x*v2.z;
  double Z=v1.x*v2.y-v1.y*v2.x;
  vector nv(X, Y, Z);
  return nv;
}

vector operator + (const vector& v1, const vector& v2){
  vector nv(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
  return nv;
}

vector operator - (const vector& v1, const vector& v2){
  vector nv(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
  return nv;
}
//****************input-output operators**************
ostream& operator << (ostream& out, const vector& v){
  out<<"<"<<v.x<<" "<<v.y<<" "<<v.z<<">";
  return out;
}

istream& operator >> (istream& in, vector& v){
  char c1, c2;
  in>>c1>>v.x>>v.y>>v.z>>c2;
  return in;
}
