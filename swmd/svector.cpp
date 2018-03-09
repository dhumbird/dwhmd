//svector.cpp
//function listing for svector class

#include "svector.h"

char c1;

void svector::EulerTrans(double phi, double theta, double psi){
   double q0, q1, q2, q3;
   
   q0=cos(0.5*theta)*cos(0.5*(phi+psi));
   q1=sin(0.5*theta)*cos(0.5*(phi-psi));
   q2=sin(0.5*theta)*sin(0.5*(phi-psi));
   q3=cos(0.5*theta)*sin(0.5*(phi+psi));
   
   double dx=(q0*q0+q1*q1-q2*q2-q3*q3)*x
     +(2*(q1*q2+q0*q3))*y
     +(2*(q1*q3-q0*q2))*z;
   double dy=(2*(q1*q2-q0*q3))*x
     +(q0*q0-q1*q1+q2*q2-q3*q3)*y
     +(2*(q2*q3+q0*q1))*z;
   double dz=(2*(q1*q3+q0*q2))*x
     +(2*(q2*q3-q1*q0))*y
     +(q0*q0-q1*q1-q2*q2+q3*q3)*z;
   
   x=dx; y=dy; z=dz;
}

void svector::Rec2Cyl(){
  double d=sqrt(x*x+y*y);
  y=atan(y/x);
  x=d; //note that x is now r and y is theta (!!!RADIANS!!!)
}

void svector::Cyl2Rec(){
  double d=x*cos(y);
  y=x*sin(y);
  x=d;
}

void svector::Rec2Sph(){
  //rho, theta, phi
  double d=mag();
  z=acos(z/d);
  y=atan(y/x);
  x=d;
}

void svector::Sph2Rec(){
  double d=x;
  x=d*sin(z)*cos(y);
  y=d*sin(z)*sin(y);
  z=d*cos(z);
}
//****************input-output functions**************
void Print(svector v, ostream& out){
  out<<v.x<<" "<<v.y<<" "<<v.z;
}

ostream& operator << (ostream& out, const svector& v){
  out<<"<"<<v.x<<" "<<v.y<<" "<<v.z<<">";
  return out;
}

istream& operator >> (istream& in, svector& v){
  in>>c1>>v.x>>v.y>>v.z>>c1;
  return in;
}

