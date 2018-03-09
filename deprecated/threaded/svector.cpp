//svector.cpp
//function listing for svector class

#include "svector.h"

char c1;

//**************************************************************************
double svector::minsqmag(double& Lx, double& Ly, double& Lz){
  minimg(Lx, Ly, Lz);
  return sqmag();
}
//*************************************************************************
void svector::minimg(double& Lx, double& Ly, double& Lz){ //min image xform in 3-D
  //if a length is negative, that is a flag that the box is open in that
  //direction...

  if (Lx>0)
    if (fabs(x/Lx)>0.5){
      if (x<0) x+=Lx;
      else x-=Lx;
    }
  if (Ly>0)
    if (fabs(y/Ly)>0.5){
      if (y<0) y+=Ly;
      else y-=Ly;
    }
  if(Lz>0)
    if (fabs(z/Lz)>0.5){
      if (z<0) z+=Lz;
      else z-=Lz;
    }
}

//**************coordinate transformations***************
void svector::EulerTrans(double phi, double theta, double psi){
   double q0, q1, q2, q3;
   double A[3][3];

   q0=cos(0.5*theta)*cos(0.5*(phi+psi));
   q1=sin(0.5*theta)*cos(0.5*(phi-psi));
   q2=sin(0.5*theta)*sin(0.5*(phi-psi));
   q3=cos(0.5*theta)*sin(0.5*(phi+psi));
   
   A[0][0]=q0*q0+q1*q1-q2*q2-q3*q3;
   A[0][1]=2*(q1*q2+q0*q3);
   A[0][2]=2*(q1*q3-q0*q2);
   A[1][0]=2*(q1*q2-q0*q3);
   A[1][1]=q0*q0-q1*q1+q2*q2-q3*q3;
   A[1][2]=2*(q2*q3+q0*q1);
   A[2][0]=2*(q1*q3+q0*q2);
   A[2][1]=2*(q2*q3-q1*q0);
   A[2][2]=q0*q0-q1*q1-q2*q2+q3*q3;
   
   double dx=A[0][0]*x+A[0][1]*y+A[0][2]*z;
   double dy=A[1][0]*x+A[1][1]*y+A[1][2]*z;
   double dz=A[2][0]*x+A[2][1]*y+A[2][2]*z;
   
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

