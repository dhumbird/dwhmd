//trj2r3d

#include <fstream>
#include "../svector.h"
#include "../errors.h"
#include "../strfun.h"
#include <string>
#include "../consts.h"

int main (int argc, char* argv[]){
  string arg;
  string infile="temp_000000.cfg";
  int i, s;
  double Lscale=0;
  float theta=-90; float psi=0; float phi=0; //Euler angles (deg)
  string shadow="F";
  float eyepos=4.0; 
  bool nowalls=0; bool bg=0;
  float tsres=0.1; float tsrad=4; float tszres=0.1;
  svector Offset(0,0,0);
  bool drawtop=0; bool hideatom=0;
  float boxscale=4; int sizex=20; int sizey=20;
  svector light(1, 1, 1);
  double Lx, Ly, Lz;

  if (argc==1){
    cerr<<"supply some more arguments.\n";
    exit(1);
  }
  
  for (i=1; i<argc; i++){
    arg=argv[i];
    if (sfind(arg,".trj")) infile=arg;
    else if (arg=="-ang"){
      theta=atof(argv[++i]);
      psi=atof(argv[++i]);
      phi=atof(argv[++i]);
    }
    else if (arg=="-nowalls") nowalls=1;
    else if (arg=="-bg") bg=1;
    else if (arg=="-size"){
      sizex=atoi(argv[++i]);
      sizey=atoi(argv[++i]);
    }
    else if (arg=="-shadow") shadow="T";
    else if (arg=="-eyepos") eyepos=atof(argv[++i]);
    else if (arg=="-boxscale") boxscale=atof(argv[++i]);
    else if (arg=="-offset"){
      Offset=svector(atof(argv[++i]), atof(argv[++i]), atof(argv[++i]));
    }
    else if (arg=="-light"){
      light=svector(atof(argv[++i]), atof(argv[++i]), atof(argv[++i]));
    }
    else CmdError(arg.c_str());
  }

  Lx=65; Ly=65; Lz=42;

    
  //pick the dimension to scale in
  Lscale=fabs(Lx);
  if (fabs(Ly)>Lscale) Lscale=fabs(Ly);
  if (fabs(Lz)>Lscale) Lscale=fabs(Lz);
  boxscale*=Lscale;
  cerr<<"boxscale "<<boxscale<<endl;
  //image descriptors***************************************
  float scale=2.5; //times normal size
  double A[3][3]; //rotation matrix
  
  double q0=cos(0.5*theta/180*PI)*cos(0.5*(phi+psi)/180*PI);
  double q1=sin(0.5*theta/180*PI)*cos(0.5*(phi-psi)/180*PI);
  double q2=sin(0.5*theta/180*PI)*sin(0.5*(phi-psi)/180*PI);
  double q3=cos(0.5*theta/180*PI)*sin(0.5*(phi+psi)/180*PI);
  
  A[0][0]=q0*q0+q1*q1-q2*q2-q3*q3;
  A[0][1]=2*(q1*q2+q0*q3);
  A[0][2]=2*(q1*q3-q0*q2);
  A[1][0]=2*(q1*q2-q0*q3);
  A[1][1]=q0*q0-q1*q1+q2*q2-q3*q3;
  A[1][2]=2*(q2*q3+q0*q1);
  A[2][0]=2*(q1*q3+q0*q2);
  A[2][1]=2*(q2*q3-q1*q0);
  A[2][2]=q0*q0-q1*q1-q2*q2+q3*q3;
  
  cout<<"beeyatch.r3d"<<endl;
  cout<<sizex<<" "<<sizey<<"     tiles in x,y \n";
  cout<<"16 16     pixels (x,y) per tile \n";
  cout<<"4         anti-aliasing level 4; 3x3->2x2 \n";
  if (bg) cout<<"1 1 1";
  else cout<<" 0 0 0";
  cout<<"     background color\n";
  cout<<shadow<<"         shadows cast(T/F) \n";
  cout<<"25        Phong power \n";
  cout<<"0.25      secondary light contribution \n";
  cout<<"0.05      ambient light contribution \n";
  cout<<"0.25      specular reflection component \n";
  cout<<eyepos<<"       eye position \n";
  cout<<light.x<<" "<<light.y<<" "<<light.z;
  cout<<"     main light source position \n";
  
  //view matrix
  cout<<A[0][0]<<" "<<A[0][1]<<" "<<A[0][2]<<" 0"<<endl;
  cout<<A[1][0]<<" "<<A[1][1]<<" "<<A[1][2]<<" 0"<<endl;
  cout<<A[2][0]<<" "<<A[2][1]<<" "<<A[2][2]<<" 0"<<endl;
  cout<<"0 0 0 "<<1.0/scale<<endl;
  
  cout<<"3         mixed objects \n";
  cout<<"*        (free format triangle and plane descriptors) \n";
  cout<<"*        (free format sphere descriptors) \n";
  cout<<"*        (free format cylinder descriptors) \n";
  
  ifstream fin(infile.c_str(), ios::in);
  svector V1, V2;
  int gar;
  while(!fin.eof()){
    fin>>gar>>V1>>V2;
    cout<<"5"<<endl;
    cout<<V1.x<<" "<<V1.y<<" "<<V1.z<<" "<<0.1<<" "<<V2.x<<" "<<V2.y<<" "<<V2.z<<" 0 0 1 0"<<endl;
  }
  
}
