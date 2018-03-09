//cfg2r3d

#include "r3d_util.h"
#include "../errors.h"

int main (int argc, char* argv[]){

  cout.precision(6);
  cout.setf(ios::fixed);

  string arg;
  string infile="temp_000000.cfg";
  int i, s;
  double Lscale=0;
  float theta=-90; float psi=0; float phi=0; //Euler angles (deg)
  string shadow="F";
  float eyepos=4.0; float radius=0.6;
  bool kecolor=0; bool nowalls=0; bool slice=0; bool bg=0;
  float colormax=-1000; float colormin=1000; float EkDiff=0;
  rgb RGB; float tsres=0.1; float tsrad=4; float tszres=0.1;
  svector Offset(0,0,0);
  bool drawtop=0; bool hideatom=0;
  float boxscale=4; int sizex=20; int sizey=20;
  svector light(1, 1, 1);
  double Lx, Ly, Lz;
  bool fcolor=0; bool ebcolor=0;

  if (argc==1){
    cerr<<"supply some more arguments.\n";
    exit(1);
  }
  
  for (i=1; i<argc; i++){
    arg=argv[i];
    if (arg.find(".cfg")!=string::npos){
      infile=arg;
    }
    else if (arg=="-ang"){
      theta=atof(argv[++i]);
      psi=atof(argv[++i]);
      phi=atof(argv[++i]);
    }
    else if (arg=="-nowalls") nowalls=1;
    else if (arg=="-slice") slice=1;
    else if (arg=="-bg") bg=1;
    else if (arg=="-size"){
      sizex=atoi(argv[++i]);
      sizey=atoi(argv[++i]);
    }
    else if (arg=="-tsres") tsres=atof(argv[++i]);
    else if (arg=="-tsrad") tsrad=atof(argv[++i]);
    else if (arg=="-rad") radius=atof(argv[++i]);
    else if (arg=="-topsurf") drawtop=1;
    else if (arg=="-hide") hideatom=1;
    else if (arg=="-if") infile=argv[++i];
    else if (arg=="-shadow") shadow="T";
    else if (arg=="-kecolor") kecolor=1;
    else if (arg=="-ebcolor") ebcolor=1;
    else if (arg=="-fcolor") fcolor=1;
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

  config cfg(infile);
  Lx=cfg.Lx; Ly=cfg.Ly; Lz=cfg.Lz;
  
  if (kecolor){
    cfg.Ek_sort_asc();
    colormin=cfg.begin->Ek();
    colormax=(cfg.end-1)->Ek();
  }

  if (fcolor){
    cfg.TimeStepInit();
    cfg.ReNeighbor();
    cfg.ForceEval(0);
    cfg.F_sort_asc();
    colormin=cfg.begin->F.mag();
    colormax=(cfg.end-1)->F.mag();
  }

  if (ebcolor){
    string ebfile=cfg.name+"_"+time2string(cfg.t)+".eb";
    if (FileExists(ebfile)){
      int a; float eb;
      fstream fin(ebfile.c_str(), ios::in);
      while (!fin.eof()){
	fin>>a>>eb;
	cfg.atomix(a)->P.x=eb;
	colormin=colormin <? eb;
	colormax=colormax >? eb;
      }
      colormin=-9;
      colormax=-1;
    }
    else{
      cerr<<"You need to create the eb profile for "<<infile<<endl;
      exit(1);
    }
  }

  //pick the dimension to scale in
  Lscale=fabs(Lx);
  if (fabs(Ly)>Lscale) Lscale=fabs(Ly);
  if (fabs(Lz)>Lscale) Lscale=fabs(Lz);
  boxscale*=Lscale;

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
  
  cout<<cfg.name<<"_"<<time2string(cfg.t)<<".r3d"<<endl;
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
  
  //draw spheres (type 2)
  if (!hideatom){
    vector<particle>::iterator I=cfg.begin;
    for (I=cfg.begin; I!=cfg.end; I++){
      if (!slice || (slice && I->R.y > 0)){
	cout<<"2 \n";
	Print(((I->R - Offset)/boxscale), cout);
	cout<<" "<<radius/boxscale<<" "; //radius
	if(I->id==14){
	  if(!I->is_fixed){
	    if(kecolor){
	      rgb kec;
	      kec.Set01((I->Ek()-colormin)/(colormax-colormin));
	      cout<<kec;
	    }
	    if (fcolor){
	      rgb fc;
	      fc.Set01((I->F.mag()-colormin)/(colormax-colormin));
	      cout<<fc;
	    }
	    if (ebcolor && (I->P.x != 0)){
	      rgb fc;
	      fc.Set01((I->P.x-colormin)/(colormax-colormin));
	      cout<<fc;
	    }
	    else cout<<"0 0 1"; //RGB color
	  }
	  else
	    cout<<"0 0 0.297"; //RGB color
	}
	if(I->id==18){
	  cout<<"1 1 0";
	}
	cout<<endl;
      }
    }
  }
  if(!nowalls){
    cout<<"8"<<endl;
    cout<<"20.0 -1.0   -1. -1. -1.   0.8  0 0 0 0"<<endl;
    
    svector Q1(-fabs(Lx)/2/boxscale,fabs(Ly)/2/boxscale,fabs(Lz)/2/boxscale);
    svector Q2(fabs(Lx)/2/boxscale,fabs(Ly)/2/boxscale,fabs(Lz)/2/boxscale);
    svector Q3(-fabs(Lx)/2/boxscale,-fabs(Ly)/2/boxscale,fabs(Lz)/2/boxscale);
    svector Q4(fabs(Lx)/2/boxscale,-fabs(Ly)/2/boxscale,fabs(Lz)/2/boxscale);
    svector Q5(-fabs(Lx)/2/boxscale,fabs(Ly)/2/boxscale,-fabs(Lz)/2/boxscale);
    svector Q6(fabs(Lx)/2/boxscale,fabs(Ly)/2/boxscale,-fabs(Lz)/2/boxscale);
    svector Q7(-fabs(Lx)/2/boxscale,-fabs(Ly)/2/boxscale,-fabs(Lz)/2/boxscale);
    svector Q8(fabs(Lx)/2/boxscale,-fabs(Ly)/2/boxscale,-fabs(Lz)/2/boxscale);
    
    rgb wc(0,1,0);
    //***** + yz wall***********
    Plane(Q2, Q4, Q6, Q8, wc);
    
    //***** - yz wall***********
    Plane(Q1, Q3, Q5, Q7, wc);
    
    //***** + xz wall***********
    Plane(Q1, Q2, Q5, Q6, wc);
    
    //***** - xz wall***********
    Plane(Q3, Q4, Q7, Q8, wc);
    
  }
  cout<<"9"<<endl; //ends transparent objects
  
  if (drawtop) TopSurf(cfg, boxscale, tsres, tsrad, tszres);
}
