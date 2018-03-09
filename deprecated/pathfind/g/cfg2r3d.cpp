//cfg2r3d

#include "r3d_util.h"
#include "../errors.h"
#include "sphere.h"
#include "plane.h"
#include "cylinder.h"
#include "triangle.h"

void ReadScene(config&, map<int, sphere>&, bool, string, bool, bool, float, 
	       float, float, bool, vector<plane>&, vector<plane>&, plane&, 
	       float, vector<cylinder>&, vector<triangle>&, float);

int main (int argc, char* argv[]){

  cout.precision(6);
  cout.setf(ios::fixed);

  string arg;
  string infile="temp_00000000.cfg";
  double Lscale=0;
  float theta=90*PI/180; float psi=0; float phi=0; //Euler angles (deg)
  string shadow="F";
  float eyepos=4.0; float radius=1;
  bool kecolor=0; bool nowalls=0; float slice=1000;
  float colormax=9; float colormin=0;
  svector Offset(0,0,0);
  bool showmask=1;
  float boxscale=4; int sizex=600; int sizey=600;
  svector light(1, 1, 1);
  double Lx, Ly, Lz;
  bool showtop=0;
  float copyx=0; bool colordef=0;
  bool keself=0;
  bool nolines=0;
  map<int,sphere> spheres;
  map<int, sphere>::iterator I;
  vector<plane> walls, mask;
  vector<cylinder> cyls;
  vector<cylinder>::iterator C;
  plane boxtop;
  vector<plane>::iterator W;
  rgb bg(0, 0, 0);
  bool fill=0;
  vector<triangle> triangles, topsurf;
  vector<triangle>::iterator T;
  float atri=1000;
  string colorfile;
  bool zcolor=0;
  float tsres=0.1; float tszres=0.1; bool drawtop=0; bool hideatom=0;
    
  if (argc==1){
    cerr<<"supply some more arguments.\n";
    exit(1);
  }
  
  for (int i=1; i<argc; i++){
    arg=argv[i];
    if (sfind(arg,"cfg")){
      infile=arg;
    }
    else if (arg=="-ang"){
      theta=atof(argv[++i])*PI/180;;
      psi=atof(argv[++i])*PI/180;
      phi=atof(argv[++i])*PI/180;
    }
    else if (arg=="-boxtop") showtop=1;
    else if (arg=="-nomask") showmask=0;
    else if (arg=="-zcolor") zcolor=1;
    else if (arg=="-nowalls") nowalls=1;
    else if (arg=="-nolines") nolines=1;
    else if (arg=="-slice") slice=atof(argv[++i]);
    else if (arg=="-bg"){
      bg.r=atof(argv[++i]);
      bg.g=atof(argv[++i]);
      bg.b=atof(argv[++i]);
    }
    else if (arg=="-size"){
      sizex=atoi(argv[++i]);
      sizey=atoi(argv[++i]);
    }
    else if (arg=="-fill") fill=1;
    else if (arg=="-copyx") copyx=atof(argv[++i]);
    else if (arg=="-tsres") tsres=atof(argv[++i]);
    else if (arg=="-rad") radius=atof(argv[++i]);
    else if (arg=="-shadow") shadow="T";
    else if (arg=="-topsurf") drawtop=1;
    else if (arg=="-hide") hideatom=1;
    else if (arg=="-kecolorself") keself=1;
    else if (arg=="-kecolor"){
      kecolor=1;
      colormin=atof(argv[++i]);
      colormax=atof(argv[++i]);
    }
    else if (arg=="-eyepos") eyepos=atof(argv[++i]);
    else if (arg=="-boxscale") boxscale=atof(argv[++i]);
    else if (arg=="-offset")
      Offset=svector(atof(argv[++i]), atof(argv[++i]), atof(argv[++i]));
    else if (arg=="-light")
      light=svector(atof(argv[++i]), atof(argv[++i]), atof(argv[++i]));
    else if (arg=="-colordef"){
      colordef=1;
      colorfile=argv[++i];
    }
    else if (arg=="-atri")
      atri=atof(argv[++i]);
    else CmdError(arg.c_str());
  }
  //*********************************
  string imagename;
  svector L;
  {
    config cfg(infile);
    L.x=cfg.Lx; L.y=cfg.Ly; L.z=cfg.Lz;
    if (!cfg.b_masked) showmask=0;
    ReadScene(cfg, spheres, colordef, colorfile, kecolor, keself, 
	      colormin, colormax, radius, zcolor, walls, mask, boxtop, copyx, 
	      cyls, triangles, atri);
    if (drawtop) TopSurf(cfg, topsurf, tsres, radius, tszres);
    imagename=cfg.name+"_"+time2string(cfg.t)+".r3d";
  }
  //pick the dimension to scale in
  Lscale=L.x;
  if (L.y>Lscale) Lscale=L.y;
  if (L.z>Lscale) Lscale=L.z;
  boxscale*=Lscale;
  L.EulerTrans(phi, theta, psi);
  L/=boxscale;
  float scale=2.5; //times normal size
  for (I=spheres.begin(); I!=spheres.end(); I++)
    I->second.Transform(phi, theta, psi, Offset, boxscale);
  for (W=walls.begin(); W!=walls.end(); W++)
    W->Transform(phi, theta, psi, Offset, boxscale);
  for (W=mask.begin(); W!=mask.end(); W++)
    W->Transform(phi, theta, psi, Offset, boxscale);
  boxtop.Transform(phi, theta, psi, Offset, boxscale);
  for (C=cyls.begin(); C!=cyls.end(); C++)
    C->Transform(phi, theta, psi, Offset, boxscale);
  for (T=triangles.begin(); T!=triangles.end(); T++)
    T->Transform(phi, theta, psi, Offset, boxscale);
  for (T=topsurf.begin(); T!=topsurf.end(); T++)
    T->Transform(phi, theta, psi, Offset, boxscale);
  //start std output*******************************************
  //IMAGE HEADER
  //file name
  cout<<imagename<<endl;
  //tiles in x,y
  cout<<sizex<<" "<<sizey<<endl;
  //pixels (x,y) per tile
  cout<<"1 1"<<endl;
  //anti-aliasing level 4; 3x3->2x2
  cout<<"4"<<endl;
  //background color
  cout<<bg<<endl;
  //shadows cast(T/F);
  cout<<shadow<<endl;
  //phong power
  cout<<"25"<<endl;
  //secondary light contribution
  cout<<"0.25"<<endl;
  //ambient light contribution
  cout<<"0.05"<<endl;
  //specular reflection component
  cout<<"0.25"<<endl;
  //eye position
  cout<<eyepos<<endl;
  //main light source position
  cout<<light.x<<" "<<light.y<<" "<<light.z<<endl;
  //view matrix
  cout<<"1 0 0 0"<<endl;
  cout<<"0 1 0 0"<<endl;
  cout<<"0 0 1 0"<<endl;
  cout<<"0 0 0 "<<1.0/scale<<endl;
  //mixed object mode
  cout<<"3\n*\n*\n*\n";

  //draw spheres (type 2)
  if (!fill){
    if (!hideatom)
      for (I=spheres.begin(); I!=spheres.end(); I++){
	if (I->second.R.z < slice){
	  I->second.cout_r3d();
	  if (copyx!=0)
	    if (I->second.R.x < (-Lx/2+copyx))
	      I->second.cout_r3d();
	}
      }
  }
    
  if (!nowalls)
    for (W=walls.begin(); W!=walls.end(); W++)
      W->cout_r3d();
  if (showmask)
    for (W=mask.begin(); W!=mask.end(); W++)
      W->cout_r3d();
  if (showtop)
    boxtop.cout_r3d();
  if (!nolines)
    for (C=cyls.begin(); C!=cyls.end(); C++)
      C->cout_r3d();
  for (T=triangles.begin(); T!=triangles.end(); T++)
    T->cout_r3d();
  if (drawtop)
    for (T=topsurf.begin(); T!=topsurf.end(); T++)
    T->cout_r3d();
}

