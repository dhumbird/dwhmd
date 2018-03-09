//cfg2r3d

#include <vector>
#include <map>
#include "rcfg.h"
#include "r3dfun.h"
#include "miscfun.h"
#include "d_sphere.h"
#include "d_plane.h"
#include "d_cylinder.h"
#include "d_triangle.h"

int main (int argc, char* argv[]){

  cout.precision(6);
  cout.setf(ios::fixed);

  string arg;
  string infile="temp_0000-000.cfg";
  double Lscale=0;
  float theta=90*PI/180; float psi=0; float phi=0; //Euler angles (deg)
  string shadow="T";
  float eyepos=4.0; float radius=1;
  bool kecolor=0; bool nowalls=1; float slice=-1000;
  float colormax=9; float colormin=0;
  svector Offset(0,0,0);
  float boxscale=4; int sizex=600; int sizey=600;
  svector light(0, 1, 1);
  bool showtop=0;
  bool colordef=0; string colorfile; 
  bool keself=0; short hide=-1;
  bool nolines=0;
  rgb bg(0, 0, 0);
  float atri=1000;
  bool drawbonds=0;
  bool zcolor=0; float b_rad=0.1;
  float bextend=0; bool asswipe=0;
  bool hideall=0;
  float beam_rad=0;

  if (argc==1){
    cerr<<"supply some more arguments.\n";
    exit(1);
  }
  
  for (int i=1; i<argc; i++){
    arg=argv[i];
    if (sfind(arg,"cfg")) infile=arg;
    else if (arg=="-ang"){
      theta=atof(argv[++i])*PI/180;;
      psi=atof(argv[++i])*PI/180;
      phi=atof(argv[++i])*PI/180;
    }
    else if (arg=="-hide") hide=atoi(argv[++i]);
    else if (arg=="-hideall") hideall=1;
    else if (arg=="-boxtop") showtop=1;
    else if (arg=="-zcolor") zcolor=1;
    else if (arg=="-walls") nowalls=0;
    else if (arg=="-beam") beam_rad=atof(argv[++i]);
    else if (arg=="-nolines") nolines=1;
    else if (arg=="-slice") slice=atof(argv[++i]);
    else if (arg=="-brad") b_rad=atof(argv[++i]);
    else if (arg=="-small"){
      sizex=300;
      sizey=300;
    }
    else if (arg=="-bg"){
      bg.r=atof(argv[++i]);
      bg.g=atof(argv[++i]);
      bg.b=atof(argv[++i]);
    }
    else if (arg=="-size"){
      sizex=atoi(argv[++i]);
      sizey=atoi(argv[++i]);
    }
    else if (arg=="-rad") radius=atof(argv[++i]);
    else if (arg=="-shadow") shadow="T";
    else if (arg=="+shadow") shadow="F";
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
    else if (arg=="-bonds") drawbonds=1;
    else if (arg=="-atri") atri=atof(argv[++i]);
    else if (arg=="-bextend") bextend=atof(argv[++i]);
    else if (arg=="-tpplot") asswipe=1;
    else CmdError(arg.c_str());
  }
  //*********************************
  string imagename;
  svector L;
  vector<triangle> zmarker;
  vector<triangle>::iterator T;
  vector<sphere> spheres;
  vector<sphere>::iterator I;
  vector<plane> walls, tpplot;
  vector<cylinder> boundary, bonds;
  vector<cylinder>::iterator C;
  cylinder beam;
  plane boxtop;
  vector<plane>::iterator W;
  {
    rcfg cfg(infile, radius);
    L=svector(cfg.Lx, cfg.Ly, cfg.Lz);
    if (!hideall){
      if (colordef) ReadColorFile(cfg, colorfile);
      if (kecolor) KEColor(cfg, colormin, colormax);
      if (zcolor) ZColor(cfg);
      FixSpheres(cfg, spheres, hide);
    }
    if (!nowalls) FixWalls(cfg, walls, bextend);
    if (!nolines) FixBoundary(cfg, boundary, bextend);
    if (atri!=1000) ZMarker(cfg, zmarker, atri);
    if (drawbonds) BondConnect(cfg, bonds, hide, b_rad);
    if (asswipe) AssWipe(cfg, tpplot, 0.1, 0.1, radius);
    if (beam_rad){
      beam.R1=svector(0,0,50);
      beam.R2=svector(0,0,cfg.the_bottom);
      beam.Color(1,.5,0);
      beam.alpha=0.9;
      beam.radius=beam_rad;
      
    }
    imagename=cfg.name+"_"+time2string(cfg.t)+".r3d";
  }
  //pick the dimension to scale in
  Lscale=L.x;
  //  if (L.y>Lscale) Lscale=L.y;
  //if (L.z>Lscale) Lscale=L.z;
  boxscale*=Lscale;
  L.EulerTrans(phi, theta, psi);
  L/=boxscale;
  float scale=2.5; //times normal size
  //start std output*******************************************
  //IMAGE HEADER
  cout<<imagename<<endl;  //file name
  cout<<sizex<<" "<<sizey<<endl;  //tiles in x,y
  cout<<"1 1"<<endl;  //pixels (x,y) per tile
  cout<<"4"<<endl;  //anti-aliasing level 4; 3x3->2x2
  cout<<bg<<endl;  //background color
  cout<<shadow<<endl;  //shadows cast(T/F);
  cout<<"25"<<endl;  //phong power
  cout<<"0.25"<<endl;  //secondary light contribution
  cout<<"0.05"<<endl;  //ambient light contribution
  cout<<"0.25"<<endl;  //specular reflection component
  cout<<eyepos<<endl;  //eye position
  //main light source position
  cout<<light.x<<" "<<light.y<<" "<<light.z<<endl;
  //view matrix
  cout<<"1 0 0 0"<<endl;
  cout<<"0 1 0 0"<<endl;
  cout<<"0 0 1 0"<<endl;
  cout<<"0 0 0 "<<1.0/scale<<endl;
  cout<<"3\n*\n*\n*\n";  //mixed object mode
  //cout<<"16\nFOG 1 0.7 0 0\n";
  //draw spheres (type 2)
  //self
  if (!spheres.empty()){
    for (I=spheres.begin(); I!=spheres.end(); I++){
      if (I->R.z > slice){
	I->Transform(phi, theta, psi, Offset, boxscale);
  	I->cout_r3d();
      }
    }
  }
  if (!walls.empty())
    for (W=walls.begin(); W!=walls.end(); W++){
      W->Transform(phi, theta, psi, Offset, boxscale);
      W->cout_r3d();
    }
  if (!tpplot.empty())
    for (W=tpplot.begin(); W!=tpplot.end(); W++){
      W->Transform(phi, theta, psi, Offset, boxscale);
      W->cout_r3d();
    }
  if (!boundary.empty())
    for (C=boundary.begin(); C!=boundary.end(); C++){
      C->Transform(phi, theta, psi, Offset, boxscale);
      C->cout_r3d();
    }
  if (!bonds.empty())
    for (C=bonds.begin(); C!=bonds.end(); C++){
      C->Transform(phi, theta, psi, Offset, boxscale);
      C->cout_r3d();
    }
  if (beam_rad){
    beam.Transform(phi, theta, psi, Offset, boxscale);
    beam.cout_r3d();
  }
  if (!zmarker.empty())
    for (T=zmarker.begin(); T!=zmarker.end(); T++){
      T->Transform(phi, theta, psi, Offset, boxscale);
      T->cout_r3d();
  }
}

