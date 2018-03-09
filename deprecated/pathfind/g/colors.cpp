
#include "../config.h"
#include "rgb.h"
#include <map>
#include <fstream>
#include "sphere.h"
#include "plane.h"
#include "cylinder.h"
#include "triangle.h"

void ReadScene(config& cfg, map<int, sphere>& spheres, bool colordef,
	       string colorfile, bool kecolor, bool keself, float colormin, 
	       float colormax, float radscale, bool zcolor,
	       vector<plane>& walls, vector<plane>& mask, plane& boxtop, 
	       float copyx,
	       vector<cylinder>& cylinders,
	       vector<triangle>& triangles, float atri){
  
  //***********************spheres****************************************
  for (VPI I=cfg.begin; I<cfg.end; I++){
    spheres[I->ix].R=I->R;
    spheres[I->ix].alpha=0;
    //default colors
    switch (I->id){
    case 14:
      if (!I->is_fixed)	spheres[I->ix].Color(0, 0, 1);
      else spheres[I->ix].Color(0, 0, 0.297);
      spheres[I->ix].radius=radscale*1.2;
      break;
    case 9:
      spheres[I->ix].Color(1, 1, 1);
      spheres[I->ix].radius=radscale*0.9;
      break;
    case 17:
      spheres[I->ix].Color(0, 1, 0);
      break;
    case 18:
      spheres[I->ix].radius=radscale*1.5;
      spheres[I->ix].Color(0, 1, 0);
      break;
    case 6:
      spheres[I->ix].radius=radscale;
      spheres[I->ix].Color(1, 0.5, 0);
      break;
    }
  }
  //color file overrides
  if (colordef){
    int ix; float r,g,b; float alpha;
    if (!FileExists(colorfile)){
      cerr<<"That color file not found! "<<colorfile<<endl;
      exit(1);
    }
    ifstream cdef(colorfile.c_str(), ios::in);
    while (!cdef.eof()){
      cdef>>ix>>r>>g>>b>>alpha;
      if (cfg.atomix(ix)!=cfg.end){
	spheres[ix].Color(r,g,b);
	spheres[ix].alpha=alpha;
      }
    }
  }
  //kinetic energy overrides (for silicon only)
  if (keself){
    kecolor=1;
    cfg.Ek_sort_asc();
    colormin=cfg.begin->Ek();
    colormax=(cfg.end-1)->Ek();
  }
  if (kecolor)
    for (VPI I=cfg.begin; I<cfg.end; I++)
      if (I->id==cfg.mat)
	spheres[I->ix].color.Set01((I->Ek()-colormin)/(colormax-colormin));
  if (zcolor){
    colormin=-cfg.Lz/2;
    colormax=cfg.Lz/2;
    for (VPI I=cfg.begin; I<cfg.end; I++)
      if (I->id==cfg.mat)
	spheres[I->ix].color.Set01((I->R.z-colormin)/(colormax-colormin));
  }
  //***********************planes***************************************
  vector<plane>::iterator W;
  svector Q1(-cfg.Lx/2, cfg.Ly/2, cfg.Lz/2);
  svector Q2(cfg.Lx/2, cfg.Ly/2, cfg.Lz/2);
  svector Q3(-cfg.Lx/2, -cfg.Ly/2, cfg.Lz/2);
  svector Q4(cfg.Lx/2, -cfg.Ly/2, cfg.Lz/2);
  svector Q5(-cfg.Lx/2, cfg.Ly/2, -cfg.Lz/2);
  svector Q6(cfg.Lx/2, cfg.Ly/2, -cfg.Lz/2);
  svector Q7(-cfg.Lx/2, -cfg.Ly/2, -cfg.Lz/2);
  svector Q8(cfg.Lx/2, -cfg.Ly/2, -cfg.Lz/2);
  rgb wc(0,1,0);
  walls.push_back(plane(Q2, Q4, Q6, Q8, wc));
  walls.push_back(plane(Q1, Q3, Q5, Q7, wc));
  walls.push_back(plane(Q1, Q2, Q5, Q6, wc));
  walls.push_back(plane(Q3, Q4, Q7, Q8, wc));
  for (W=walls.begin(); W!=walls.end(); W++)
    W->alpha=0.8;
  if (cfg.b_masked){
    mask.push_back(plane(svector(cfg.my_mask->finish.x,cfg.Ly/2,cfg.Lz),
			 svector(cfg.my_mask->finish.x,-cfg.Ly/2,cfg.Lz),
			 svector(cfg.my_mask->finish.x,-cfg.Ly/2,-cfg.Lz/2),
			 svector(cfg.my_mask->finish.x,cfg.Ly/2,-cfg.Lz/2),
			 wc));
    mask.push_back(plane(svector(cfg.my_mask->start.x,cfg.Ly/2,cfg.Lz),
			 svector(cfg.my_mask->start.x,-cfg.Ly/2,cfg.Lz),
			 svector(cfg.my_mask->start.x,-cfg.Ly/2,-cfg.Lz/2),
			 svector(cfg.my_mask->start.x,cfg.Ly/2,-cfg.Lz/2),
			 wc));
    if (copyx>0){
      mask.push_back(plane(svector(cfg.my_mask->start.x+cfg.Lx,cfg.Ly/2,cfg.Lz),
			   svector(cfg.my_mask->start.x+cfg.Lx,-cfg.Ly/2,cfg.Lz),
			   svector(cfg.my_mask->start.x+cfg.Lx,-cfg.Ly/2,-cfg.Lz/2),
			   svector(cfg.my_mask->start.x+cfg.Lx,cfg.Ly/2,-cfg.Lz/2),
			   wc));
    }
    for (W=mask.begin(); W!=mask.end(); W++)
      W->alpha=0.8;
  }
  boxtop=plane(svector(-(cfg.Lx+copyx)/2, cfg.Ly/2, cfg.Lz/2),
	       svector((cfg.Lx+copyx)/2, cfg.Ly/2, cfg.Lz/2),
	       svector(-(cfg.Lx+copyx)/2, -cfg.Ly/2, cfg.Lz/2),
	       svector((cfg.Lx+copyx)/2, -cfg.Ly/2, cfg.Lz/2),wc);
  boxtop.alpha=0.8;

  //************************cylinders***********************************

  vector<cylinder>::iterator C;
  cylinders.push_back(cylinder(Q1, Q2));
  cylinders.push_back(cylinder(Q1, Q3));
  cylinders.push_back(cylinder(Q1, Q5));
  cylinders.push_back(cylinder(Q2, Q4));
  cylinders.push_back(cylinder(Q2, Q6));
  cylinders.push_back(cylinder(Q3, Q4));
  cylinders.push_back(cylinder(Q3, Q7));
  cylinders.push_back(cylinder(Q4, Q8));
  cylinders.push_back(cylinder(Q5, Q6));
  cylinders.push_back(cylinder(Q5, Q7));
  cylinders.push_back(cylinder(Q6, Q8));
  cylinders.push_back(cylinder(Q7, Q8));
  for (C=cylinders.begin(); C!=cylinders.end(); C++){
    C->radius=0.1;
    C->Color(0,1,0);
    C->alpha=0;
  }
  
  //*******************triangles****************************************
  if (atri!=1000){
    svector T1(cfg.Lx/2+1, 0, atri);
    svector T2(cfg.Lx/2+6, 0 , atri-1);
    svector T3(cfg.Lx/2+6, 0, atri+1);
    triangles.push_back(triangle(T1, T2, T3, rgb(0,0,0)));
  }
}
