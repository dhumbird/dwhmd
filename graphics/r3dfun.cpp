#include "r3dfun.h"

void ReadColorFile(rcfg& cfg, string colorfile){
  string field1,r,g,b,alpha; VAI s;
  if (!FileExists(colorfile)){
    cerr<<"That color file not found! "<<colorfile<<endl;
    exit(1);
  }
  ifstream cdef(colorfile.c_str(), ios::in);
  while (!cdef.eof()){
    cdef>>field1>>r>>g>>b>>alpha;
    if (field1=="*"){
      for (VAI I=cfg.begin; I<cfg.end; I++){
	if (r!="*") I->color.r=atof(r.c_str());
	if (g!="*") I->color.g=atof(g.c_str());
	if (b!="*") I->color.b=atof(b.c_str());
	if (alpha!="*") I->alpha=atof(alpha.c_str());
      }
    }
    if (sfind(field1, "%")){
      if ((s=cfg.atomix(atoi((field1.replace(field1.find("%"),1,"")).c_str())))
	  !=cfg.end){
	for (SAI I=s->nlist.begin(); I!=s->nlist.end(); I++){
	  if (r!="*") (*I)->color.r=atof(r.c_str());
	  if (g!="*") (*I)->color.g=atof(g.c_str());
	  if (b!="*") (*I)->color.b=atof(b.c_str());
	  if (alpha!="*") (*I)->alpha=atof(alpha.c_str());
	}
      }
    }
    else{
      int match=0;
      if (field1=="Si" || field1=="si" || field1=="SI") match=14;
      if (field1=="C" || field1=="c") match=6;
      if (field1=="F" || field1=="f") match=9;
      if (field1=="Cl" || field1=="cl") match=17;
      if (match!=0){
	for (VAI I=cfg.begin; I<cfg.end; I++){
	  if (I->id==match){
	    if (r!="*") I->color.r=atof(r.c_str());
	    if (g!="*") I->color.g=atof(g.c_str());
	    if (b!="*") I->color.b=atof(b.c_str());
	    if (alpha!="*") I->alpha=atof(alpha.c_str());
	  }
	}
      }
      else{
	if ((s=cfg.atomix(atoi(field1.c_str())))!=cfg.end){
	  if (r!="*") s->color.r=atof(r.c_str());
	  if (g!="*") s->color.g=atof(g.c_str());
	  if (b!="*") s->color.b=atof(b.c_str());
	  if (alpha!="*") s->alpha=atof(alpha.c_str());
	}
      }
    }
  }
}

void KEColor(rcfg& cfg, float colormin, float colormax){
  for (VAI I=cfg.begin; I<cfg.end; I++)
    I->color.Set01((I->Ek()-colormin)/(colormax-colormin));
}

void ZColor(rcfg& cfg){
  float colormin=-cfg.Lz/2;
  float colormax=cfg.Lz/2;
  for (VAI I=cfg.begin; I<cfg.end; I++){
    I->color.Set01((I->R.z-colormin)/(colormax-colormin));
    //cerr<<I->ix<<" "<<I->color<<" 0"<<endl;
  }
}

void FixSpheres(rcfg& cfg, vector<sphere>& spheres, short hide){
  for (VAI s=cfg.begin; s<cfg.end; s++){
    if (s->id!=hide){
      sphere sph(s->R);
      sph.color=s->color;
      sph.radius=s->radius;
      sph.alpha=s->alpha;
      spheres.push_back(sph);
    }
  }
}

void FixWalls(rcfg& cfg, vector<plane>& walls, float bextend){
  vector<plane>::iterator W;
  double z_high=cfg.the_bottom+cfg.Lz+bextend;
  double z_low=cfg.the_bottom;
  svector Q1(-cfg.Lx/2, cfg.Ly/2, z_high);
  svector Q2(cfg.Lx/2, cfg.Ly/2, z_high);
  svector Q3(-cfg.Lx/2, -cfg.Ly/2, z_high);
  svector Q4(cfg.Lx/2, -cfg.Ly/2, z_high);
  svector Q5(-cfg.Lx/2, cfg.Ly/2, z_low);
  svector Q6(cfg.Lx/2, cfg.Ly/2, z_low);
  svector Q7(-cfg.Lx/2, -cfg.Ly/2, z_low);
  svector Q8(cfg.Lx/2, -cfg.Ly/2, z_low);
  rgb wc(0,1,0);
  walls.push_back(plane(Q2, Q4, Q6, Q8, wc));
  walls.push_back(plane(Q1, Q3, Q5, Q7, wc));
  walls.push_back(plane(Q1, Q2, Q5, Q6, wc));
  walls.push_back(plane(Q3, Q4, Q7, Q8, wc));
  for (W=walls.begin(); W!=walls.end(); W++) W->alpha=0.8;
}

void FixBoundary(rcfg& cfg, vector<cylinder>& cylinders, float bextend){
  vector<cylinder>::iterator C;
  double z_high=cfg.the_bottom+cfg.Lz+bextend;
  double z_low=cfg.the_bottom;
  svector Q1(-cfg.Lx/2, cfg.Ly/2, z_high);
  svector Q2(cfg.Lx/2, cfg.Ly/2, z_high);
  svector Q3(-cfg.Lx/2, -cfg.Ly/2, z_high);
  svector Q4(cfg.Lx/2, -cfg.Ly/2, z_high);
  svector Q5(-cfg.Lx/2, cfg.Ly/2, z_low);
  svector Q6(cfg.Lx/2, cfg.Ly/2, z_low);
  svector Q7(-cfg.Lx/2, -cfg.Ly/2, z_low);
  svector Q8(cfg.Lx/2, -cfg.Ly/2, z_low);
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
}

void ZMarker(rcfg& cfg, vector<triangle>& triangles, float zpos){
  svector T1(cfg.Lx/2+1, 0, zpos);
  svector T2(cfg.Lx/2+6, 0, zpos-1);
  svector T3(cfg.Lx/2+6, 0, zpos+1);
  triangles.push_back(triangle(T1, T2, T3, rgb(0,0,0)));
}

void BondConnect(rcfg& cfg, vector<cylinder>& cylinders, short hide, 
		 float b_rad){
  cylinder cc;
  rgb bcolor;
  for (VAI I=cfg.begin; I<cfg.end; I++){
    for (SAI j0=I->nlist.begin(); j0!=I->nlist.end(); j0++){
      if (I->id!=hide && (*j0)->id!=hide){
	if (*j0 > &(*I)){
	  //I side
	  svector Ij= I->R-(*j0)->R;
	  svector Rhat=Ij/Ij.mag();
	  if (Ij.minsqmag(cfg.Lx,cfg.Ly) == Ij.sqmag()){
	    cc=cylinder(I->R-Rhat*I->radius, I->R-Ij/2, b_rad);
	    cc.color=I->color;
	    cc.alpha=I->alpha;
	    cylinders.push_back(cc);
	    cc=cylinder((*j0)->R+Rhat*(*j0)->radius, (*j0)->R+Ij/2, b_rad);
	    cc.color=(*j0)->color;
	    cc.alpha=(*j0)->alpha;
	    cylinders.push_back(cc);
	  }
	  else{ 
	    svector RR=(I->R - (*j0)->R);
	    cylinder C1, C2;
	    C1.R1=I->R; C1.R2=(*j0)->R; C1.color=bcolor; C1.radius=b_rad;
	    C1.alpha=(*j0)->alpha;
	    C2.R1=(*j0)->R; C2.R2=I->R; C2.color=bcolor; C2.radius=b_rad;
	    C2.alpha=(*j0)->alpha;
	    if (fabs(RR.x) > cfg.Lx/2){
	      if (I->R.x < (*j0)->R.x){
		C1.R2.x-=cfg.Lx;
		C2.R2.x+=cfg.Lx;
	      }
	      else{
		C1.R2.x+=cfg.Lx;
		C2.R2.x-=cfg.Lx;
	      }
	    }
	    if (fabs(RR.y) > cfg.Ly/2){
	      if (I->R.y < (*j0)->R.y){
		C1.R2.y-=cfg.Ly;
		C2.R2.y+=cfg.Ly;
	      }
	      else{
		C1.R2.y+=cfg.Ly;
		C2.R2.y-=cfg.Ly;
	      }
	    }
	    if (C1.R2.x > cfg.Lx/2) C1.R2.x=cfg.Lx/2;
	    if (C1.R2.x < -cfg.Lx/2) C1.R2.x=-cfg.Lx/2;
	    if (C1.R2.y > cfg.Ly/2) C1.R2.y=cfg.Ly/2;
	    if (C1.R2.y < -cfg.Ly/2) C1.R2.y=-cfg.Ly/2;
	    if (C2.R2.x > cfg.Lx/2) C2.R2.x=cfg.Lx/2;
	    if (C2.R2.x < -cfg.Lx/2) C2.R2.x=-cfg.Lx/2;
	    if (C2.R2.y > cfg.Ly/2) C2.R2.y=cfg.Ly/2;
	    if (C2.R2.y < -cfg.Ly/2) C2.R2.y=-cfg.Ly/2;
	    cylinders.push_back(C1);
	    cylinders.push_back(C2);
	  }
 	}
      }
    }
  }
}

void AssWipe(rcfg& cfg, vector<plane>& tp, float res, float zres, float radscale){
  float xyres=res;
  int Nx, Ny, i, j; 
  
  string ts_file=cfg.name+"_"+time2string(cfg.t)+".ts";
  if (!FileExists(ts_file)){
    Nx=(int)(cfg.Lx/xyres)+1;
    Ny=(int)(cfg.Ly/xyres)+1;
    float box[Nx][Ny];
    float maxz=cfg.the_top+5;
    for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) box[i][j]=maxz;
    
    //start stepping down. decrease the z-value in box[i][j] until it is within
    //sqrt(radsq) of any atom.
    svector R;
    bool flag;
    VAI j0;
    for (i=0; i<Nx; i++){
      cerr<<"\r topsurf: "<<(int)((float)i/Nx*100)<<" % done";
      for (j=0; j<Ny; j++){
	flag=1;
	R.x=i*xyres-cfg.Lx/2;
	R.y=cfg.Ly/2-j*xyres;
	while(flag){
	  box[i][j]-=zres;
	  R.z=box[i][j];
	  for (j0=cfg.begin; j0!=cfg.end; j0++)
	    if(flag)
	      flag=(R - j0->R).minsqmag(cfg.Lx, cfg.Ly) > 
		pow(j0->radius*radscale,2);
	    else break;
	}
      }
    }
    cerr<<endl;
    ofstream fout(ts_file.c_str(), ios::out);
    fout<<Nx<<" "<<Ny<<" "<<xyres<<endl;
    for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) fout<<box[i][j]<<endl;
  }
  
  ifstream fin(ts_file.c_str(), ios::in);
  fin>>Nx>>Ny>>xyres;
  float box[Nx][Ny];
  for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) fin>>box[i][j];

  rgb pc(1, 1, 1);
  //output as raster3d triangles
  for (int i=0; i<Nx-1; i++)
    for (int j=0; j<Ny-1; j++){
      svector v1((i*xyres-cfg.Lx/2), (cfg.Ly/2-j*xyres), box[i][j]);
      svector v2((i*xyres-cfg.Lx/2), (cfg.Ly/2-(j+1)*xyres), box[i][j+1]);
      svector v3(((i+1)*xyres-cfg.Lx/2),(cfg.Ly/2-(j+1)*xyres), box[i+1][j+1]);
      svector v4(((i+1)*xyres-cfg.Lx/2),(cfg.Ly/2-j*xyres), box[i+1][j]);
      tp.push_back(plane(v1, v2, v3, v4, pc));
    }
  for (vector<plane>::iterator W=tp.begin(); W!=tp.end(); W++)
    W->alpha=0;
}
