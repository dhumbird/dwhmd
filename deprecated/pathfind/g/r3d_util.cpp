#include "r3d_util.h"

void Plane(svector v1, svector v2, svector v3, svector v4, rgb c){
  //infallibly creates two triangles from four points

  map<float, short> rmap;
  rmap[(v1-v2).sqmag()]=12;
  rmap[(v1-v3).sqmag()]=13;
  rmap[(v1-v4).sqmag()]=14;
  rmap[(v2-v3).sqmag()]=23;
  rmap[(v2-v4).sqmag()]=24;
  rmap[(v3-v4).sqmag()]=34;
  map<float, short>::iterator m=rmap.end(); m--;
  switch (m->second){
  case 12:
    cout<<"1"<<endl; PrintV(v1); PrintV(v2); PrintV(v3); cout<<c<<endl;
    cout<<"1"<<endl; PrintV(v1); PrintV(v2); PrintV(v4); cout<<c<<endl;
    break;
  case 13:
    cout<<"1"<<endl; PrintV(v1); PrintV(v3); PrintV(v4); cout<<c<<endl;
    cout<<"1"<<endl; PrintV(v1); PrintV(v3); PrintV(v2); cout<<c<<endl;
    break;
  case 14:
    cout<<"1"<<endl; PrintV(v1); PrintV(v4); PrintV(v3); cout<<c<<endl;
    cout<<"1"<<endl; PrintV(v1); PrintV(v4); PrintV(v2); cout<<c<<endl;
    break;
  case 23:
    cout<<"1"<<endl; PrintV(v2); PrintV(v3); PrintV(v1); cout<<c<<endl;
    cout<<"1"<<endl; PrintV(v2); PrintV(v3); PrintV(v4); cout<<c<<endl;
    break;
  case 24:
    cout<<"1"<<endl; PrintV(v2); PrintV(v4); PrintV(v1); cout<<c<<endl;
    cout<<"1"<<endl; PrintV(v2); PrintV(v4); PrintV(v3); cout<<c<<endl;
    break;
  case 34:
    cout<<"1"<<endl; PrintV(v3); PrintV(v4); PrintV(v1); cout<<c<<endl;
    cout<<"1"<<endl; PrintV(v3); PrintV(v4); PrintV(v2); cout<<c<<endl;
    break;
  }
}
//*********************************************************************
void PrintV(svector& V){
  cout<<V.x<<" "<<V.y<<" "<<V.z<<" ";
}
//**********************************************************************
void TopSurf(config& cfg, vector<triangle>& topsurf, float xyres, 
	     float radscale, float zres){
  
  int Nx, Ny, i, j; 
  float radsq=pow(1.2*radscale, 2);
  //set up the input file
  string ts_file=cfg.name+"_"+time2string(cfg.t)+".ts";
  if (!FileExists(ts_file)){
    //set up the grid
    Nx=(int)(cfg.Lx/xyres)+1;
    Ny=(int)(cfg.Ly/xyres)+1;
    float box[Nx][Ny];
    float maxz=cfg.MaxZ()+5;
    for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) box[i][j]=maxz;
  
    //start stepping down. decrease the z-value in box[i][j] until it is within
    //sqrt(radsq) of any atom.
    svector R;
    bool flag;
    set<subcell*>::iterator c;
    set<particle*>::iterator j0;
    for (i=0; i<Nx; i++){
      cerr<<"\r topsurf: "<<(int)((float)i/Nx*100)<<" % done";
      for (j=0; j<Ny; j++){
	flag=1;
	R.x=i*xyres-cfg.Lx/2;
	R.y=cfg.Ly/2-j*xyres;
	while(flag){
	  box[i][j]-=zres;
	  R.z=box[i][j];
	  int x=(int)((R.x + cfg.Lx/2)/cfg.Lcx)%cfg.Nx;
	  int y=(int)(fabs((R.y - cfg.Ly/2)/cfg.Lcy))%cfg.Ny;
	  int z=(int)(fabs((R.z - cfg.the_top)/cfg.Lcz))%cfg.Nz;
	  int a=x+cfg.Nx*y+cfg.Nx*cfg.Ny*z;
	  for (c=cfg.cells[a].nbegin(); c!=cfg.cells[a].nend(); c++)
	    for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	      if(flag)
		flag=(R - (*j0)->R).minsqmag(cfg.Lx, cfg.Ly) > radsq;
	      else break;
	}
      }
    }
    cerr<<endl;
    //print out the output for next time
    ofstream fout(ts_file.c_str(), ios::out);
    fout<<Nx<<" "<<Ny<<" "<<xyres<<endl;
    for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) fout<<box[i][j]<<endl;
  }
  
  ifstream fin(ts_file.c_str(), ios::in);
  fin>>Nx>>Ny>>xyres;
  float box[Nx][Ny];

  for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) fin>>box[i][j];

  rgb pc(0,0,1);
  //output as raster3d triangles
  for (int i=0; i<Nx-1; i++)
    for (int j=0; j<Ny-1; j++){
      svector v1(i*xyres-cfg.Lx/2, cfg.Ly/2-j*xyres, box[i][j]);
      svector v2(i*xyres-cfg.Lx/2, cfg.Ly/2-(j+1)*xyres, box[i][j+1]);
      svector v3((i+1)*xyres-cfg.Lx/2, cfg.Ly/2-(j+1)*xyres, box[i+1][j+1]);
      svector v4((i+1)*xyres-cfg.Lx/2, cfg.Ly/2-j*xyres, box[i+1][j]);
      topsurf.push_back(triangle(v1, v2, v4, pc));
      topsurf.push_back(triangle(v3, v2, v4, pc));
    }
}
