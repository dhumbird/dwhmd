#include "../config.h"

void FindAmorph(int argc, char * argv[]){
  string outfile,arg;
  int mono=0;
  int Ix=8; int Iy=8; int Iz=6;
  bool profile=0;
  float cutoff=0.9;
  float every=-1;
  bool trans=1; bool print=0;
  vector<string> filelist;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if(sfind(arg,"cfg")) filelist.push_back(arg);
    else if(arg=="-I"){ Ix=atoi(argv[++i]); Iy=Ix; Iz=Ix;}
    else if(arg=="-Ix") Ix=atoi(argv[++i]);
    else if(arg=="-Iy") Iy=atoi(argv[++i]);
    else if(arg=="-Iz") Iz=atoi(argv[++i]);
    else if (arg=="-every") every=atof(argv[++i]);
    else if (arg=="-mono") mono=atoi(argv[++i]);
    else if(arg=="-profile") profile=1;
    else if(arg=="-cutoff") cutoff=atof(argv[++i]);
    else if(arg=="-notrans") trans=0;
    else if (arg=="-print") print=1;
  }
  
  vector<svector> cryst;
  vector<svector>::iterator vsi;
  VPI i;
  float rho=0.049852723;
  double a=pow(8/rho, ONE_THIRD);
  //float norm=0.25*a;
  float norm=sqrt(3.0)/4.0*a;
  {
    config cfg(filelist.front());
    Iz=int(cfg.Lz_full/a)+1;
    if (mono==0) mono=cfg.NFixed/2;
    outfile=cfg.name+".amo";
  
    //make the ghost lattice******************************************
    cryst.clear();
    for (short x=0; x<Ix; x++)
      for (short y=0; y<Iy; y++)
	for (short z=0; z<Iz; z++){
	  cryst.push_back(svector((x+0.25)*a, (y+0.75)*a, (z+0.75)*a));
	  cryst.push_back(svector((x+0.75)*a, (y+0.25)*a, (z+0.75)*a));
	  cryst.push_back(svector(x*a, (y+0.5)*a, (z+0.5)*a));
	  cryst.push_back(svector((x+0.5)*a, (y+1)*a, (z+0.5)*a));
	  cryst.push_back(svector((x+0.25)*a, (y+0.25)*a, (z+0.25)*a));
	  cryst.push_back(svector((x+0.75)*a, (y+0.75)*a, (z+0.25)*a));
	  cryst.push_back(svector(x*a, y*a, z*a));     
	  cryst.push_back(svector((x+0.5)*a, (y+0.5)*a, z*a));
	}
    //reset ghost lattice origin**************************************
    for (vsi=cryst.begin(); vsi!=cryst.end(); vsi++){
      vsi->x-=cfg.Lx/2;
      vsi->y-=cfg.Ly/2;
      vsi->z-=cfg.Lz/2;
    }
  }
  
  //finally, produce the order profile************************************
  if (profile && filelist.size()==1){
    config cfg(filelist.front());
    float my_top=cfg.AvgTop();
    float res=0.25; //a factor times a
    for (float z=0; z<Iz; z+=res){
      float mean=0;
      int count=0;
      for (i=cfg.begin; i<cfg.end; i++){
	if (i->R.z >= z*a-cfg.Lz/2-0.0001 && i->R.z < (z+res)*a-cfg.Lz/2-0.0001){
	  float d=10000;
	  for (vsi=cryst.begin(); vsi!=cryst.end(); vsi++){
	    d = d <? (i->R - *vsi).minsqmag(cfg.Lx,cfg.Ly);
	  }
	  mean+=d;
	  count++;
	}
      }
      mean = (count>0)?sqrt(mean/count)/norm:0;
      if (mean!=0 && my_top-((z)*a-cfg.Lz/2)>0)
	printf("%f\t%f\n", my_top-((z)*a-cfg.Lz/2), mean);
    }
  }
  else{
    int runs=filelist.size();
    fstream fout(outfile.c_str(), ios::out);
    char outbuf[256];
    int inc=1; float res=0.25; //a length in ang.
    if (every != -1) inc=(int)(mono*every);
    for (int run=0; run<runs; run+=inc){
      if (run<runs){
	float xi=10;
	config cfg(filelist[run]);
	float my_top=cfg.AvgTop();
	float z=my_top;
	float minz=cfg.MinZ();
	while (z > minz && xi > cutoff*1.001){
	  float mean=0;
	  int count=0;
	  for (i=cfg.begin; i<cfg.end; i++)
	    if (i->R.z >= z-a/4 && i->R.z < z+a/4){
	      float d=10000;
	      for (vsi=cryst.begin(); vsi!=cryst.end(); vsi++)
		d = d <? (i->R - *vsi).minsqmag(cfg.Lx,cfg.Ly);
	      mean+=d;
	      count++;
	    }
	  xi = (count>0)?sqrt(mean/count)/norm:0;
	  z-=res;
	}
	z+=res;
	if (z > minz)
	  sprintf(outbuf,"%f\t%f\n", (float)run/mono, my_top-z);
	else sprintf(outbuf,"%f\t%f\n", (float)run/mono, 0);
	if (!print) fout<<outbuf;
	else cout<<outbuf;
      }
    }
  }
}


