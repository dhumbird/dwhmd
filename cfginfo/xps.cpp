#include "i_config.h"

extern map<short,double> F_MIN, F_MAX;

//***************************************************************************
void InfoXPS(int argc, char * argv[]){
  string outfile;
  int mono=1;
  float every=-1;
  vector<string>filelist;
  string arg, orig, scfg;
  bool _stdout=0;
  bool parse=0; string logfile;
  bool abs=0; bool ignore=0;
  svector R; double r,f,f1;
  for (int q=2; q<argc; q++){
    arg=argv[q];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-mono") mono=atoi(argv[++q]);
    else if (arg=="-every") every=atof(argv[++q]);
    else if (arg=="-stdout") _stdout=1;
    else if (arg=="-orig") orig=(argv[++q]);
    else if (arg=="-abs") abs=1;
    else if (arg=="-ignore") ignore=1;
    else if (arg=="-parse"){
      parse=1;
      logfile=argv[++q];
    }
  }
  config cfg(scfg);
  //reduce neighbor list
  double cut=0;
  for (VAI i=cfg.begin; i<cfg.end; i++){
    for (SAI j=i->nlist.begin(); j!=i->nlist.end(); j++){
      switch (i->id + (*j)->id){
      case 12: 
	//cut=2.092; break;
	cut=3; break;
      case 7: 
	cut=2.25; break;
      case 2: 
	cut=1.21; break;
      case 20: 
	cut=4.141225; break;
      case 28: 
	cut=6.682; break;
      case 23: 
	cut=3.1006; break;
      case 15: 
	cut=1.941; break;
      }
      if ((i->R - (*j)->R).minsqmag(cfg.Lx, cfg.Ly) > cut){
	i->nlist.erase(j);
      }
    }
    if (i->id==6 || i->id ==14){
      while (i->nlist.size() > 4){
	float lmax=-10; SAI jmax;
	for (SAI j=i->nlist.begin(); j!=i->nlist.end(); j++){
	  if ((i->R - (*j)->R).minsqmag(cfg.Lx, cfg.Ly) > lmax){
	    lmax=(i->R - (*j)->R).minsqmag(cfg.Lx, cfg.Ly);
	    jmax=j;
	  }
	}
	i->nlist.erase(jmax);
      }
    }
    if (i->id==1){
      while (i->nlist.size() > 1){
	float lmax=-10; SAI jmax;
	for (SAI j=i->nlist.begin(); j!=i->nlist.end(); j++){
	  if ((i->R - (*j)->R).minsqmag(cfg.Lx, cfg.Ly) > lmax){
	    lmax=(i->R - (*j)->R).minsqmag(cfg.Lx, cfg.Ly);
	    jmax=j;
	  }
	}
	i->nlist.erase(jmax);
      }
    }
  }
  float res=2.5;
  float lambda=25;
  float ttop=-1000;
  float bottom=1000;
 
  for (VAI i=cfg.begin; i<cfg.end; i++){
    if (i->R.z > ttop) ttop=i->R.z;
    if (i->R.z < bottom) bottom=i->R.z;
  }
  int Nx=(int)(cfg.Lx/res)+1;
  int Ny=(int)(cfg.Ly/res)+1;
  res = cfg.Lx/(float)Nx;
  float top[Nx][Ny];
  for (int x=0; x<Nx; x++)
    for (int y=0; y<Ny; y++)
      top[x][y]=-1000;
  double Lx2=cfg.Lx/2;
  double Ly2=cfg.Ly/2;
  double nf=0; double nc=0;
  float CF, CF2, CF3, CC, CCF, SiC, SiF, SiSi;
  CF=CF2=CF3=CC=CCF=SiF=SiC=SiSi=0;
  
  for (VAI i=cfg.begin; i<cfg.end; i++){
    int x=(int)((Lx2-i->R.x)/res);
    int y=(int)((Ly2-i->R.y)/res);
    if (i->R.z > top[x][y]) top[x][y]=i->R.z;
  }
  map<short, int> Nmap;
  for (VAI i=cfg.begin; i<cfg.end; i++){
    Nmap[6]=Nmap[9]=Nmap[14]=0;
    int x=(int)((Lx2-i->R.x)/res);
    int y=(int)((Ly2-i->R.y)/res);
    //double depth=top[x][y] - i->R.z;
    double depth=ttop - i->R.z;
    double c=exp(-depth/lambda);
    for (SAI j=i->nlist.begin(); j!=i->nlist.end(); j++){
      int type=i->id+(*j)->id;
      Nmap[(*j)->id]++;
      if (i->ix < (*j)->ix){
	if (type==12){
	  bool ccf=0;
	  for (SAI k=(*j)->nlist.begin(); k!=(*j)->nlist.end(); k++){
	    if ((*k)->id==9) ccf=1;
	  }
	  for (SAI k=i->nlist.begin(); k!=i->nlist.end(); k++){
	    if ((*k)->id==9) ccf=1;
	  }
	  if (ccf) CCF+=c;
	  else CC+=c;
	}
	else if (type==20){
	  //  cerr<<i->ix<<" "<<(*j)->ix<<" "<<depth<<" "<<c<<endl;
	  SiC+=c;
	}
	else if (type==28){
	  SiSi+=c;
	}
      }
    }
    if (i->id==6){
      nc+=1;
      if (Nmap[9] == 3) CF3+=c;
      if (Nmap[9] == 2) CF2+=c;
      if (Nmap[9] == 1) CF+=c;
    }
    if (i->id==14){
      if (Nmap[9]>0){
	SiF+=c*Nmap[9];
      }
    }
    if (i->id==9) nf+=1;
  }
  //double f2c=(CF+2*CF2+3*CF3)/(CC+SiC+CCF+CF+CF2+CF3);
  double f2c=nf/nc;
  cout<<"C-Si       "<<SiC<<endl;
  cout<<"C-C        "<<CC<<endl;
  cout<<"C-CFx      "<<CCF<<endl;
  cout<<"CF         "<<CF<<endl;
  cout<<"CF2        "<<CF2<<endl;
  cout<<"CF3        "<<CF3<<endl;
  cout<<"Si-Si      "<<SiSi<<endl;
  cout<<"Si-F       "<<SiF<<endl;
  cout<<"C-C/Si     "<<SiC+CC<<endl;
  cout<<"F/C        "<<f2c<<endl;
}
