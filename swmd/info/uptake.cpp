#include "../config.h"

//***************************************************************************
void InfoUptake(int argc, char * argv[]){
  string outfile;
  VPI i;
  int mono=0;
  float every=-1;
  vector<string>filelist;
  SPI s, s2;
  string arg, orig;
  bool stdout=0;
  for (int j=2; j<argc; j++){
    arg=argv[j];
    if (sfind(arg, ".cfg")) filelist.push_back(arg);
    else if (arg=="-mono") mono=atoi(argv[++j]);
    else if (arg=="-every") every=atof(argv[++j]);
    else if (arg=="-stdout") stdout=1;
    else if (arg=="-orig") orig=(argv[++j]);
  }
  int runs=filelist.size();
  int Nmax_orig;
  int N_orig;
  {
    config cfg(filelist.front());
    if (mono==0) mono=cfg.NFixed/2;
    outfile=cfg.name+".upt";
    Nmax_orig=cfg.Nmax;
    N_orig=cfg.N;
    for (i=cfg.begin; i<cfg.end; i++)
      if (i->id==9) N_orig--;
  }
  fstream fout;
  if (!stdout) fout.open(outfile.c_str(), ios::out);
  int inc=1;
  if (every != -1) inc=(int)(mono*every);
  for (int run=0; run<runs; run+=inc){
    int si_up=0; int hal_up=0; int si_sput=0;
    int sput=0;
    int sif=0; int sif2=0; int sif3=0;
    if (run<runs){
      config cfg(filelist[run]);
      for (i=cfg.begin; i<cfg.end; i++){
	if (i->id==14){
	  if (i->ix >= Nmax_orig) si_up++;
	  else sput++;
	  int nbs= 4 <? i->nlist.size();
	  int fnear=0;
	  for (int k=0; k < nbs; k++){
	    float nearest=1000;
	    for (s=i->nlist.begin(); s!=i->nlist.end(); s++)
	      if (((*s)->R - i->R).minsqmag(cfg.Lx, cfg.Ly) < nearest){
		nearest=((*s)->R - i->R).minsqmag(cfg.Lx, cfg.Ly);
		s2=s;
	      }
	    if ((*s2)->id==9) fnear++;
	    i->nlist.erase(s2);
	  }
	  if (fnear==1) sif++;
	  if (fnear==2) sif2++;
	  if (fnear==3) sif3++;
	}
	if (i->id==9) hal_up++;
      }
      si_sput=N_orig-sput;
      char outbuf[256];
      sprintf(outbuf,"%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n",(float)run/mono, 
	      (float)si_up/mono, (float)hal_up/mono,
	      (float)si_sput/mono, (float)sif/mono, (float)sif2/mono,
	      (float)sif3/mono);
      if (!stdout) fout<<outbuf;
      else cout<<outbuf;
    }
  }
}
