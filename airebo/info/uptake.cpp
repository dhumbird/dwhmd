#include "../config.h"

//***************************************************************************
void InfoUptake(int argc, char * argv[]){
  string outfile;
  VAI i,j,k;
  int mono=1;
  float every=-1;
  vector<string>filelist;
  VNI s, s2;
  string arg, orig;
  bool stdout=0;
  bool parse=0; string logfile;
  bool abs=0;
  for (int q=2; q<argc; q++){
    arg=argv[q];
    if (sfind(arg, ".cfg")) filelist.push_back(arg);
    else if (arg=="-mono") mono=atoi(argv[++q]);
    else if (arg=="-every") every=atof(argv[++q]);
    else if (arg=="-stdout") stdout=1;
    else if (arg=="-orig") orig=(argv[++q]);
    else if (arg=="-abs") abs=1;
    else if (arg=="-parse"){
      parse=1;
      logfile=argv[++q];
    }
  }
  int runs=filelist.size();
  int Nmax_orig;
  int C_orig=0; int N_orig;
  int H_orig=0;
  {
    config cfg(filelist.front());
    if (mono==0) mono=cfg.NFixed/2;
    outfile=cfg.name+".upt";
    Nmax_orig=cfg.Nmax;
    for (i=cfg.begin; i<cfg.end; i++){
      if (i->id==1) H_orig++;
      else C_orig++;
    }
    N_orig=C_orig;
  }
  map<string,float> fluence;
  if (parse){
    string m, name; int run;
    ifstream mdlog(logfile.c_str(), ios::in);
    while (!mdlog.eof()){
      getline(mdlog, m);
      if (sfind(m, "---run"))
	run=atoi(((m.substr(m.find("run ")+4, 10))
		  .substr(0, (m.substr(m.find("run ")+4, 10))
			  .find("-"))).c_str());
      if (sfind(m, "saved")){
	name=m.substr(m.find(":")+2, 1000);
	fluence[name]=(float)run/mono;
      }
    }
  }
  else{
    int run=0;
    for (vector<string>::iterator vsi=filelist.begin(); vsi!=filelist.end();
	 vsi++){
      fluence[*vsi]=(float)run/mono;
      run++;
    }
  }
  fstream fout;
  if (!stdout) fout.open(outfile.c_str(), ios::out);
  int inc=1;
  if (every != -1) inc=(int)(mono*every);
  for (int run=0; run<runs; run+=inc){
    int c_up=0; int h_abs=0; int c_sput=0;
    int sput=0; int h_tot=0; int h_new=0;
    int ch=0; int ch2=0; int ch3=0;
    int dangle=0;
    if (run<runs){
      config cfg(filelist[run]);
      for (i=cfg.begin; i<cfg.end; i++){
	if (i->id==6){
	  if (i->ix >= Nmax_orig) c_up++;
	  else sput++;
	  if (i->NH > 0){
	    if (i->NH > 0.3 && i->NH < 1.3) ch++;
	    else if (i->NH > 1.3 && i->NH < 2.3) ch2++;
	    else if (i->NH > 2.3) ch3++;
	  }
	  if (!i->is_fixed && i->Nt < 3.3) dangle++;
	}
	if (i->id==1){
	  h_tot++;
	  if (i->ix >= Nmax_orig) h_new++;
	}
      }
      c_sput=N_orig-sput;
      h_abs=(h_tot-h_new);
      float field4;
      if (!abs) field4=(float)c_sput/mono;
      else field4=(float)h_abs/mono;
      char outbuf[256];
      sprintf(outbuf,
	      "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
	      fluence[filelist[run]], 
	      (float)c_up/mono, (float)h_tot/mono, field4, 
	      (float)ch/mono, (float)ch2/mono, (float)ch3/mono, 
	      (float)dangle/mono);
      if (!stdout) fout<<filelist[run]<<" "<<outbuf;
      else cout<<filelist[run]<<" "<<outbuf;
    }
  }
}
