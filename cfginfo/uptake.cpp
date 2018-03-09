#include "i_config.h"

extern map<short,double> F_MIN, F_MAX;

struct saved_run{
  int run;
  float fluence;
  string filename;
  int pushes;
};

//***************************************************************************
void InfoUptake(int argc, char * argv[]){
  string outfile;
  int mono=1;
  float every=-1;
  vector<saved_run> saved_runs;
  string arg, orig;
  bool stdout=0;
  bool parse=0; string logfile="md.log";
  bool abs=0; bool ignore=0;
  int etchstart=0;
  float fby=1; bool dep=0;
  vector<string> filelist;
  svector R; double r,f,f1;
  for (int q=2; q<argc; q++){
    arg=argv[q];
    if (sfind(arg, ".cfg")) filelist.push_back(arg);
    else if (arg=="-mono") mono=atoi(argv[++q]);
    else if (arg=="-every") every=atof(argv[++q]);
    else if (arg=="-stdout") stdout=1;
    else if (arg=="-orig") orig=(argv[++q]);
    else if (arg=="-abs") abs=1;
    else if (arg=="-dep") dep=1;
    else if (arg=="-ignore") ignore=1;
    else if (arg=="-fby") fby=atof(argv[++q]);
    else if (arg=="-etch") etchstart=atoi(argv[++q]);
    else if (arg=="-parse"){
      parse=1;
      logfile=argv[++q];
    }
  }

  short sub_mat=14;
  int Nmax_orig; int Sub_max=0;
  int Nt_orig, Nsub_orig;
  map<short,int> Norig;
  double Lz_orig=0;

  if (!filelist.empty()){
    int count=1;
    for (vector<string>::iterator f=filelist.begin(); f!=filelist.end(); f++){
      saved_run s;
      s.filename=*f;
      s.run=count++;
      s.fluence=(float)s.run/mono*fby;
      s.pushes=0;
      saved_runs.push_back(s);
    }
  }
  else{
    string m, name; int arun;
    ifstream mdlog(logfile.c_str(), ios::in);
    int totalpush=0;
    while (!mdlog.eof()){
      getline(mdlog, m);
      if (sfind(m, "---run")){
	arun=atoi(((m.substr(m.find("run ")+4, 10))
		   .substr(0, (m.substr(m.find("run ")+4, 10))
			   .find("-"))).c_str());
      }
      if (sfind(m, "(14) pushed")){
	totalpush++;
      }
      if (sfind(m, "saved")){
	name=m.substr(m.find(":")+2, 1000);
	saved_run s;
	s.filename=name;
	s.run=arun;
	s.fluence=(float)s.run/mono*fby;
	s.pushes=totalpush;
	saved_runs.push_back(s);
      }
    }
  }
  
  string header="";
  Norig[1]=Norig[6]=Norig[9]=Norig[14]=Norig[17]=0;
  {
    config cfg((saved_runs.begin())->filename);
    if (mono==0) mono=cfg.NFixed/2;
    for (VAI i=cfg.begin; i<cfg.end; i++){
      if (i->is_fixed){ sub_mat=i->id; break;}
    }
    outfile=cfg.name+".upt";
    Nmax_orig=cfg.Nmax;
    for (VAI i=cfg.begin; i<cfg.end; i++){
      Norig[i->id]++;
      if (i->id==sub_mat && i->ix > Sub_max) Sub_max=i->ix;
    }
    Nt_orig=cfg.N;
    Nsub_orig=Norig[sub_mat];
    Lz_orig=cfg.Lz;
    for (int q=0; q<(saved_runs.begin())->filename.length(); q++)
      header+=" ";
  }

  fstream fout;
  if (!stdout) fout.open(outfile.c_str(), ios::out);
  int inc=1;
  if (every != -1) inc=(int)every;
  
  header+="    Fl.";
  header+="   C_up ";
  header+="  F_up ";
  header+="  Etch";
  header+="    SiF";
  header+="   SiF2";
  header+="   SiF3";
  header+="     CF";
  header+="    CF2";
  header+="    CF3";
  header+="     CC";
  header+="    CSi";
  if (!stdout) fout<<header<<endl;
  else cout<<header<<endl;
  vector<saved_run>::iterator s=saved_runs.begin();
  bool kill=0;
  while (s<saved_runs.end()){  
    int tpush=0;
    int c_up=0; int h_abs=0; int sub_sput=0;
    int sput=0; int h_tot=0; int h_new=0;
    int ch=0; int ch2=0; int ch3=0;
    int cf=0; int cf2=0; int cf3=0;
    int sif=0; int sif2=0; int sif3=0;
    int cc=0; int csi=0;
    int dangle=0; int f_tot=0; int sub_up=0;
    map<short, int> Nmap; int Nt;
    config cfg(s->filename);
    if (cfg.Lz != Lz_orig){
      for (VAI i=cfg.begin; i<cfg.end; i++){
	if (i->id==sub_mat && i->ix > Sub_max){
	  Norig[i->id]++;
	  Sub_max=i->ix;
	}
      }
      Nsub_orig=Norig[sub_mat];
      Lz_orig=cfg.Lz;
    }
    for (VAI i=cfg.begin; i<cfg.end; i++){
      Nmap[1]=Nmap[6]=Nmap[9]=Nmap[14]=Nmap[17]=0;
      Nt=0;
      for (SAI j=i->nlist.begin(); j!=i->nlist.end(); j++){
	int type=i->id+(*j)->id;
	Nmap[(*j)->id]++;
	Nt++;
	if ((*j)->ix > i->ix){
	  if (type==12) cc++;
	  if (type==20) csi++;
	}
      }
      if (i->id==sub_mat){
	if (i->ix > Sub_max) sub_up++;
	else sput++;
      }
      if (i->id==6){
	c_up++;
	if (!i->is_fixed && Nt < 4) dangle+=(4-Nt);
	if (Nmap[1] > 0){
	  if (!ignore || i->ix > Sub_max){
	    if (Nmap[1] == 1) ch++;
	    else if (Nmap[1] == 2) ch2++;
	    else if (Nmap[1] == 3) ch3++;
	  }
	}
	if (Nmap[9] > 0){
	  if (!ignore || i->ix > Sub_max){
	    if (Nmap[9] ==1) cf++;
	    else if (Nmap[9] == 2) cf2++;
	    else if (Nmap[9] == 3) cf3++;
	  }
	}
      }
      if (i->id == 14 && Nmap[9] > 0){
	if (!ignore || i->ix > Sub_max){
	  if (Nmap[9] == 1) sif++;
	  else if (Nmap[9] == 2) sif2++;
	  else if (Nmap[9] == 3) sif3++;
	}
	if (!i->is_fixed && Nt < 4) dangle+=(4-Nt);
      }
      if (i->id == 14 && Nmap[17] > 0){
	if (!ignore || i->ix > Sub_max){
	  if (Nmap[17] == 1) sif++;
	  else if (Nmap[17] == 2) sif2++;
	  else if (Nmap[17] == 3) sif3++;
	}
	if (!i->is_fixed && Nt < 4) dangle+=(4-Nt);
      }
      if (i->id==1){
	h_tot++;
	if (i->ix >= Nmax_orig) h_new++;
      }
      if (i->id==9||i->id==17) f_tot++;
    }
    sub_sput=Nsub_orig-sput-s->pushes;
    h_abs=(h_tot-h_new);
    float field2,field4;
    
    if (ignore) h_tot=ch+2*ch2+3*ch3;
    if (ignore) f_tot=cf+2*cf2+3*cf3+sif+2*sif2+3*sif3;

    if (dep) field2=(float)sub_up/mono;
    else field2=(float)c_up/mono;

    if (!abs) field4=(float)(sub_sput+etchstart)/mono;
    else field4=(float)h_abs/mono;
    char outbuf[256];
    sprintf(outbuf,
	    "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
	    s->fluence, 
	    (float)field2, (float)f_tot/mono, field4, 
	    (float)sif/mono, (float)sif2/mono, (float)sif3/mono, 
	    (float)cf/mono, (float)cf2/mono, (float)cf3/mono,
	    (float)cc/mono, (float)csi/mono);
    if (!stdout) fout<<s->filename<<" "<<outbuf;
    else cout<<s->filename<<" "<<outbuf;
    if (kill) break;
    if (s == saved_runs.end()-1) break;
    else if (s+inc >= saved_runs.end()){
      kill=1;
      s=saved_runs.end()-1;
    }
    else  s+=inc;
  }
}
