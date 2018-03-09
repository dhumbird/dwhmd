#include "i_config.h"

//***************************************************************************
void InfoTraj(int argc, char * argv[]){
  string outfile;
  int mono=1;
  float every=-1;
  vector<string>filelist;
  string arg, orig;
  bool stdout=0;
  bool parse=0; string logfile;
  bool abs=0; bool ignore=0;
  float fby=1;
  bool all=0;
  svector R; double r,f,f1;
  int ix;
  for (int q=2; q<argc; q++){
    arg=argv[q];
    if (sfind(arg, ".cfg")) filelist.push_back(arg);
    else if (arg=="-ix") ix=atoi(argv[++q]);
    else if (arg=="-mono") mono=atoi(argv[++q]);
    else if (arg=="-every") every=atof(argv[++q]);
    else if (arg=="-stdout") stdout=1;
    else if (arg=="-parse"){
      parse=1;
      logfile=argv[++q];
    }
    else if (arg=="-all") all=1;
  }
  
  if (all && parse){
    string m, name;
    ifstream mdlog(logfile.c_str(), ios::in);
    while (!mdlog.eof()){
      getline(mdlog, m);
      if (sfind(m, "saved")){
	name=m.substr(m.find(":")+2, 1000);
	filelist.push_back(name);
      }
    }
  }
  int runs=filelist.size();
  
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
  outfile=itoa(ix,5)+".trj";
  if (!stdout) fout.open(outfile.c_str(), ios::out);
  int inc=1;
  if (every != -1) inc=(int)(mono*every);
 
  string header="";
  for (int q=0; q<(filelist.front()).length(); q++)
    header+=" ";
  header+="    Fl.";
  header+="   x    ";
  header+="  y    ";
  header+="  z";
  if (!stdout) fout<<header<<endl;
  else cout<<header<<endl;
  VAI s;
  char outbuf[256];
  for (int run=0; run<runs; run+=inc){
    if (run<runs){
      config cfg(filelist[run]);
      if ((s=cfg.atomix(ix)) != cfg.end){
	sprintf(outbuf,
	      "%6.2f %6.2f %6.2f %6.2f\n",
	      fluence[filelist[run]]*fby, s->R.x, s->R.y, s->R.z);
	if (!stdout) fout<<filelist[run]<<" "<<outbuf;
	else cout<<filelist[run]<<" "<<outbuf;
      }
    }
  }
}
