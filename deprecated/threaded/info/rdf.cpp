#include "rdf.h"

int RadDistFun(int argc, char * argv[]){

  string scfg="temp_000000.cfg";
  float runtime=1; float dr=0.1;
  double dt=0;
  int cfgevery=0;
  bool printe=0; bool tso=0;
  
  int i=0;
  string arg;

  for (i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-cfgevery") cfgevery=atoi(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-printe") printe=1;
    else if (arg=="-tso") tso=1;
    else if (arg=="-dr") dr=atof(argv[++i]);
    else CmdError(arg.c_str());
  }
  
  config cfg(scfg);
  VPI I, J;

  //set up inputs
  if (dt!=0) cfg.Setdt(dt);

  double start=cfg.t;
  double finish=start+runtime;
  double rij=0;
  int bin=0;
  int num_bins=(int)(fabs(cfg.Lx)/dr/2)+1;
  int hist[num_bins];
  double g[num_bins];
  for (i=0; i<num_bins; i++){
    hist[i]=0;
    g[i]=0.0;
  }
  
  cerr<<"Calculating RDF for "<<cfg.name<<"."<<endl;

  while (cfg.t<finish+cfg.dt){
    if (tso) cfg.dtOptimize();
    if (cfg.t!=start) cfg.FirstVV();
    if (cfgevery>0) if ((int)(cfg.t*1000)%cfgevery==0 || cfg.t>=finish)
      cfg.Dump(0);
    cfg.TimeStepInit();
    for (I=cfg.begin; I!=cfg.end-1; I++){
      for (J=I+1; J!=cfg.end; J++){
	rij=(I->R - J->R).minsqmag(cfg.Lx, cfg.Ly, cfg.Lz);
	bin=(int)(sqrt(rij)/dr);
	if (bin<num_bins)
	  if (!I->is_fixed && !J->is_fixed) hist[bin]+=2;
	if (rij < (I->rc >? J->rc)){
	  I->nlist.insert(J);
	  J->nlist.insert(I);
	}
      }
    }
    cfg.ForceEval(0);
    cfg.SecondVV();
    if(printe) cfg.PrintE();
    cfg.t+=cfg.dt;
  }
  cfg.t-=cfg.dt;

  ofstream fout;
  string rdf_file=cfg.name;
  rdf_file+=".rdf";
  fout.open(rdf_file.c_str(), ios::out);
  fout.fill(' ');
  fout.precision(12);
  fout.setf(ios::fixed);

  double c=4.0*PI*(double)cfg.N/fabs(cfg.Lx*cfg.Ly*cfg.Lz)/3.0;      
  for (i=0; i<num_bins; i++){
    double rlower=(double)i*dr;
    double rupper=rlower+dr;
    double nideal=c*(pow(rupper,3)-pow(rlower,3));
    g[i]=(double)hist[i]/((finish-start)/dt)/
      (double)(cfg.N-cfg.NFixed)/nideal;
    fout<<(double)i*dr<<"\t"<<g[i]<<endl;
  }
  cerr<<"RDF complete. Output file: "<<rdf_file<<"."<<endl;
  return(0);
}
