#include "config.h"
#include "errors.h"

int main(int argc, char * argv[]){

  string scfg="temp_0000-000.cfg";
  float runtime=1; float dr=0.1;
  double dt=0;
  bool printe=0; bool tso=0;
  
  int i=0;
  string arg;

  for (i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-printe") printe=1;
    else if (arg=="-tso") tso=1;
    else if (arg=="-dr") dr=atof(argv[++i]);
    else CmdError(arg.c_str());
  }
  
  config cfg(scfg);
  
  //set up inputs
  if (dt!=0) cfg.Setdt(dt);
  else dt=cfg.dt;

  double start=cfg.t;
  double finish=start+runtime;
  double rij=0;
  int bin=0;
  int num_bins=(int)(fabs(cfg.Lx)/dr/2)+1;
  vector<int> histSiSi(num_bins), histSiF(num_bins), histFF(num_bins);
  vector<double> gSiSi(num_bins), gSiF(num_bins), gFF(num_bins);
  for (i=0; i<num_bins; i++){
    histSiSi[i]=histSiF[i]=histFF[i]=0;
    gSiSi[i]=gSiF[i]=gFF[i]=0.0;
  }
  double lambda=1;

  vector<subcell>::iterator cells_it;
  set<subcell*>::iterator c;
  SPI i0,j0;
  while (cfg.t<finish+cfg.dt){
    if (tso) cfg.dtOptimize();
    if (cfg.u!=0) cfg.FirstVV();
    cfg.TimeStepInit();
    for (cells_it=cfg.cells.begin(); cells_it!=cfg.cells.end(); cells_it++)
      for (i0=cells_it->abegin(); i0!=cells_it->aend(); i0++)  
	for (c=cells_it->nbegin(); c!=cells_it->nend(); c++)
	  for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	    if (*i0 < *j0){
	      rij=((*i0)->R - (*j0)->R).minsqmag(cfg.Lx, cfg.Ly);
	      bin=(int)(sqrt(rij)/dr);
	      if (bin<num_bins && !(*i0)->is_fixed && !(*j0)->is_fixed){
		if ((*i0)->id + (*j0)->id==28)
		  histSiSi[bin]+=2;
		else if ((*i0)->id + (*j0)->id==23)
		  histSiF[bin]+=2;
		else
		  histFF[bin]+=2;
	      }
	      if (rij< ((*i0)->rc >? (*j0)->rc)){
		(*i0)->nlist.insert(*j0);
		(*j0)->nlist.insert(*i0);
	      }
	    }
    cfg.ForceEval();
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
  char outbuf[256];
  double C=4.0*PI*(double)(cfg.N-cfg.NFixed)
    /fabs(cfg.Lx*cfg.Ly*cfg.Lz)/3.0;
//    float sumhSiSi, sumhSiF, sumhFF;
//    sumhSiSi=sumhSiF=sumhFF=0;
//    for (i=0; i<num_bins; i++){
//      sumhSiSi+=histSiSi[i];
//      sumhSiF+=histSiF[i];
//      sumhFF+=histFF[i];
//    }
//    sumhSiSi*=dr; sumhSiF*=dr; sumhFF*=dr;
  for (i=0; i<num_bins; i++){
    double rlower=i*dr;
    double rupper=rlower+dr;
    double nideal=C*(pow(rupper,3)-pow(rlower,3));
    gSiSi[i]=(double)histSiSi[i]/((finish-start)/dt)/
      (double)(cfg.N-cfg.NFixed)/nideal;
    gSiF[i]=(double)histSiF[i]/((finish-start)/dt)/
      (double)(cfg.N-cfg.NFixed)/nideal;
    gFF[i]=(double)histFF[i]/((finish-start)/dt)/
      (double)(cfg.N-cfg.NFixed)/nideal;
    sprintf(outbuf,"%6.2f\t%e\t%e\t%e\n",i*dr, gSiSi[i], gSiF[i], gFF[i]);
    fout<<outbuf;
  }
  cerr<<"RDF complete. Output file: "<<rdf_file<<"."<<endl;
  return(0);
}
x
