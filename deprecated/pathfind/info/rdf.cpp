#include "../config.h"
#include "../errors.h"

int RadDistFun(int argc, char * argv[]){

  string scfg="temp_000000.cfg";
  float runtime=10; float dr=0.1;
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
  vector<int> hist1(num_bins);
  vector<int> hist2(num_bins);
  vector<double> g1(num_bins);
  vector<double> g2(num_bins);
  for (i=0; i<num_bins; i++){
    hist1[i]=0;
    hist2[i]=0;
    g1[i]=0.0;
    g2[i]=0.0;
  }
  
  cerr<<"Calculating RDF for "<<cfg.name<<"."<<endl;
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
		  hist1[bin]+=2;
		else
		  hist2[bin]+=2;
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
  for (i=0; i<num_bins; i++){
    double rlower=i*dr;
    double rupper=rlower+dr;
    double nideal=C*(pow(rupper,3)-pow(rlower,3));
    g1[i]=(double)hist1[i]/((finish-start)/dt)/
      (double)(cfg.N-cfg.NFixed)/nideal;
    g2[i]=(double)hist2[i]/((finish-start)/dt)/
      (double)(cfg.N-cfg.NFixed)/nideal;
    sprintf(outbuf,"%6.2f\t%e\t%e\n",i*dr, g1[i], g2[i]);
    fout<<outbuf;
  }
  cerr<<"RDF complete. Output file: "<<rdf_file<<"."<<endl;
  return(0);
}

void InfoIndivRDF(int argc, char * argv[]){

  string scfg="temp_00000000.cfg";
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
  VPI I;
  double rcmax=0;
  for (I=cfg.begin; I!=cfg.end; I++)
    rcmax=rcmax >? I->rc;
  
  int num_bins=(int)(cfg.Lx/dr/2)+1;
  //int num_bins=(int)(sqrt(rcmax)/dr)+1;
  map<VPI, vector<float> > g;
  map<VPI, float> gmax, gsum;
  map<VPI, int> gmax_loc;
  for (I=cfg.begin; I!=cfg.end; I++)
    for (i=0; i<num_bins; i++){
      g[I].push_back(0);
      gmax[I]=gsum[I]=0;
      gmax_loc[I]=0;
    }
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
	      if ((*i0)->id==cfg.mat && (*j0)->id==cfg.mat){
		bin=(int)(sqrt(rij)/dr);
		if (bin<num_bins){
		  g[*i0][bin]+=1;
		  g[*j0][bin]+=1;
		}
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
  rdf_file+=".irdf";
  fout.open(rdf_file.c_str(), ios::out);
  fout.fill(' ');
  fout.precision(6);
  fout.setf(ios::fixed);
  double C=4.0*PI*(double)(cfg.N-cfg.NFixed)/fabs(cfg.Lx*cfg.Ly*cfg.Lz)/3.0;
  for (I=cfg.begin; I!=cfg.end; I++){
    fout<<I->ix;
    for (i=0; i<num_bins; i++){
      double rlower=i*dr;
      double rupper=rlower+dr;
      double nideal=C*(pow(rupper,3)-pow(rlower,3));
      g[I][i]/=((finish-start)/dt*(double)(cfg.N-cfg.NFixed)*nideal);
      if (g[I][i]>gmax[I]){
	gmax[I]=g[I][i];
	gmax_loc[I]=i;
	gsum[I]+=g[I][i]*nideal;
      }
      fout<<","<<g[I][i];
    }
    fout<<endl;
  }
  for (I=cfg.begin; I!=cfg.end; I++)
    if (!I->is_fixed)
      cout<<I->ix<<"\t"<<gmax[I]<<"\t"<<gmax_loc[I]<<"\t"<<gsum[I]<<endl;
  
}
