#include "../config.h"
#include <list>
#include <map>

void HalBondPartner(int argc, char * argv[]){
  list<string> filelist;
  set<int> exclude;
  string arg;
  bool uniq=0;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) filelist.push_back(arg);
    else if (arg=="-exclude") exclude.insert(atoi(argv[++i]));
    else if (arg=="-nuniq") uniq=1;
  }
  if (uniq){
    float hist[20][2];
    for (int i=0; i<20; i++)
      for (int j=0; j<2; j++)
	hist[i][j]=0;
    map<int, list<int> > partners;
    {
      config cfg(filelist.back());
      for (VPI i=cfg.begin; i<cfg.end; i++)
	if (i->id==9 || i->id==17)
	  exclude.insert(i->ix);
    }
    list<string>::iterator f;
    for (f=filelist.begin(); f!=filelist.end(); f++){
      config cfg(*f);
      for (VPI i=cfg.begin; i<cfg.end; i++)
	if (i->id==9 || i->id==17)
	  if (!binary_search(exclude.begin(), exclude.end(), i->ix)){
	    int partner=-1;
	    double d1=1000; double d2=1000;
	    for (SPI s=i->nlist.begin(); s!=i->nlist.end(); s++){
	      d2=d1;
	      d1 = d1 <? ((*s)->R - i->R).minmag(cfg.Lx,cfg.Ly);
	      if (d2>d1) partner=(*s)->ix;
	    }
	    partners[i->ix].push_back(partner);
	  }
    }
    map<int, list<int> >::iterator m;
    list<int>::iterator v;
    for (m=partners.begin(); m!=partners.end(); m++){
      int bin=int(m->second.size()/50);
      hist[bin][0]++;
      m->second.unique();
      hist[bin][1]+=m->second.size();
    }
    for (int i=0; i<20; i++)
      if (hist[i][0]>0)
	cout<<i*50<<","<<hist[i][0]<<","<<hist[i][1]/hist[i][0]<<endl;
  }
  else{
    int count=0;
    list<string>::iterator f;
    for (f=filelist.begin(); f!=filelist.end(); f++){
      config cfg(*f);
      for (VPI i=cfg.begin; i<cfg.end; i++)
	if (i->id==9 || i->id==17)
	  if (!binary_search(exclude.begin(), exclude.end(), i->ix)){
	    int partner=-1;
	    double d1=1000; double d2=1000;
	    for (SPI s=i->nlist.begin(); s!=i->nlist.end(); s++){
	      d2=d1;
	      d1 = d1 <? ((*s)->R - i->R).minmag(cfg.Lx,cfg.Ly);
	      if (d2>d1) partner=(*s)->ix;
	    }
	    cout<<*f<<"\t"<<i->ix<<"\t"<<partner<<endl;
	  }
      count++;
    }
  }
}

