#include "../config.h"
#include <list>
//***********************************************************************
void MSqD(int argc, char * argv[]){
  float res=1;
  string outfile;
  set<int> remain;
  VPI i;
  int mono=128;
  list<string>filelist;
  string arg;
  for (int j=2; j<argc; j++){
    arg=argv[j];
    if (sfind(arg, ".cfg")) filelist.push_back(arg);
    else if (arg=="-mono") mono=atoi(argv[++j]);
  }
  int runs=filelist.size();
  {
    config cfg(filelist.back());
    outfile=cfg.name+".msd";
    for (i=cfg.begin; i<cfg.end; i++)
      if (!i->is_fixed)
	remain.insert(i->ix);
  }

  float my_top;
  map<int, set<int> > atom_map;
  {
    config cfg(filelist.front());
    my_top=cfg.AvgTop();
    for (i=cfg.begin; i<cfg.end; i++)
      if (binary_search(remain.begin(), remain.end(), i->ix))
	atom_map[(int)((my_top-i->R.z)/res)].insert(i->ix);
  }
  
  
  fstream fout(outfile.c_str(), ios::out);
  list<string>::iterator f;
  map<int, set<int> >::iterator m;
  set<int>::iterator s;
  fout<<" ";
  for (m=atom_map.begin(); m!=atom_map.end(); m++)
    if (m->first>=0) fout<<","<<m->first*res;
  fout<<endl;
  int run=-1;
  for (f=filelist.begin(); f!=filelist.end(); f++){
    run++;
    if ((run*4)%mono==0){
      fout<<(float)run/(float)mono;
      config cfg(*f);
      for (m=atom_map.begin(); m!=atom_map.end(); m++){
	if (m->first>=0){
	  float avg=0;
	  for (s=m->second.begin(); s!=m->second.end(); s++){
	    i=cfg.atomix(*s);
	    svector P(i->xpass*cfg.Lx, i->ypass*cfg.Ly, 0);
	    avg+=(i->R+P-i->O).mag();
	  }
	  fout<<","<<avg/m->second.size();
	}
      }
      fout<<endl;
    }
  }
}
//*********************************************************************
void ZDiff(int argc, char * argv[]){
  string lastfile=argv[argc-1];
  string firstfile=argv[2];
  string outfile;
  set<int> remain;
  map<int, int>parts;
  map<int, int>::iterator mii;
  double z_bottom=0; double z_top=0;
  vector<particle>::iterator i;
  int nparts=10;

  { //which atoms are remaining?
    config cfg(lastfile);
    outfile=cfg.name+".zdiff";
    for (i=cfg.begin; i<cfg.end; i++)
      if (!i->is_fixed && i->id==cfg.mat)
	remain.insert(i->ix);
  }
  fstream fout(outfile.c_str(), ios::out);
  { //how to set up z-partitions?
    config cfg(firstfile);
    cfg.Rz_sort_asc();
    i=cfg.begin;
    while (i<cfg.end && i->is_fixed){
      z_bottom=i->R.z;
      i++;
    }
    //i should now be lowest non-fixed atom
    z_top=cfg.the_top;
    while (i<cfg.end){
      if (binary_search(remain.begin(), remain.end(), i->ix))
	parts[i->ix]=int((i->R.z-z_bottom)/(z_top-z_bottom)*nparts);
      i++;
    }
  }
  for (int j=2; j<argc; j++){
    fout<<j-2;
    config cfg((string)argv[j]);
    for (int k=0; k<nparts; k++){
      double d=0; int count=0;
      for (i=cfg.begin; i<cfg.end; i++){
	if (parts[i->ix]==k && 
	    binary_search(remain.begin(), remain.end(), i->ix)){
	  svector P(i->xpass*cfg.Lx, i->ypass*cfg.Ly, 0);
	  d+=(i->R+P-i->O).sqmag();
	  count++;
	}
      }
      if (count!=0) fout<<" "<<d/count;
      else fout<<" 0";
    }
    fout<<endl;  
  }
}
  
      
