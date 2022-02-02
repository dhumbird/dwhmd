#include "config.h"

//*******************************************************************
void config::Load(string cfgfile)
{
  if (FileExists(cfgfile))
  {
    int critical=0;
    int commands=0;
    ifstream fin(cfgfile.c_str(), ios::in);
    string s="#";
    string arg;
    bool kill_me=0;
    while (!kill_me && !fin.eof())
    {
      getline(fin,s);
      if (sfind(s,"#"))
      {
        istringstream ist(s);
        ist>>arg;
        if (arg=="#N")
        {
          ist>>N; critical++;
        }
      	else if (sfind(arg, "ixed"))
        {
          ist>>NFixed; critical++;
        }
        else if (sfind(arg, "max"))
        {
          ist>>Nmax; critical++;
        }
      	else if (sfind(arg, "emperature"))
        {
          ist>>T; critical++;
        }
      	else if (sfind(arg, "imensions"))
        {
          ist>>Lx>>Ly>>Lz; critical++;
        }
        else if (sfind(arg, "ime"))
        {
          ist>>t; critical++;
        }
        else if (sfind(arg, "#-")) kill_me=1;
        else
        {
          cerr<<"I don't understand the argument "<<arg<<" in "<<cfgfile<<endl;
          exit(1);
        }
      }
    }
    if (critical < 6)
    {
      cerr<<"Hmm, something seems to be missing from "<<cfgfile<<endl;
      exit(1);
    }
    resize(N);
    for (VAI i=begin; i<end; i++)
    { 
      fin>>*i;
      if (i->id == 18)
      {
        if (inert == -1) inert=i->ix;
        else
        {
          cerr<<"Dude, there is more than one inert in this cfg!\n";
          exit(1);
        }
      }
    }
    fin.close();
    SetProps();
    name=cfgfile.substr(0, cfgfile.find("_"));
  }
  else
  {
    cerr<<"Can't open specified cfg file: "<<cfgfile<<endl;
    exit(1);
  } 
}
//************************************************************************
void config::Dump(){
  string out_cfg=name+"_"+time2string(t)+".cfg";
  fstream fout(out_cfg.c_str(), ios::out);
  fout.fill(' '); fout.precision(12); fout.setf(ios::fixed);
  fout<<"%"<<DateTime()<<" "<<out_cfg<<endl;
  fout<<"%file format revision 05232002"<<endl;
  fout<<"#N "<<N<<endl;
  fout<<"#Nfixed "<<NFixed<<endl;
  fout<<"#Nmax "<<Nmax<<endl;
  fout<<"#temperature "<<T<<endl;
  fout<<"#dimensions "<<Lx<<" "<<Ly<<" "<<Lz<<endl;
  fout<<"#time "<<t<<endl;
  fout<<"#----------"<<endl;
  for (VAI i=begin; i<end; i++) fout<<*i<<endl;
  cerr<<"* Cfg saved. File: "<<out_cfg<<endl;
}

void config::DumpFlist(){
  TimeStepInit();
  ReNeighbor();
  ForceEval();
  printf("%6s %10s %10s %10s %10s\n","index","Ek","F.x","F.y","F.z");
  for (VAI i=begin; i<end; i++)
  {
    printf("%6d %10.5f %10.5f %10.5f %10.5f\n",i->ix,i->Ek(),i->F.x,i->F.y,i->F.z);
    //cout<<i->ix<<" "<<i->Ek()<<" "<<i->F.x<<" "<<i->F.y<<" "<<i->F.z<<endl;
  } 

}