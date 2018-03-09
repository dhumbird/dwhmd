class cluster{
 public:
  set<atom*> members;
  bool cl_is_fixed;
  double fmax;
  cluster(){cl_is_fixed=0; fmax=0;}
  void insert(atom* p){members.insert(p); p->p_clust=this;}
  SAI begin(){return members.begin();}
  SAI end(){return members.end();}
  int size(){return members.size();}

};
