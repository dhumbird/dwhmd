class cluster{
 public:
  set<atom*> members;
  bool cl_is_fixed;
  cluster(){cl_is_fixed=0;}
  void insert(atom* p){members.insert(p); p->p_clust=this;}
  SAI begin(){return members.begin();}
  SAI end(){return members.end();}
  int size(){return members.size();}
};
