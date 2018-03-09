class cluster{
 public:
  set<particle*> members;
  bool cl_is_fixed;
  cluster(){cl_is_fixed=0;}
  void insert(particle* p){members.insert(p); p->p_clust=this;}
  SPI begin(){return members.begin();}
  SPI end(){return members.end();}
  int size(){return members.size();}
};
