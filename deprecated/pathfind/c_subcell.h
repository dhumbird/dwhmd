#ifndef SUBCELL_H
#define SUBCELL_H

//this is a class defining the subcells used in the linked-cell method.

class subcell{
 public:
  set<particle*> atom;
  set<subcell*> next;
  set<subcell*>::iterator nbegin(){return next.begin();}
  set<subcell*>::iterator nend(){return next.end();}
  SPI abegin(){return atom.begin();}
  SPI aend(){return atom.end();}
  void clear(){atom.clear();}
  void insert(particle* p){atom.insert(p); p->my_cell=this;}
  void erase(particle *p){atom.erase(p);}
  void ninsert(subcell* p){next.insert(p);}
  int size(){return atom.size();}
};


#endif

