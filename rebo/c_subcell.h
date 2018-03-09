#ifndef SUBCELL_H
#define SUBCELL_H

//this is a class defining the subcells used in the linked-cell method.

class subcell{
 public:
  set<atom*> atom_list;
  set<subcell*> next;
  set<subcell*>::iterator nbegin(){return next.begin();}
  set<subcell*>::iterator nend(){return next.end();}
  SAI abegin(){return atom_list.begin();}
  SAI aend(){return atom_list.end();}
  void clear(){atom_list.clear();}
  void insert(atom* p){atom_list.insert(p); p->my_cell=this;}
  void erase(atom *p){atom_list.erase(p);}
  void ninsert(subcell* p){next.insert(p);}
  int size(){return atom_list.size();}
};


#endif

