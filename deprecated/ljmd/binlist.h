//binlist.h
//class listing for binlist class

#ifndef BINLIST_H
#define BINLIST_H 1

class binlist{

 public:
  int bins;
  int* histo;
  double min;
  double max;
  double inc;

  binlist(double MIN, double MAX, double INC);
  ~binlist();
  void toss(double r);
  void display(double c);

};

#endif
