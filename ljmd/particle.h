//particle.h
//class listing for particle class

#ifndef PARTICLE_H
#define PARTICLE_H 1

class particle{

 public:
  vector F;
  vector R;
  vector V;
  short* nlist;
  double m;

  particle();
  ~particle();
  void listinit(int N);
  friend ostream& operator << (ostream& out, particle& p);
  friend istream& operator >> (istream& in, particle& p);

};

#endif
