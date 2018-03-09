//vector.h
//class listing for vector class

#ifndef VECTOR_H
#define VECTOR_H 1

#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>


class vector{
 private:
  double x;
  double y;
  double z;

 public:
  vector(); //generic constructor
  vector(double X, double Y, double Z); //specific constructor
  ~vector(); //destructor
  void clear();
  double sqmag()const;
  double mag()const;
  void minimg(const double L);
  void minimg(double Lx, double Ly, double Lz);
  void set(double X, double Y, double Z);
  double X()const;
  double Y()const;
  double Z()const;
  void operator = (const vector& v);
  void operator *= (const double d);
  void operator /= (const double d);
  void operator += (const vector& v);
  void operator -= (const vector& v);
  void operator += (const double d);
  void operator -= (const double d);
  friend double operator ^ (const vector& v1, const vector& v2);  //dot product
  friend vector operator * (const vector& v1, const vector& v2);  //cross prod
  friend vector operator * (const vector& v, const double d);
  friend vector operator * (const double d, const vector& v);
  friend vector operator / (const vector& v, const double d);
  friend vector operator + (const vector& v1, const vector& v2);
  friend vector operator - (const vector& v1, const vector& v2);
  friend ostream& operator << (ostream& out, const vector& v); //output
  friend istream& operator >> (istream& in, vector& v); //input
};

#endif
