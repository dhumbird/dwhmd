//rgb.h
//a color class

#ifndef RGB_H
#define RGB_H

#include <cmath>
#include <iostream>
#include "../consts.h"

class rgb{
 public:
  float r, g, b;
  rgb(){r=0; g=0; b=0;}
  rgb(float x1, float x2, float x3){r=x1; g=x2; b=x3;}
  void clear(){r=0; g=0; b=0;}
  void set(float f1, float f2, float f3){r=f1; g=f2; b=f3;}
  void operator = (rgb c){set(c.r, c.g, c.b);}
  void Set01(float);
  friend ostream& operator << (ostream&, const rgb&);
};

#endif
