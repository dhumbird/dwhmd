//rendertab.h
//rendering properties

#ifndef RENDERTAB_H
#define RENDERTAB_H

#include <map>

class RGB{
 public:
  float r, g, b;
  RGB(){r=0; g=0; b=0;}
  void clear(){r=0; g=0; b=0;}
  void set(float f1, float f2, float f3){r=f1; g=f2; b=f3;}
  void operator = (RGB c){set(c.r, c.g, c.b);}
};

typedef struct renderinfo{
  RGB color;
  float rad;
};

void RenderDef();


#endif
