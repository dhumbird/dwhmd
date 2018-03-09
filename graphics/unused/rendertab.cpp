#include "rendertab.h"

map<short, renderinfo> rendertab;
bool ren_is_def=0;

void RenderDef(){

  ren_is_def=1;  
 
  //silicon (14)
  rendertab[14].color.set(0,0,1);
  rendertab[14].rad=0.6;
  
  //argon (18)
  rendertab[18].color.set(0,1,1);
  rendertab[18].rad=0.6;
}
