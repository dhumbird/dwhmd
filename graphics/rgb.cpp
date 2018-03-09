#include "rgb.h"

void rgb::Set01(float f){
  f=f*TWO_THIRDS+ONE_THIRD;
  if(f>=0.0 && f<ONE_SIXTH){
    r=1; g=0; b=6*f;}
  if(f>=ONE_SIXTH && f<ONE_THIRD){
    r=1-6*(f-ONE_SIXTH); g=0; b=1;}
  if(f>=ONE_THIRD && f<ONE_HALF){
    r=0; g=6*(f-ONE_THIRD); b=1;}
  if(f>=ONE_HALF && f<TWO_THIRDS){
    r=0; g=1; b=1-6*(f-ONE_HALF);}
  if(f>=TWO_THIRDS && f<FIVE_SIXTHS){
    r=6*(f-TWO_THIRDS); g=1; b=0;}
  if(f>=FIVE_SIXTHS){
    r=1; g=1-6*(f-FIVE_SIXTHS); b=0;}
  if(f>=1){
    r=1; g=0; b=0;}
}
//*************************************************************************
ostream& operator << (ostream& out, const rgb& v){
  out<<v.r<<" "<<v.g<<" "<<v.b<<" ";
  return out;
}
//*************************************************************************
istream& operator >> (istream& in, rgb& v){
  in>>v.r>>v.g>>v.b;
  return in;
}
