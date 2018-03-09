//c_info.cpp
//member functions that return interesting things about the config.

#include "config.h"

//**************************************************************************
double config::MaxZ(){
  d1=-1000;
  for (i=begin; i<end; i++) d1 = i->R.z >? d1;
  return d1;
}
//**************************************************************************
double config::MinZ(){
  d1=1000;
  for (i=begin; i<end; i++) d1 = i->R.z <? d1;
  return d1;
}
//***************************************************************************

