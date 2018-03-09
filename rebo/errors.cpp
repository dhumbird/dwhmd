#include "errors.h"

void CmdError(const char * arg){
  cerr<<"error. "<<arg;
  cerr<<" is not recognized."<<endl;
  exit(0);
}
