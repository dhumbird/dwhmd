#include <stdlib.h>
#include <iostream.h>

main(){

  float r;
  int i;

  srand(120);

  for (i=0; i<10; i++){
    r=rand()/float(RAND_MAX);
    cout<<r<<endl;
  }
};

