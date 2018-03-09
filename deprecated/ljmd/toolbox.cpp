//toolbox.cpp
//general-use functions

//****************************************************************
void ReNeighbor(int N, particle p[], double rl){
  int i,j;
  double rijsq;
  
  //re-neighbor
  for (i=0; i<N-1; i++){
    for (j=i+1; j<N; j++){
      p[i].nlist[j]=0;
      if (i!=j){
	rijsq=(p[i].R-p[j].R).sqmag();
	if (rijsq<rl) p[i].nlist[j]=1;
      }
    }
  }
}
