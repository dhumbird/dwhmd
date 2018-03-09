#include "config.h"

double config::AvgTop(){
  //calculates the average height of the cell by forming a grid over the top,
  //finding the highest atom under each grid point, then averaging over z.

  int nx=(int)(Lx)+1;
  int ny=(int)(Ly)+1;
  double box[nx][ny];
  double my_top=0;
  for (int m=0; m<nx; m++) for (int n=0; n<ny; n++) box[m][n]=the_top+5;
  svector R;
  bool flag;
  for (int m=0; m<nx; m++){
    for (int n=0; n<ny; n++){
      flag=1;
      R.x=m-Lx/2;
      R.y=Ly/2-n;
      while(flag){
	box[m][n]-=1;
	R.z=box[m][n];
	int x=(int)((R.x + Lx/2)/Lcx)%Nx;
	int y=(int)(fabs((R.y - Ly/2)/Lcy))%Ny;
	int z=(int)(fabs((R.z - the_top)/Lcz))%Nz;
	int a=x+Nx*y+Nx*Ny*z;
	for (c=cells[a].nbegin(); c!=cells[a].nend(); c++)
	  for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	    if (flag)
	      if ((R - (*j0)->R).minsqmag(Lx, Ly) > (*j0)->rc)
		flag=1;
	      else{
		flag=0;
		box[m][n]=(*j0)->R.z;
		break;
	      }
      }
    }
  }
  for (int m=0; m<nx; m++) for (int n=0; n<ny; n++) my_top+=box[m][n];
  return my_top/(nx*ny);
}
