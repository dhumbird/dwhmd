//c_lcell.cpp
//config's member functions for the linked-cell method.
#include "config.h"
//*************************************************************************
void config::Partition(){
  cells.clear();
  double rc_max=0;
  int a=0;
  Lz_full=0;
  float zmax=-1000; float zmin=1000;
  for (VPI p=begin; p!=end; p++){
    rc_max = rc_max >? p->rc;
    zmax = zmax >? p->R.z;
    zmin = zmin <? p->R.z;
  }
  zmin-=0.0001;
  the_top = zmax;
  the_top+=0.0001;
  Lz_full= Lz >? the_top-zmin;
  rc_max=sqrt(rc_max);
  Nx=(int)(Lx/rc_max);
  Ny=(int)(Ly/rc_max);
  Nz=(int)(Lz_full/rc_max);
  Lcx=Lx/Nx;
  Lcy=Ly/Ny;
  Lcz=Lz_full/Nz;
  Nc=Nx*Ny*Nz;
  cells.resize(Nc);
  int xm,xp,ym,yp,zm,zp;
  for (int x=0; x<Nx; x++)
    for (int y=0; y<Ny; y++)
      for (int z=0; z<Nz; z++){
	a=x+Nx*y+Nx*Ny*z;
	xp=(x+1)%Nx; xm=(x-1+Nx)%Nx;
	yp=(y+1)%Ny; ym=(y-1+Ny)%Ny;
	zp=(z+1)%Nz; zm=(z-1+Nz)%Nz;
	
	//same z-plane
	cells[a].ninsert(&cells[x + Nx*(y + Ny*z)]);
	cells[a].ninsert(&cells[x + Nx*(yp + Ny*z)]);
	cells[a].ninsert(&cells[x + Nx*(ym + Ny*z)]);
	cells[a].ninsert(&cells[xp + Nx*(y + Ny*z)]);
	cells[a].ninsert(&cells[xp + Nx*(yp + Ny*z)]);
	cells[a].ninsert(&cells[xp + Nx*(ym + Ny*z)]);
	cells[a].ninsert(&cells[xm + Nx*(y + Ny*z)]);
	cells[a].ninsert(&cells[xm + Nx*(yp + Ny*z)]);
	cells[a].ninsert(&cells[xm + Nx*(ym + Ny*z)]);
	if (z > 0){
	  //plane above
	  cells[a].ninsert(&cells[x + Nx*(y + Ny*zm)]);
	  cells[a].ninsert(&cells[x + Nx*(yp + Ny*zm)]);
	  cells[a].ninsert(&cells[x + Nx*(ym + Ny*zm)]);
	  cells[a].ninsert(&cells[xp + Nx*(y + Ny*zm)]);
	  cells[a].ninsert(&cells[xp + Nx*(yp + Ny*zm)]);
	  cells[a].ninsert(&cells[xp + Nx*(ym + Ny*zm)]);
	  cells[a].ninsert(&cells[xm + Nx*(y + Ny*zm)]);
	  cells[a].ninsert(&cells[xm + Nx*(yp + Ny*zm)]);
	  cells[a].ninsert(&cells[xm + Nx*(ym + Ny*zm)]);
	}
	if (z < Nz-1){
	  //plane below
	  cells[a].ninsert(&cells[x + Nx*(y + Ny*zp)]);
	  cells[a].ninsert(&cells[x + Nx*(yp + Ny*zp)]);
	  cells[a].ninsert(&cells[x + Nx*(ym + Ny*zp)]);
	  cells[a].ninsert(&cells[xp + Nx*(y + Ny*zp)]);
	  cells[a].ninsert(&cells[xp + Nx*(yp + Ny*zp)]);
	  cells[a].ninsert(&cells[xp + Nx*(ym + Ny*zp)]);
	  cells[a].ninsert(&cells[xm + Nx*(y + Ny*zp)]);
	  cells[a].ninsert(&cells[xm + Nx*(yp + Ny*zp)]);
	  cells[a].ninsert(&cells[xm + Nx*(ym + Ny*zp)]);
	}
      }
  //allocate particles
  for (cells_it=cells.begin(); cells_it!=cells.end(); cells_it++)
    cells_it->clear();
  for (VPI p=begin; p!=end; p++) WhichCell(p->R)->insert(&(*p));
}
	     
