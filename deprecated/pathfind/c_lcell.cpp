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
  for (p=begin; p!=end; p++){
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
  for (int x=0; x<Nx; x++)
    for (int y=0; y<Ny; y++)
      for (int z=0; z<Nz; z++){
	a=x+Nx*y+Nx*Ny*z;
	//same z-plane
	cells[a].ninsert(&cells[x + Nx*(y + Ny*z)]);
	cells[a].ninsert(&cells[x + Nx*(((y+1)%Ny) + Ny*z)]);
	cells[a].ninsert(&cells[x + Nx*(((y-1+Ny)%Ny) + Ny*z)]);
	cells[a].ninsert(&cells[(x+1)%Nx + Nx*(y + Ny*z)]);
	cells[a].ninsert(&cells[(x+1)%Nx + Nx*(((y+1)%Ny) + Ny*z)]);
	cells[a].ninsert(&cells[(x+1)%Nx + Nx*(((y-1+Ny)%Ny) + Ny*z)]);
	cells[a].ninsert(&cells[(x-1+Nx)%Nx + Nx*(y + Ny*z)]);
	cells[a].ninsert(&cells[(x-1+Nx)%Nx + Nx*(((y+1)%Ny) + Ny*z)]);
	cells[a].ninsert(&cells[(x-1+Nx)%Nx + Nx*(((y-1+Ny)%Ny) + Ny*z)]);
	if (z > 0){
	  //plane above
	  cells[a].ninsert(&cells[x+Nx*(y+Ny*((z-1+Nz)%Nz))]);
	  cells[a].ninsert(&cells[x+Nx*(((y+1)%Ny)+Ny*((z-1+Nz)%Nz))]);
	  cells[a].ninsert(&cells[x+Nx*(((y-1+Ny)%Ny)+Ny*((z-1+Nz)%Nz))]);
	  cells[a].ninsert(&cells[(x+1)%Nx+Nx*(y+Ny*((z-1+Nz)%Nz))]);
	  cells[a].ninsert(&cells[(x+1)%Nx+Nx*(((y+1)%Ny)+Ny*((z-1+Nz)%Nz))]);
	  cells[a].ninsert
	    (&cells[(x+1)%Nx+Nx*(((y-1+Ny)%Ny)+Ny*((z-1+Nz)%Nz))]);
	  cells[a].ninsert(&cells[(x-1+Nx)%Nx+Nx*(y+Ny*((z-1+Nz)%Nz))]);
	  cells[a].ninsert
	    (&cells[(x-1+Nx)%Nx+Nx*(((y+1)%Ny)+Ny*((z-1+Nz)%Nz))]);
	  cells[a].ninsert
	    (&cells[(x-1+Nx)%Nx+Nx*(((y-1+Ny)%Ny)+Ny*((z-1+Nz)%Nz))]);
	}
	if (z < Nz-1){
	  //plane below
	  cells[a].ninsert(&cells[x + Nx*(y + Ny*((z+1)%Nz))]);
	  cells[a].ninsert(&cells[x + Nx*(((y+1)%Ny) + Ny*((z+1)%Nz))]);
	  cells[a].ninsert(&cells[x + Nx*(((y-1+Ny)%Ny) + Ny*((z+1)%Nz))]);
	  cells[a].ninsert(&cells[(x+1)%Nx + Nx*(y + Ny*((z+1)%Nz))]);
	  cells[a].ninsert
	    (&cells[(x+1)%Nx + Nx*(((y+1)%Ny) + Ny*((z+1)%Nz))]);
	  cells[a].ninsert
	    (&cells[(x+1)%Nx + Nx*(((y-1+Ny)%Ny) + Ny*((z+1)%Nz))]);
	  cells[a].ninsert(&cells[(x-1+Nx)%Nx + Nx*(y + Ny*((z+1)%Nz))]);
	  cells[a].ninsert
	    (&cells[(x-1+Nx)%Nx + Nx*(((y+1)%Ny) + Ny*((z+1)%Nz))]);
	  cells[a].ninsert
	    (&cells[(x-1+Nx)%Nx+Nx*(((y-1+Ny)%Ny)+Ny*((z+1)%Nz))]);
	}
      }
  //allocate particles
  for (cells_it=cells.begin(); cells_it!=cells.end(); cells_it++)
    cells_it->clear();
  for (p=begin; p!=end; p++) WhichCell(p->R)->insert(p);
}
	     
