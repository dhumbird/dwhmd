//c_lcell.cpp
//config's member functions for the linked-cell method.
#include "config.h"
//*************************************************************************
void config::Partition(){
  double rc_max=0;
  sLz=0;
  if(Lz > 0){
    for (p=begin; p!=end; p++) rc_max = rc_max >? p->rc;
    sLz=Lz;
    sLz2=Lz/2;
  }
  else{
    float zmax=-1000; float zmin=1000;
    for (p=begin; p!=end; p++){
      rc_max = rc_max >? p->rc;
      zmax = zmax >? p->R.z;
      zmin = zmin <? p->R.z;
    }
    zmin-=0.0001;
    sLz=zmax-zmin;
    sLz2=zmax;
  }
  rc_max=sqrt(rc_max);
  Nx=(int)(Lx/rc_max);
  Ny=(int)(Ly/rc_max);
  Nz=(int)(sLz/rc_max);
  Lcx=Lx/Nx;
  Lcy=Ly/Ny;
  Lcz=sLz/Nz;
  Nc=Nx*Ny*Nz;
  cells.resize(Nc);
  for (a=0; a<Nc; a++){
    cells[a].Lx=Lx;
    cells[a].Ly=Ly;
    cells[a].Lz=Lz;
  }
  for (x=0; x<Nx; x++)
    for (y=0; y<Ny; y++)
      for (z=0; z<Nz; z++){
	a=x+Nx*y+Nx*Ny*z;
	//same z-plane
	cells[a].ninsert(&cells[x + Nx*y + Nx*Ny*z]);
	cells[a].ninsert(&cells[x + Nx*((y+1)%Ny) + Nx*Ny*z]);
	cells[a].ninsert(&cells[x + Nx*((y-1+Ny)%Ny) + Nx*Ny*z]);
	cells[a].ninsert(&cells[(x+1)%Nx + Nx*y + Nx*Ny*z]);
	cells[a].ninsert(&cells[(x+1)%Nx + Nx*((y+1)%Ny) + Nx*Ny*z]);
	cells[a].ninsert(&cells[(x+1)%Nx + Nx*((y-1+Ny)%Ny) + Nx*Ny*z]);
	cells[a].ninsert(&cells[(x-1+Nx)%Nx + Nx*y + Nx*Ny*z]);
	cells[a].ninsert(&cells[(x-1+Nx)%Nx + Nx*((y+1)%Ny) + Nx*Ny*z]);
	cells[a].ninsert(&cells[(x-1+Nx)%Nx + Nx*((y-1+Ny)%Ny) + Nx*Ny*z]);
	if (Lz > 0 || z > 0){
	  //plane above
	  cells[a].ninsert(&cells[x+Nx*y+Nx*Ny*((z-1+Nz)%Nz)]);
	  cells[a].ninsert(&cells[x+Nx*((y+1)%Ny)+Nx*Ny*((z-1+Nz)%Nz)]);
	  cells[a].ninsert(&cells[x+Nx*((y-1+Ny)%Ny)+Nx*Ny*((z-1+Nz)%Nz)]);
	  cells[a].ninsert(&cells[(x+1)%Nx+Nx*y+Nx*Ny*((z-1+Nz)%Nz)]);
	  cells[a].ninsert(&cells[(x+1)%Nx+Nx*((y+1)%Ny)+Nx*Ny*((z-1+Nz)%Nz)]);
	  cells[a].ninsert
	    (&cells[(x+1)%Nx+Nx*((y-1+Ny)%Ny)+Nx*Ny*((z-1+Nz)%Nz)]);
	  cells[a].ninsert(&cells[(x-1+Nx)%Nx+Nx*y+Nx*Ny*((z-1+Nz)%Nz)]);
	  cells[a].ninsert
	    (&cells[(x-1+Nx)%Nx+Nx*((y+1)%Ny)+Nx*Ny*((z-1+Nz)%Nz)]);
	  cells[a].ninsert
	    (&cells[(x-1+Nx)%Nx+Nx*((y-1+Ny)%Ny)+Nx*Ny*((z-1+Nz)%Nz)]);
	}
	if (Lz > 0 || z < Nz-1){
	  //plane below
	  cells[a].ninsert(&cells[x + Nx*y + Nx*Ny*((z+1)%Nz)]);
	  cells[a].ninsert(&cells[x + Nx*((y+1)%Ny) + Nx*Ny*((z+1)%Nz)]);
	  cells[a].ninsert(&cells[x + Nx*((y-1+Ny)%Ny) + Nx*Ny*((z+1)%Nz)]);
	  cells[a].ninsert(&cells[(x+1)%Nx + Nx*y + Nx*Ny*((z+1)%Nz)]);
	  cells[a].ninsert
	    (&cells[(x+1)%Nx + Nx*((y+1)%Ny) + Nx*Ny*((z+1)%Nz)]);
	  cells[a].ninsert
	    (&cells[(x+1)%Nx + Nx*((y-1+Ny)%Ny) + Nx*Ny*((z+1)%Nz)]);
	  cells[a].ninsert(&cells[(x-1+Nx)%Nx + Nx*y + Nx*Ny*((z+1)%Nz)]);
	  cells[a].ninsert
	    (&cells[(x-1+Nx)%Nx + Nx*((y+1)%Ny) + Nx*Ny*((z+1)%Nz)]);
	  cells[a].ninsert
	    (&cells[(x-1+Nx)%Nx+Nx*((y-1+Ny)%Ny)+Nx*Ny*((z+1)%Nz)]);
	}
      }
  Allocate();
  cellset1.clear(); cellset2.clear();
  cells_it=cells.begin(); 
  while (cells_it!=cells.end()){
    cellset1.insert(&(*(cells_it++)));
    if (cells_it!=cells.end()) cellset2.insert(&(*(cells_it++)));
  }
}
//************************************************************************
void config::Allocate(){
  for (cells_it=cells.begin(); cells_it!=cells.end(); cells_it++)
    cells_it->clear();
  for (p=begin; p!=end; p++)
    cells[(int)((p->R.x + Lx/2)/Lcx)%Nx +
	 Nx*((int)(fabs((p->R.y - Ly/2)/Lcy))%Ny) +
	 Nx*Ny*((int)(fabs((p->R.z - sLz2)/Lcz))%Nz)].insert(p);
}
