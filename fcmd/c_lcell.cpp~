//c_lcell.cpp
//config's member functions for the linked-cell method.
#include "config.h"
extern map<short,double> RC_SQ;
//*************************************************************************
void config::Partition(){
  cells.clear();
  double rc_max=3;
  int a=0;
  Lz_full=0;
  the_top=-1000; the_bottom=1000;
  for (VAI p=begin; p!=end; p++){
    the_top = mfmax(the_top, p->R.z);
    the_bottom = mfmin(the_bottom, p->R.z);
  }
  the_bottom -= 0.0001;
  the_top += 0.0001;
  Lz_full= mfmax(Lz, the_top - the_bottom);
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
  for (vscit=cells.begin(); vscit!=cells.end(); vscit++)
    vscit->clear();
  for (VAI p=begin; p!=end; p++) WhichCell(p->R)->insert(&(*p));
}
	     
