set term postscript eps enhanced "Helvetica" 22
set encoding iso_8859_1
set size 0.7,1.0
set out "sif.eps"
rh(r)=r/2.0951
Asw=21.23414138
Bsw=0.5695476433
csw(r)=exp(1.3/(rh(r)-1.8))
vsw(r)=rh(r)<1.8?2.1678*Asw*(Bsw*((rh(r))**(-3.))-((rh(r))**(-2.)))*csw(r):0
A=37412.28
B=925.846
l=5.4875
m=2.7437
r1=1.83922
r2=2.13922
a(r)=pi*((r-(r2+r1)/2)/(r2-r1))
f(r)=r<r1?1:r<r2?0.5-9./16.*sin(a(r))-1./16.*sin(3*a(r)):0
v(r)=f(r)*(A*exp(-l*r)-B*exp(-m*r))
set ylabel "V/eV" 0.5,0.0
set xlabel "r_{Si--F}/ {\305}" 0.0,0.25
set key spacing 1.2
plot [r=1:3.5] [-6:16] vsw(r) t"SW" w l 1, v(r) t"T" w l 2, 0 not w l -1 
