#
# plot potential energy vs r for the following configuration
#
#  (Si)---(F)- - -(F)
#    1     2       3
#
#  where r(12) is fixed at 1.6 A, and atoms are collinear.
#  
set term post eps enhanced "Helvetica" 22
set encoding iso_8859_1
set out "fsif.eps"
set size 0.7,1.0
rh(r)=r/2.0951	# conversion of A to sigmas
ep=2.1678	# conversion of epsilons to eV
r12=1.6		# set Si--F bond length
Asw=21.23414138
Bsw=0.5695476433
csw(r)=exp(1.3/(rh(r)-1.8))
vsf(r)=rh(r)<1.8?ep*Asw*(Bsw*(rh(r))**(-3.)-(rh(r))**(-2.))*csw(r):0
cff(r)=exp(0.0818182/(rh(r)-2.086182))
Aff=0.52276
Bff=0.112771
vff(r)=rh(r)<2.086182?ep*Aff*(Bff*((rh(r))**(-8.0))-((rh(r))**(-4.0)))*cff(r):0

r12=1.6;
r13=1.6;
v12=vsf(r12)
v13=vsf(r13)
M_PI  = 3.14159265358979323846
rad(th)=M_PI/180.0*(th)
rff(th)=sqrt(r12*r12+r13*r13-2*r12*r13*cos(rad(th)))
v23(th)=vff(rff(th))

h213(th)=(rh(r12)<1.8&&rh(r13)<1.8)?\
    ep*(24*(cos(rad(th))-cos(rad(103)))**2-3.2)*\
    exp(1./(rh(r12)-1.8)+1./(rh(r13)-1.8)):0

h123=rh(r12)<1.8?ep*3.5*exp(1.0/(rh(r12)-1.8)+1.0/(rh(r13)-1.8)):0
h132=rh(r13)<1.8?ep*3.5*exp(1.0/(rh(r12)-1.8)+1.0/(rh(r13)-1.8)):0

print "h123 ", h123
print "h132 ", h132
print rh(r12)

v123(th)=v23(th)+h213(th)+h123+h132+v12+v13

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
set xlabel "{{/Symbol q}_{F--Si--F}/ ^{o}" 0.0,0.25
set key spacing 1.2
#set xtics 90, 10, 140
plot [th=90:180] [-12:-5] v123(th) t"SW" w l, "fsif-t.dat" t"T" w l
#0 not w l -1
#plot [r=0:10] rh(r12+r) t"vff" w l
!awk -f bbox.awk fsif.eps > tmp
!mv tmp fsif.eps
