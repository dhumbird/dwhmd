#
# plot potential energy vs r for the following configuration
#
#  (Si)---(F)- - -(F)
#    1     2       3
#
#  where r(12) is fixed at 1.6 A, and atoms are collinear.
#  the variable "r" denotes r23.
#  
#set term post eps enhanced "Helvetica" 22
set encoding iso_8859_1
#set out "siff.eps"
#set size 0.7,1.0

sig=2.0951
ep=2.1678	# conversion of epsilons to eV
r12=1.6/sig	# set Si--F bond length to 1.6 Angstroms

r13(r)=r12+r    # r13 is r12 + r23

Asf=21.23414138
Bsf=0.5695476433
csf(r)=r<1.8?exp(1.3/(r-1.8)):0
vsf(r)=Asf*(Bsf*(r**(-3.0))-(r**(-2.0)))*csf(r)

v12=vsf(r12)

v13(r)=vsf(r13(r))

cff(r)=r<2.086182?exp(0.579495/(r-2.086182)):0
Aff=0.52276
Bff=0.112771
vff(r)=Aff*(Bff*(r**(-8.0))-(r**(-4.0)))*cff(r)
v23(r)=vff(r)

h213(r)=(r12<1.8&&r13(r)<1.8)?\
   ((24*(1.22495**2))-3.2)*exp(1.0/(r12-1.8)+1.0/(r13(r)-1.8)):0
h123(r)=(r<1.8&&r12<1.8)?\
    3.5*exp(1.0/(r12-1.8)+1.0/(r-1.8)):0
h132(r)=(r13(r)<1.8&&r<1.8)?\
    3.5*exp(1.0/(r13(r)-1.8)+1.0/(r-1.8)):0

v123(r)=v23(r)+v13(r)+h213(r)+h123(r)+h132(r)

set ylabel "V/eV" 0.5,0.0
set xlabel "r_{SiF--F}/ {\305}" 0.0,0.25
set key spacing 1.2

#plot [r=0:5] [-6:4] vsf(r),  vff(r)
plot [r=0:4] [-2.5:2.5] ep*v123(r/sig) t"sum" w l -1, \
ep*v23(r/sig), ep*v13(r/sig), ep*h213(r/sig), ep*h123(r/sig), ep*h132(r/sig)
#plot [r=0:4] [-0.5:2.5] v123(r) t"SW","siff-t.dat" u 1:($2+5.723) t"T" w l,0 not w l -1
#plot [r=0:4] [-6.3:-3.3]  "siff-t.dat" u 1:2 t"T" w p
