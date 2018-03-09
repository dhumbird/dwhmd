#include "mathfun.h"

void Gpoly(const double cosO, int id, double* G, double* Gprime, double* gam,
	   double* gamprime){
  double cosO2=cosO*cosO;
  double cosO3=cosO2*cosO;
  double cosO4=cosO3*cosO;
  double cosO5=cosO4*cosO;
  if (id==1){
    if (cosO >= -0.5){
      *G=
	19.065031149937783
	+2.017732531534021*cosO
	-2.566444502991983*cosO2
	+3.291353893907436*cosO3
	-2.653536801884563*cosO4
	+0.837650930130006*cosO5;
      *Gprime=
	2.017732531534021
	-2*2.566444502991983*cosO
	+3*3.291353893907436*cosO2
	-4*2.653536801884563*cosO3
	+5*0.837650930130006*cosO4;
    }
    else if (cosO >= -FIVE_SIXTHS){
      *G=
	16.956325544514659
	-21.059084522755980*cosO
	-102.394184748124742*cosO2
	-210.527926707779059*cosO3
	-229.759473570467513*cosO4
	-94.968528666251945*cosO5;
      *Gprime=
	-21.059084522755980
	-2*102.394184748124742*cosO
	-3*210.527926707779059*cosO2
	-4*229.759473570467513*cosO3
	-5*94.968528666251945*cosO4;
    }
    else{
      *G=
	270.467795364007301
	+1549.701314596994564*cosO
	+3781.927258631323866*cosO2
	+4582.337619544424228*cosO3
	+2721.538161662818368*cosO4
	+630.658598136730774*cosO5;
      *Gprime=
	1549.701314596994564
	+2*3781.927258631323866*cosO
	+3*4582.337619544424228*cosO2
	+4*2721.538161662818368*cosO3
	+5*630.658598136730774*cosO4;
    }
    *gam=*gamprime=0;
  }
  else{
    if (cosO >= -ONE_THIRD){
      *G=
	0.3754490870000
	+1.407252749388*cosO
	+2.255103926323*cosO2
	+2.028902219952*cosO3
	+1.426981217906*cosO4
	+0.5063107994308*cosO5;
      *Gprime=
	1.407252749388
	+2*2.255103926323*cosO
	+3*2.028902219952*cosO2
	+4*1.426981217906*cosO3
	+5*0.5063107994308*cosO4;
      *gam=
	0.2718558000000
	+0.4892727456293*cosO
	-0.4328199017473*cosO2
	-0.5616795197048*cosO3
	+1.270874966906*cosO4
	-0.03750409108350*cosO5;
      *gamprime=
	0.4892727456293
	-2*0.4328199017473*cosO
	-3*0.5616795197048*cosO2
	+4*1.270874966906*cosO3
	-5*0.03750409108350*cosO4;
    }
    else if (cosO >= -0.5){
      *G=
	0.6900668660000
	+5.460691360000*cosO
	+23.01345680000*cosO2
	+54.91519344000*cosO3
	+68.62037040000*cosO4
	+34.70897779200*cosO5;
      *Gprime=
	5.460691360000
	+2*23.01345680000*cosO
	+3*54.91519344000*cosO2
	+4*68.62037040000*cosO3
	+5*34.70897779200*cosO4;
      *gam=*gamprime=0;
    }
    else{
      *G=
	0.2817216000000
	+1.062912000000*cosO
	+2.136736000000*cosO2
	+2.533952000000*cosO3
	+1.554736000000*cosO4
	+0.3863296000000*cosO5;
      *Gprime=
	1.062912000000
	+2*2.136736000000*cosO
	+3*2.533952000000*cosO2
	+4*1.554736000000*cosO3
	+5*0.3863296000000*cosO4;
      *gam=*gamprime=0;
    }
  }
}

void PolySwitch(double t, double* S, double* Sprime){
  if (t <= 0){
    *S=1;
    *Sprime=0;
  }
  else if (t <= 1){
    *S = 1; *Sprime=0;
    double tt=t*t; //t^2
    *Sprime -= 30*tt;
    tt*=t; //t^3
    *S -= 10*tt;
    *Sprime += 60*tt;
    tt*=t; //t^4
    *S += 15*tt;
    *Sprime -=30*tt;
    *S-=6*tt*t;
  }
  else{
    *S=*Sprime=0;
  }
}
