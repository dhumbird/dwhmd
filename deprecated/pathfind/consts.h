//consts.h

#ifndef CONSTS_H
#define CONSTS_H
//***************general******************
#define KB 8.6173857e-5
#define PI 3.1415926559
#define ONE_HALF 0.5
const double ONE_THIRD=1.0/3.0;
const double ONE_SIXTH=1.0/6.0;
const double TWO_THIRDS=2.0/3.0;
const double FIVE_SIXTHS=5.0/6.0;
const double COS103=cos(103.0/180.0*PI);
//******* SW Si/Si & Si/Si/Si ***********
#define SW_SI_SI_SI_LAMBDA 21.0
#define SW_SI_SI_A 7.049556277
#define SW_SI_SI_B 0.6022245584
#define SW_SI_SI_EPS 2.16966
#define SW_SI_SIGMA 2.0951
const double SW_SI_SIGMA2 = pow(SW_SI_SIGMA,2);
const double SW_SI_SIGMA3 = pow(SW_SI_SIGMA,3);
const double SW_SI_B_SIGMA4=SW_SI_SI_B*pow(SW_SI_SIGMA, 4);
const double SW_SI_a=1.8*SW_SI_SIGMA;
const double SW_SI_SI_SI_GAMMA = 1.2 * SW_SI_SIGMA;
const double SW_SI_SI_SI_GEPS = SW_SI_SI_SI_GAMMA * SW_SI_SI_EPS;
const double SW_SI_SI_EPS2 = SW_SI_SI_EPS*2;
const double SW_SI_SI_EPSA = SW_SI_SI_EPS * SW_SI_SI_A;


#endif



