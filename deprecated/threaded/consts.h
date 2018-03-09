//consts.h

#ifndef CONSTS_H
#define CONSTS_H

#define KB 8.6173857e-5
#define PI 3.1415926559
#define ONE_HALF 0.5
#define LAMBDA_SW 21.0
#define A_SW 7.049556277
#define B_SW 0.6022245584
#define EPS_SW 2.16966
#define SIGMA_SW 2.0951
#define a_MO_AR_SI 0.09734

const double ONE_THIRD=1.0/3.0;
const double ONE_SIXTH=1.0/6.0;
const double TWO_THIRDS=2.0/3.0;
const double FIVE_SIXTHS=5.0/6.0;

const double B_SIGMA4_SW=B_SW*pow(SIGMA_SW, 4);
const double a_SW=1.8*SIGMA_SW;
const double a_sq_SW=a_SW*a_SW;
const double GAMMA_SW=1.2*SIGMA_SW;
const double GEPS_SW=GAMMA_SW*EPS_SW;
const double EPS2_SW=EPS_SW*2;
const double EPSA_SW=EPS_SW*A_SW;

const double EPS_MO_AR_SI=14.39965*14*18;

#endif



