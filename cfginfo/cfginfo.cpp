#include "cfginfo.h"

int main (int argc, char * argv[]){

  string arg=argv[1];

//   if (arg=="-ekmax") EkMaximum(argc, argv);
if (arg=="-ekdist") EkDistribution(argc, argv);
//   else if (arg=="-ekz") EkZ(argc, argv);
//   else if (arg=="-spdist") SpeedDistribution(argc, argv);
//   else if (arg=="-vdist") VelocityDistribution(argc, argv);
//   else if (arg=="-zdist") ZDistribution(argc, argv);
//   else if (arg=="-rdist") RDistribution(argc, argv);
//   else if (arg=="-msqd") MSqD(argc, argv);
//   else if (arg=="-ndens") NDensity(argc, argv);
//   else if (arg=="-zdiff") ZDiff(argc, argv);
 else if (arg=="-tellme") InfoTellme(argc, argv);
//   else if (arg=="-amorph") FindAmorph(argc,argv);
//   else if (arg=="-nbrs") InfoNeighborDist(argc,argv);
//   else if (arg=="-rdf") RadDistFun(argc, argv);
 else if (sfind(arg, "-u")) InfoUptake(argc,argv);
 else if (sfind(arg, "-x")) InfoXPS(argc,argv);
 else if (sfind(arg, "-dep")) InfoDepth(argc,argv); 
 else if (arg=="-traj") InfoTraj(argc,argv);
//   else if (arg=="-coord") InfoCoord(argc,argv);
  else cout<<arg<<" not a recognized command.\n"; 
};

