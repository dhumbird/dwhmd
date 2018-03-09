#include "rcfg.h"
#include "rgb.h"
#include <fstream>
#include <utility>
#include <algorithm>
#include "d_sphere.h"
#include "d_plane.h"
#include "d_cylinder.h"
#include "d_triangle.h"

void ReadColorFile(rcfg&, string);
void KEColor(rcfg&, float, float);
void ZColor(rcfg&);
void FixSpheres(rcfg&, vector<sphere>&, short);
void FixWalls(rcfg&, vector<plane>&, float);
void FixBoundary(rcfg&, vector<cylinder>&, float);
void ZMarker(rcfg&, vector<triangle>&, float);
void BondConnect(rcfg&, vector<cylinder>&, short, float);
void AssWipe(rcfg&, vector<plane>&, float, float, float);
