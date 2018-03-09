#ifndef R3D_UTIL_H
#define R3D_UTIL_H

#include "../config.h"
#include "../svector.h"
#include "../sysfun.h"
#include "rgb.h"
#include <map>
#include <set>

void TopSurf(config&, float, float, float, float);
void Plane(svector&, svector&, svector&, svector&, rgb&);
void PrintV(svector&);

#endif
