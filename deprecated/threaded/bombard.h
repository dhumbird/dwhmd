//bombard.h
//Perform ion bombardment runs

#include "errors.h"
#include "config.h"
#include "consts.h"
#include "sysfun.h"
#include "strfun.h"
#include <fstream>

void AddMultiIons(int, config * cfg);
void DelMultiIons(int, config * cfg);
int IonBombard(int argc, char * argv[]);
