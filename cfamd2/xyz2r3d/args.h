/* xyz2r3d.c */
#ifndef ARGS_H
#define ARGS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "et.h"
#include "point.h"
#include "r3d_utils.h"
#include "colormap.h"

#define IGNORE_BADWORDS 0
#define TRAP_BADWORDS   1
#define NO_RECURSION    2
#define MAXCHAR		255

void usage (FILE * fp);
void arg_handler (int argc, char * argv[], short bwc);
void parseline (char * ln, char * argv[], int * argc);
void echo_args (int argc, char * argv[]);
#endif
