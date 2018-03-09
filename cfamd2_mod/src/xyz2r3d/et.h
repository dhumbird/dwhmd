/*
 * md  series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:  et 
 * 
 * Module Description:  et defines and initializes the periodic
 * table data structure used by other modules in the series.
 * 
 */

#ifndef ET_H
#define ET_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "colormap.h"

#define MAXNUMELEMENTS      100

typedef struct ELEMENT
{   /* record of data on an element */
    char    name[15];       /* name of element */
    char    sym[3];         /* symbol */
    double  r;	    /* radius (for render) */
    struct
    { 
	double r, g, b; 
    } color;		    /* color (for render) */
} element;	

void et_initialize (void);
int et_s2i (char * sym);
char * et_i2s (int i);
#endif

