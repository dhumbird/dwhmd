#include <stdio.h>
#include <stdlib.h>

typedef struct POINT
{
    double x, y, z;
} pt;

typedef struct ATOM
{
    int Z; /* symbol */
    pt r;  /* position */
} atom;

#define MAXNUMATOMS 5000
atom cfg[MAXNUMATOMS];
int nAtom_=0;
pt Lr_, half_Lr_;
pt cnrs[8];
static char ln[255];
void xyz_scan (FILE * fp)
{   
    char * p;
    int nc=0, i=0;
    char dumsym[5];
    if (!fp) return;
    fgets(ln, 255, fp);
    sscanf(ln, "%i", &nAtom_);
    fgets(ln, 255, fp); /* skip line */
    while (fgets(ln, 255, fp))
    {
	if (ln[0]=='H') 
	{
	    if (nc==8) {printf("Error: too many corners\n");exit(0);}
	    p=ln+1; while (isspace(*p)) p++;
	    sscanf(p, "%lf %lf %lf", 
		&(cnrs[nc].x), &(cnrs[nc].y), &(cnrs[nc].z));
	    nc++;
	}
	else
	{
	    if (i==MAXNUMATOMS) {printf("Error: too many atoms\n");exit(0);}
	    sscanf(ln, "%s %lf %lf %lf", 
		dumsym, &(cfg[i].r.x), &(cfg[i].r.y), &(cfg[i].r.z));
	    cfg[i].Z=et_s2i(dumsym);
	    i++;
	}
    }
    nAtom_=i;
}

void main (int argc, char * argv[])
{
    FILE * ifp = fopen(argv[1], "r");
    int i;
    
    xyz_scan(ifp);
    fclose(ifp);
    
    for (i=0;i<8;i++) printf("corner: %lf %lf %lf\n", 
	cnrs[i].x, cnrs[i].y, cnrs[i].z);
    for (i=0;i<nAtom_;i++) printf("atom: %i %s %lf %lf %lf\n", 
	cfg[i].Z, et_i2s(cfg[i].Z), cfg[i].r.x, cfg[i].r.y, cfg[i].r.z);
}
