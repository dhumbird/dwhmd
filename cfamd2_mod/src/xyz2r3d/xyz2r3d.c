/*
 * MD series 2 -- Molecular Dynamics of Plasma-Surface Chemistry
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 * Module Name:  xyz2r3d:  creates a Raster3D-format datafile from
 * the supplied xyz (Winter-style) config. 
 */

#include <stdio.h>
#include "et.h"
#include "dblmat.h"
#include "colormap.h"
#include "point.h"
#include "r3d_utils.h"

#include "cam_colors.ext"

double version_=1.00;
int build_=1;

typedef struct ATOM * atomPtr;
typedef struct ATOM
{
    int id; /* unique index */
    int Z; /* symbol */
    pt r;  /* position */
} atom;
void r3d_ball (FILE * fp, atomPtr L, double Max);
#define MAXNUMATOMS 5000
atom cfg[MAXNUMATOMS];
int nAtom_=0;
ipt per_={1, 1, 0};
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
	    cfg[i].id=i++;
	}
    }
    nAtom_=i;
    /* Determine the boxsize from the corners */
    ptPtr_scalmult(&Lr_, &(cnrs[0]), 2.0);
    ptPtr_scalmult(&half_Lr_, &Lr_, 0.5);
}


enum {NATURAL, DEPTH_RES};
char * color_modeLabels[4] =
{
    "Natural", "Depth-resolved"
};
extern float hue_max_;	/* colormap.c */
extern float hue_min_;	/* colormap.c */
color_t c0, *cp0 = &c0;
color_t c1, *cp1 = &c1;



extern char * rot_label[];
double normMax_ = 0.0;
int c_flow = -1;
int color_mode = NATURAL;
double zmin = 0.0, zmax = 0.0;
pt Shift ={0.0, 0.0, 0.0}, *shift=&Shift;
pt  per_shift = {0, 0, 0};
pt  orbit_center = {0, 0, 0},  * oc = &orbit_center;
int oc_Lr[3] = {0, 0, 0};
int rot_order[MAXROTS], nRots = 0;
double rots[MAXROTS];
double scale = 1.0;
short drawPBCs = 0, frntPbcCt = 0;
double pbcHt = 0.0, pbcClrty = 0.6;
short noatoms = 0;
color_t * pbcColor = &peagreen;
color_t * hlclrs[500];
int	hlIDs[500], nhl=0, nhlc=0;
char * dof = NULL;

/* render Header parameters and their default values */
extern double ImageBufferSize_;
extern int xTiles_, yTiles_;
extern int xPixPerTile_, yPixPerTile_;
extern int Scheme_;
extern color_t * bgClr;
extern int Shadows_;
extern int Phong_;
extern double SecLtContr_;
extern double AmbLtContr_;
extern double SpecRefComp_;
extern double Eyepos_;
extern const ptPtr MainLtPos;
extern short transparency_;
extern short RotateLightWithObject_;
short elem_only_ = 0;
int elem_[] = {0, 0, 0, 0, 0, 0};
int neo_ = 0;

extern element et_[];
extern short quiet_;

void main (int argc, char * argv[])
{
    atomPtr a = NULL, hL = NULL, b = NULL;
    FILE * fp = NULL;
    pt TMP_PT={0, 0, 0}, *tmppt=&TMP_PT;
    char * cfn;
    int i = 0,  j = 0;
    double span = 0.0, x, y, z;
    void r3d_atomballs (FILE * fp, atomPtr L, double span);
    void r3d_clratomballs (FILE * fp, atomPtr L, color_t **clrs, double span);
    ptPtr xyz = NULL, xyZ = NULL, xYz = NULL, xYZ = NULL, 
	  Xyz = NULL, XyZ = NULL, XYz = NULL, XYZ = NULL;
    
    color_init(cp0, hue_min_, hue_min_, hue_max_);
    color_init(cp1, hue_max_, hue_min_, hue_min_);

    quiet_=1;
    et_initialize();
    
    for (i = 0; i < MAXROTS; i++) rot_order[i] = 0;
    for (i = 0; i < 500; i++)
    {
	hlIDs[i] = 0;
	hlclrs[i] = NULL;
    }
    if (argc < 2)
    {
	usage(stdout);
	exit(-1);
    }

    cfn=argv[1];
    arg_handler(argc, argv, 0);
    
    i = 0;
    while (i < 500 && hlclrs[i]) i++;
    if (i == 0) hlclrs[i++] = &white;
    for (j = i; j < 500; j++) hlclrs[j] = hlclrs[i-1];
	
    fp = fopen(cfn, "r");
    if (fp)
    {
	xyz_scan(fp);
	fclose(fp);
	if (oc_Lr[0]) orbit_center.x += oc_Lr[0]*half_Lr_.x;
	if (oc_Lr[1]) orbit_center.y += oc_Lr[1]*half_Lr_.y;
	if (oc_Lr[2]) orbit_center.z += oc_Lr[2]*half_Lr_.z;
	if (drawPBCs && pbcHt == 0.0) pbcHt = Lr_.z;
	for (i=0;i<nAtom_;i++)
	    ptPtr_minimg(&(cfg[i].r), ptPtr_add(&(cfg[i].r), &(cfg[i].r), &per_shift));
	for (i=0;i<nAtom_;i++) 
	    ptPtr_add(&(cfg[i].r), ptPtr_scalmult(&TMP_PT, oc, -1.0), &(cfg[i].r));
		
	if (!normMax_)
	{
	    normMax_ = sqrt(ptPtr_sqabs(&Lr_));
	}
	span = normMax_;
	/* add the default buffer to the image span */
	span += ImageBufferSize_;
	/* find the max and min z-values */
	if ((!zmin || !zmax) && color_mode == DEPTH_RES) 
	    for (i=0;i<nAtom_;i++)
	    {
		if (cfg[i].r.z < zmin) zmin = cfg[i].r.z;
		if (cfg[i].r.z > zmax) zmax = cfg[i].r.z;
	    }
	
	if (dof)
	{
	    fp = fopen(dof, "w");
	    fprintf(fp, "# %s direction file %s\n", argv[0], dof);
	    fprintf(fp, "cfg %s\n", argv[1]);
	    fprintf(fp, "norm %.10lf\n", normMax_);
	    for (i = 0; i < MAXROTS && rot_order[i] != 0; i++)
		fprintf(fp, "%s %.10lf\n", rot_label[rot_order[i]], rots[i]);
	    fprintf(fp, "zoom %.10lf\n", scale);
	    fprintf(fp, "oc.x %.10lf\n", oc->x);
	    fprintf(fp, "oc.y %.10lf\n", oc->y);
	    fprintf(fp, "oc.z %.10lf\n", oc->z);
	    fprintf(fp, "tr.x %.10lf\n", shift->x);
	    fprintf(fp, "tr.y %.10lf\n", shift->y);
	    fprintf(fp, "tr.z %.10lf\n", shift->z);
	    if (drawPBCs) fprintf(fp, "pbc\n");
	    fprintf(fp, "mlp.x %.10lf\n", MainLtPos->x);
	    fprintf(fp, "mlp.y %.10lf\n", MainLtPos->y);
	    fprintf(fp, "mlp.z %.10lf\n", MainLtPos->z);
	    fprintf(fp, "#end of direction file %s\n", dof);
	    fclose(fp);
	}
	
	
	r3d_header(stdout, rots, rot_order, scale, shift, span);
	printf("# xyz2r3d version %.2lf by cam abrams\n", version_);
	printf("# <cfgName> = %s\n", cfn);
        printf("# Box size is %.5le %.5le %.5le\n", Lr_.x, Lr_.y, Lr_.z);
	if (rot_order[0] != 0)
	{
	    printf("# Rotation order is: ");
	    for (i = 0; i < MAXROTS && rot_order[i] != 0; i++)
		printf("%s[%.2lf] ", rot_label[rot_order[i]], rots[i]);
	    printf("\n");
	}
	if (scale) printf("# Scale is %.3lf.\n", scale);
	printf("# nAtom = %i\n", nAtom_);
	printf("# Color mode is %s\n", color_modeLabels[color_mode]);
	switch(color_mode)
	{
	    case DEPTH_RES: printf("# z-min is %.5lf\n", zmin);
			    printf("# z-max is %.5lf\n", zmax);
			    break;
	    case NATURAL:   printf("# natural coloring\n");
			    break;			    
	}
	printf("# span is %.5lf\n", span);
	if (!noatoms) 
	{
	    for (i=0;i<nAtom_;i++)
	    {
		r3d_ball (stdout, &(cfg[i]), span);
	    }
	}
    }
    
    if (drawPBCs)
    {
	xyz = ptPtr_create((half_Lr_.x-oc->x)/span, 
			 (half_Lr_.y-oc->y)/span, 
			 (-half_Lr_.z+pbcHt-oc->z)/span);
	Xyz = ptPtr_create((-half_Lr_.x-oc->x)/span, 
			 (half_Lr_.y-oc->y)/span, 
			 (-half_Lr_.z+pbcHt-oc->z)/span);
	XYz = ptPtr_create((-half_Lr_.x-oc->x)/span, 
			(-half_Lr_.y-oc->y)/span, 
			((-half_Lr_.z+pbcHt-oc->z))/span);
	XyZ = ptPtr_create((-half_Lr_.x-oc->x)/span, 
			(half_Lr_.y-oc->y)/span, 
			(-half_Lr_.z-oc->z)/span);
	XYZ = ptPtr_create((-half_Lr_.x-oc->x)/span,
			(-half_Lr_.y-oc->y)/span, 
			(-half_Lr_.z-oc->z)/span);
	xYz = ptPtr_create((half_Lr_.x-oc->x)/span, 
			(-half_Lr_.y-oc->y)/span, 
			((-half_Lr_.z+pbcHt-oc->z))/span);
	xYZ = ptPtr_create((half_Lr_.x-oc->x)/span, 
			(-half_Lr_.y-oc->y)/span, 
			(-half_Lr_.z-oc->z)/span);
	xyZ = ptPtr_create((half_Lr_.x-oc->x)/span, 
			(half_Lr_.y-oc->y)/span, 
			(-half_Lr_.z-oc->z)/span);

	printf("# pbcs:\n");
	r3d_transparency_on(stdout, pbcClrty);
	r3d_tri(stdout, XYz, XYZ, XyZ, *pbcColor);
	r3d_tri(stdout, XYz, Xyz, XyZ, *pbcColor);
	r3d_tri(stdout, Xyz, xyz, xyZ, *pbcColor);
	r3d_tri(stdout, Xyz, XyZ, xyZ, *pbcColor);
	r3d_tri(stdout, xyz, xyZ, xYZ, *pbcColor);
	r3d_tri(stdout, xyz, xYz, xYZ, *pbcColor);
	if (frntPbcCt)
	{
	    xyz->z = 0.0;
	    Xyz->z = 0.0;
	    XYz->z = 0.0;
	    xYz->z = 0.0;
	}
	r3d_tri(stdout, xYz, xYZ, XYZ, *pbcColor);
	r3d_tri(stdout, xYz, XYz, XYZ, *pbcColor);
	
	r3d_transparency_off(stdout);
	r3d_tri(stdout, XYZ, XyZ, xyZ, *pbcColor);
	r3d_tri(stdout, XYZ, xYZ, xyZ, *pbcColor);
    }
    
    printf("# Program Ends.\n");
}

color_t c, *cp = &c;
void r3d_ball (FILE * fp, atomPtr L, double Max)
{
    if (!fp || !L) return;
    color_init(cp, et_[L->Z].color.r, 
		        et_[L->Z].color.g, 
			et_[L->Z].color.b );
    if (color_mode == DEPTH_RES)
	color_mapFloat(zmin, cp0, zmax, cp1, L->r.z, cp, c_flow);
    if (!elem_only_ || (elem_only_ && 
		       (L->Z == elem_[0] || L->Z == elem_[1] || 
			L->Z == elem_[2] || L->Z == elem_[3] || 
			L->Z == elem_[4])))
    {
	fprintf(fp, "# atom %s_{%i}, r = %.4lf, rgb=(%.3lf %.3lf %.3lf)\n", 
	    et_[L->Z].sym, L->id, et_[L->Z].r, cp->r, cp->g, cp->b);
	fprintf(fp, "2\n");
	fprintf(fp, "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n", 
	    L->r.x/Max, L->r.y/Max, L->r.z/Max, 
	    et_[L->Z].r/Max, 
	    cp->r, cp->g, cp->b);
    }
}

