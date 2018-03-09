#include "args.h"
/* 
 * md series 2,  xyz2r3d 
 * (c) 1999 cameron abrams & regents of uc berkeley
 * dept of chemical engineering
 * university of california, berkeley 
 * 
 *
 * args:  the command-line file argument handling module.
 * Globals are assigned based on command-line keyword/value pairs.
 * Arguments are processed in
 * the order in which they are issued.
 *
 */

/* Here is the list of globals accessible by this module: */
/* render Header parameters and their default values */
extern double SecLtContr_;
extern double AmbLtContr_;
extern double SpecRefComp_;
extern double Eyepos_;

extern int xTiles_, yTiles_;
extern int xPixPerTile_, yPixPerTile_;
extern int Scheme_;
extern int Shadows_;
extern int Phong_;

extern color_t * bgClr;

extern const ptPtr MainLtPos;

extern short transparency_;
extern short RotateLightWithObject_;

extern char * rot_label[];
extern double normMax_;
extern int c_flow;
extern int color_mode;
extern double zmin, zmax;
extern pt Shift, *shift;
extern pt  per_shift;
extern pt  orbit_center,  *oc;
extern int oc_Lr[3];
extern int rot_order[MAXROTS], nRots;
extern double rots[MAXROTS];
extern double scale;
extern short drawPBCs, frntPbcCt;
extern double pbcHt, pbcClrty;
extern short noatoms;
extern color_t * pbcColor;
extern color_t * hlclrs[500];
extern int hlIDs[500], nhl, nhlc;
extern char * dof;

short quiet_;		    /* flag to supress headers and 
			     * messages in the output */
extern pt Lr_, half_Lr_;    /* Box size and 1/2 box size */
extern ipt per_;	    /* periodicity flags (def=110) */

extern char cfgName_[];	    /* Name of config -- not the file name */
extern char cfgCreateDate_[]; /* string cfg creation date */
extern char et_setupfile[]; /* datafile for element table data */
extern short elem_only_;
extern int elem_[];
extern int neo_;

short echo_=0;		    /* echos all parameter assignments */

char setup_infile_[MAXCHAR]="";
char setup_outfile_[MAXCHAR]="";


static char buf0[MAXCHAR]="";
void arg_handler (int argc, char * argv[], short bwc)
{
    int i=0,a=0,l=0, j=0;
    char * p;
    short found=0;
    void help (char *);
    void scan_setup (char *);
    void output_setup (char *);
        
    j = 0;
    nRots = 0;
    for (i = 2; i < argc; i++)
    {
	if (!strcmp(argv[i], "-xr"))
	{
	    if (nRots == MAXROTS)
	    {
		printf("Error -- too many rotations requested.\n");exit(0);
	    }
	    else
	    {
		rots[nRots] = atof(argv[++i]);
		rot_order[nRots] = X;
		nRots++;
	    }
	}
	else if (!strcmp(argv[i], "-yr"))
	{
	    if (nRots == MAXROTS)
	    {
		printf("Error -- too many rotations requested.\n");exit(0);
	    }
	    else
	    {
		rots[nRots] = atof(argv[++i]);
		rot_order[nRots] = Y;
		nRots++;
	    }
	}
	else if (!strcmp(argv[i], "-zr"))
	{
	    if (nRots == MAXROTS)
	    {
		printf("Error -- too many rotations requested.\n");exit(0);
	    }
	    else
	    {
		rots[nRots] = atof(argv[++i]);
		rot_order[nRots] = Z;
		nRots++;
	    }
	}
	else if (!strcmp(argv[i], "-hlids"))
	{
	    nhl = 0;
	    i++;
	    while (i<argc&&isdigit(argv[i][0]))
		hlIDs[nhl++] = atoi(argv[i++]);
	    i--;
	}
	else if (!strcmp(argv[i], "-hlclrs"))
	{
	    nhlc = 0;
	    i++;
	    for (; i < argc && argv[i][0] != '-'; i++)
	    {
		if (isdigit(argv[i][0]))
		    hlclrs[nhlc] = color_init(hlclrs[j++], 
				   atof(argv[i]), 
				   atof(argv[++i]), 
				   atof(argv[++i]));
		else
		    hlclrs[nhlc++] = colorPtr(argv[i]);
	    }
	    i--;
	}
	else if (!strcmp(argv[i], "-df")) dof = argv[++i];
	else if (!strcmp(argv[i], "-xtr")) shift->x = atof(argv[++i]);
	else if (!strcmp(argv[i], "-ytr")) shift->y = atof(argv[++i]);
	else if (!strcmp(argv[i], "-ztr")) shift->z = atof(argv[++i]);
	else if (!strcmp(argv[i], "-xs")) per_shift.x = atof(argv[++i]);
	else if (!strcmp(argv[i], "-ys")) per_shift.y = atof(argv[++i]);
	else if (!strcmp(argv[i], "-zs")) per_shift.z = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-zoom")) scale = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-norm")) normMax_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-zmin")) zmin = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-zmax")) zmax = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-xt")) xTiles_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-yt")) yTiles_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-xp")) xPixPerTile_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-yp")) yPixPerTile_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-scheme")) Scheme_ = atoi(argv[++i]);
 	else if (!strcmp(argv[i], "-phong")) Phong_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-slc")) SecLtContr_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-alc")) AmbLtContr_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-src")) SpecRefComp_ = atof(argv[++i]);
 	else if (!strcmp(argv[i], "-eyepos")) Eyepos_ = atof(argv[++i]);
  	else if (!strcmp(argv[i], "-noatoms")) noatoms = 1;
 	else if (!strcmp(argv[i], "-elemonly"))
	{
	    elem_only_ = 1;
	    neo_ = 0;
	    i++;
	    while (i < argc && argv[i][0] != '-')
		elem_[neo_++] = et_s2i(argv[i++]);
	    i--;
	}
 	else if (!strcmp(argv[i], "-pbc")) drawPBCs = 1;
 	else if (!strcmp(argv[i], "-pbcHt")) 
	{
	    drawPBCs = 1;
	    pbcHt = atof(argv[++i]);
 	}
	else if (!strcmp(argv[i], "-pbcClr")) 
	{
	    pbcClrty = atof(argv[++i]);
 	    drawPBCs = 1;
	}
	else if (!strcmp(argv[i], "-pbcCtwy")) frntPbcCt = 1;
 	else if (!strcmp(argv[i], "-pbcColor"))
	{
	    if (isdigit(argv[++i][0]))
	    {
		pbcColor = color_init(pbcColor, atof(argv[i]), atof(argv[++i]), atof(argv[++i]));
	    }
	    else pbcColor = colorPtr(argv[i]);
	}
	else if (!strcmp(argv[i], "-bg"))
	{
	    if (isdigit(argv[++i][0]))
	    {
		bgClr = color_init(bgClr, atof(argv[i]), atof(argv[++i]), atof(argv[++i]));
	    }
	    else bgClr = colorPtr(argv[i]);
	}
 	else if (!strcmp(argv[i], "-mlpos"))
	{
	    MainLtPos->x = atof(argv[++i]);
	    MainLtPos->y = atof(argv[++i]);
	    MainLtPos->z = atof(argv[++i]);
	}
 	else if (!strcmp(argv[i], "-noshadows")) 
	{
	    Shadows_ = 0;
	}
 	else if (!strcmp(argv[i], "-orbit"))
	{
	    RotateLightWithObject_ = 1;
	}
 	else if (!strcmp(argv[i], "-oc"))
	{
	    RotateLightWithObject_ = 1;
	    i++;
	    if (argv[i][0] == 'X')
	    {
		oc_Lr[0] = 1;
		sscanf(&(argv[i][1]), "%lf", &(oc->x));
	    }
	    else if (argv[i][0] == 'x')
	    {
		oc_Lr[0] = -1;
		sscanf(&(argv[i][1]), "%lf", &(oc->x));
	    }
	    else orbit_center.x = atof(argv[i]);
	    i++;
	    if (argv[i][0] == 'Y')
	    {
		oc_Lr[1] = 1;
		sscanf(&(argv[i][1]), "%lf", &(oc->y));
	    }
	    else if (argv[i][0] == 'y')
	    {
		oc_Lr[1] = -1;
		sscanf(&(argv[i][1]), "%lf", &(oc->y));
	    }
	    else orbit_center.y = atof(argv[i]);
	    i++;
	    if (argv[i][0] == 'Z')
	    {
		oc_Lr[2] = 1;
		sscanf(&(argv[i][1]), "%lf", &(oc->z));
	    }
	    else if (argv[i][0] == 'z')
	    {
		oc_Lr[2] = -1;
		sscanf(&(argv[i][1]), "%lf", &(oc->z));
	    }
	    else orbit_center.z = atof(argv[i]);
	}
 	else if (!strcmp(argv[i], "-colmode"))
	{
	    i++;
	    if (!strcmp(argv[i], "dep")) color_mode = 1;
	    else color_mode = 0;
	}
	else
	{
	    printf("Error:  Keyword %s not recognized.\n", argv[i]);
	    usage(stdout);
	    exit(-1);
	}
    }
}


void usage (FILE * fp)
{
    if (!fp) return;
    fprintf(fp, "Usage: xyz2r3d <xyzName> [control params] [> <r3dName>]\n");
    fprintf(fp, "\t<xyzName>   = name of cfg data file\n");
    fprintf(fp, "\t<r3dName>   = name of r3d data file\n");
    fprintf(fp, "\t-df <directions output file>\n");
    fprintf(fp, "Rotation Controls: Values are in degrees\n");
    fprintf(fp, "\t -xr <x-rotation> -yr <y-rotation> -zr <z-rotation>\n");
    fprintf(fp, "Translation Controls: Values are in cfg-unit length\n");
    fprintf(fp, "\t -xtr <x-translation> -ytr <y-translation> -ztr <z-translation>\n");
    fprintf(fp, "Shifting Controls: Values are in cfg-unit length\n");
    fprintf(fp, "\t -xs <x-shift> -ys <y-shift> -zs <z-shift>\n");
    fprintf(fp, "(Shifting means that the actual atom coordinates are\n");
    fprintf(fp, "translated *prior* to rendering the image.)\n");
    fprintf(fp, "Zoom and Normalizing Length Controls:\n");
    fprintf(fp, "\t -zoom <scaleFactor> -norm <normMax>\n");
    fprintf(fp, "Colormode:  Natural, KE-resolved, or Depth-resolved\n");
    fprintf(fp, "\t -colmode <nat|ke|dep|pe>\n");
    fprintf(fp, "\t -kemin <kemin,eV> -kemax <kemax,eV> -cooltrans <YES|NO> for ke colormode only\n");
    fprintf(fp, "\t -pemin <pemin,eV> -pemax <pemax,eV> for pe colormode only\n");
    fprintf(fp, "\t -zmin <zmin> -zmax <zmax> for dep colormode only\n");
    fprintf(fp, "Object controls:\n");
    fprintf(fp, "\t -noatoms does not render atoms\n");
    fprintf(fp, "\t -ionid <ionIndex> colors ion\n");
    fprintf(fp, "\t -ioncolor <r g b>\n");
    fprintf(fp, "\t -pbc renders periodic boundary planes\n");
    fprintf(fp, "\t -pbcHt <float> height extent of vertical pbc planes\n");
    fprintf(fp, "\t -pbcClr <float> clarity of pbc planes\n");
    fprintf(fp, "\t -pbcColor <rgb or color key> color of pbc planes\n");
    fprintf(fp, "\t -pbcCtwy cuts the -y pbc plane off at x for z > 0.0\n");
    fprintf(fp, "\t -velarrfor <id-list> atoms for which velocity arrows are drawn\n");
    fprintf(fp, "\t -velarrlf <float> length of velocity arrows\n");
    fprintf(fp, "\t -velarrshaftrad <float> radius of velocity arrow shafts\n");
    fprintf(fp, "\t -velarrheadw <float> width of velocity arrow heads\n");
    fprintf(fp, "\t -velarrheadl <float> length of velocity arrow heads\n");
    fprintf(fp, "\t -velarrheadl <float> length of velocity arrow heads\n");
    fprintf(fp, "\t -velarrheadnface <int> number of faces of pyramidal arrow heads\n");
    fprintf(fp, "\t -hlids <id-list> list of id's for atoms that are to\n");
    fprintf(fp, "\t                 be highlighted.\n");
    fprintf(fp, "\t -hlclrs <color-namelist> list of colors for highlighted atoms\n");
    fprintf(fp, "\t -onlyelem <elem-str> only atoms of this element are rendered.\n");
    fprintf(fp, "Raster3D render Header Macros:\n");
    fprintf(fp, "\t -xt <xTiles(32)> -yt <yTiles(32)>\n");
    fprintf(fp, "\t -xp <xPixelsPerTile(18)> -yp <yPixelsPerTile(18)>\n");
    fprintf(fp, "\t -scheme <scheme(3)>\n");
    fprintf(fp, "\t -bg (colorName or <r> <g> <b> (0 0 0))\n");
    fprintf(fp, "\t -noshadows -phong <phongPower(25)>\n");
    fprintf(fp, "\t -slc <secondaryLightContributionFraction(0.15)>\n");
    fprintf(fp, "\t -alc <ambientLightContributionFraction(0.05>\n");
    fprintf(fp, "\t -src <specularReflectionComponentFraction(0.25)>\n");
    fprintf(fp, "\t -eyepos <EYEPOS(4.0)>\n");
    fprintf(fp, "\t -mlpos <x> <y> <z> (1 1 1) (Main Light Position)\n");
    fprintf(fp, "\t -orbit or -oc <x> <y> <z> Orbit Center\n");
    fprintf(fp, "Orbit means rotate the light source direction with any\n");
    fprintf(fp, "object rotation to simulate the camera orbiting a stationary\n");
    fprintf(fp, "object.  The (x, y, z) coordinates are the orbit center.\n");
    fprintf(fp, "The -orbit keyword assumes the orbit center is 0, 0, 0.\n");
    fprintf(fp, "The -oc keyword requires exactly 3 numbers to specify\n");
    fprintf(fp, "the orbit center if the user wishes it to be something\n");
    fprintf(fp, "other than (0, 0, 0).  Instead of a floating point value, \n");
    fprintf(fp, "you may specify 'X' for +Box.x/2 and 'x' for -Box.x/2, \n");
    fprintf(fp, "and likewise for <y> and <z> orbit center coords.\n");
}
