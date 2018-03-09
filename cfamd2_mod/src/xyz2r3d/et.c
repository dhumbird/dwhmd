/*
 * md series 2
 * (c) 1999 Cam Abrams, Dept. of Chemical Engineering
 * Regents of the University of California, Berkeley
 * 
 */

#include "et.h"

element et_[]=
{
    {"Water", "O", 0.5, {0.0, 0.0, 1.0}},  
    {"Head", "N", 0.5, {0.5, 0.5, 0.5}},  
    {"Tail", "C", 0.5, {1.0, 0.0, 0.0}}, 
    {"Null", "XX", 0.0, {0, 0, 0}}, 
    {"Null", "XX", 0.0, {0, 0, 0}}, 
    {"Null", "XX", 0.0, {0, 0, 0}}, 
    {"Null", "XX", 0.0, {0, 0, 0}}, 
    {"Null", "XX", 0.0, {0, 0, 0}}, 
    {"Null", "XX", 0.0, {0, 0, 0}}, 
    {"Null", "XX", 0.0, {0, 0, 0}} 
};

int et_s2i (char * sym)
{
    int i=0;
    for (i=0;i<MAXNUMELEMENTS&&et_[i].r&&strcmp(sym, et_[i].sym);i++);
    if (!strcmp(sym, et_[i].sym)) return i;
    else return -1;
}
char * et_i2s (int i)
{
    return et_[i].sym;
}

char et_setupfile[255]="et.dat";
void et_initialize (void)
{
    int i;
    double r;
    FILE * fp = NULL;
    char ln[255], dumsym[4], dumcol[50], dumname[50];
    color_t dumcolor={0, 0, 0};
    
    fp=fopen(et_setupfile, "r");
    if (fp)
    {
	while (fgets(ln, 255, fp))
	{
	    if (ln[0]!='#')
	    {
		sscanf(ln, "%i %s %lf %s %s", 
		    &i, dumsym, &r, dumcol, dumname);
		et_[i].r=r;
		strcpy(et_[i].sym, dumsym);
		strcpy(et_[i].name, dumname);
		et_[i].color.r=colorPtr(dumcol)->r;
		et_[i].color.g=colorPtr(dumcol)->g;
		et_[i].color.b=colorPtr(dumcol)->b;
	    }
	}
	fclose(fp);
    }
}

