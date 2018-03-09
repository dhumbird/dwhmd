//strfun.h
//general string functions for filenames, etc.

#ifndef STRFUN_H
#define STRFUN_H

#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;
bool sfind(string, string);
string time2string (double);
bool FileExists(string);
void CmdError(const char *);
#endif
