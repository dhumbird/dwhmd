//strfun.h
//general string functions for filenames, etc.

#ifndef STRFUN_H
#define STRFUN_H

#include <string>
#include <cmath>
#include <utility>

bool sfind(string, string);
string itoa(int, int);
string ftoa(float, int);
string time2string (double);
pair<int, int> CheckIntRange (string);
pair<float, float> CheckFloatRange (string);

#endif
