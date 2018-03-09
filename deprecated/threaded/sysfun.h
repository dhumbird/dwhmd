//sysfun.h
//various functions involving system calls

#ifndef SYSFUN_H
#define SYSFUN_H

#include <string>
#include <unistd.h>
#include <csignal>
#include <cstdio>
#include <fstream>

string Hostname();
string DateTime();
bool FileExists(string);
int system(char *);
string EstTimeFinish(int);

#endif
