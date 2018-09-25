#ifndef GLOBAL_TESTING_H
#define GLOBAL_TESTING_H

#include "pre_define.h"
using namespace std;

bool  gl_bonf(double p[], int col, double alpha);
bool gl_chi(double p[], int col, double alpha, double mu0, double sig0);
bool gl_fish(double p[], int col, double alpha);
bool gl_sims(double p[], int col, double alpha);

#endif
