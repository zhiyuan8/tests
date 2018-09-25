#ifndef MULTIPLE_TESTING_H
#define MULTIPLE_TESTING_H

#include "pre_define.h"
using namespace std;

double  mul_bonf(double p[], int col, double alpha);
double mul_holm(double p[], int col, double alpha);
double mul_hoch(double p[], int col, double alpha);
double mul_homm(double p[], int col, double alpha);
double mul_BH(double p[], int col, double alpha);

#endif
