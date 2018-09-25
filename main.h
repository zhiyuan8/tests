// Define and initialize variables
#ifndef MAIN_H
#define MAIN_H
#include "pre_define.h"
#include <chrono>  // for high_resolution_clock

// This part should be user-defined later
INPUT input = { 0.05 /*alpha*/,
500/*row*/,   3000/*col*/,
0/*mu0*/,    1/*sig0*/,
0/*mu1*/,    1/*sig1*/,
0.5/*ratio1*/,  0.2 /*ratio2*/,
2/*rep_time*/, 10/*seed*/,
              {},{},{},{},{},{}};

P p; // define P struct

double pValue_sort[COLMAX] = { 1 };
bool global[GLMAX][ROWMAX]; // 1/0 matrix, the first index is # of methods, the second index is # of rows
double multiple[MULMAX][ROWMAX]; // matrix for rejected p-value threshold
double confusion_gl[GLMAX][GLNEED] = { 0 };// alpha, power for a total p matrix
double confusion_mul[MULMAX][MULNEED] = { 0 };// V, FDP, mean power for a total p matrix

int a;
int n_gl;
int n_mul;
int n_effect;
int n_ratio2;
int position = 0;
double generate_time = 0;
double startTime = 0;
double endTime = 0;
double testTime = 0;

auto start = std::chrono::high_resolution_clock::now();
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;

#endif
