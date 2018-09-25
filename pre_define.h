// Define stuct and global constants
#ifndef PRE_DEFINE_H
#define PRE_DEFINE_H
#include <vector>
#include <cmath> // for sqrt() and acos()
using namespace std;

static const int GLMAX = 4; // the maximal global testing methods I will use
static const int MULMAX = 5; // the maximal multiple testing methods I will use
static const int ROWMAX = 500; // The maximal row for p-value matrix
static const int COLMAX = 10000; // The maximal column for p-value matrix
static const int GLNEED = 2; // I need alpha and power for global tests
static const int MULNEED = 3; // I need FWER FDR power for multiple tests

struct INPUT
{
    double alpha; int row; int col; double mu0; double sig0;double mu1; double sig1;
    double ratio1; double ratio2; int rep_time; int seed;
    int index[COLMAX];
    vector<int> gl;
    vector<int> mul;
    vector <int> rejected;
    vector <double> MU1;
    vector <double> RATIO2;
};

struct P
{
    double data[COLMAX]; // temp data row by row
    bool identity[COLMAX];
    double value[COLMAX];
    double value_sort[COLMAX];
};

// All parameters used in CDF approximation
static const double RT2PI = sqrt(4.0*acos(0.0));
static const double SPLIT = 7.07106781186547;
static const double N0 = 220.206867912376;
static const double N1 = 221.213596169931;
static const double N2 = 112.079291497871;
static const double N3 = 33.912866078383;
static const double N4 = 6.37396220353165;
static const double N5 = 0.700383064443688;
static const double N6 = 3.52624965998911e-02;
static const double M0 = 440.413735824752;
static const double M1 = 793.826512519948;
static const double M2 = 637.333633378831;
static const double M3 = 296.564248779674;
static const double M4 = 86.7807322029461;
static const double M5 = 16.064177579207;
static const double M6 = 1.75566716318264;
static const double M7 = 8.83883476483184e-02;
#endif
