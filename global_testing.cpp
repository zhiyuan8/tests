#include <iostream>
#include <boost/math/distributions/chi_squared.hpp>  //quantile of chi-square
#include <cmath> // use log

#include "pre_define.h"
using namespace std;

bool gl_bonf(double p[], int col, double alpha)
{
    if (p[0] <= alpha / col)
    {
        return(true);
    }
    else
    {
        return(false);
    }
    return(0);
}

bool gl_chi(double p[], int col, double alpha, double mu0, double sig0)
{
    bool result = false;
    double temp = 0;
    double z = 0;
    // calculate qchisq((1 - a), df = length(z))
    boost::math::chi_squared mydist(col);
    // first calculat 	z = (y - input$mu0) / input$sig0;
    // Then calculate sum(z ^ 2) in same loop
    for (int i = 0; i < col; i++)
    {
        z = (p[i] - mu0) / sig0;
        temp = temp + z*z;
    }
    if (temp >= boost::math::quantile(mydist, 1 - alpha))
    {
        result = true;
    }
    return(result);
}

bool gl_fish(double p[], int col, double alpha)
{
    double temp = 0;
    bool result = false;
    // calculate  qchisq((1 - a), df = 2 * length(p)))
    boost::math::chi_squared mydist(2 * col);
    // calculate  sum(-2 * log(p)
    for (int i = 0; i < col; i++)
    {
        temp = temp - 2 * log(p[i]);
    }
    if (temp > boost::math::quantile(mydist, 1 - alpha))
    {
        result = true;
    }
    return(result);
}

bool gl_sims(double p[], int col, double alpha)
{
    bool result = false;
    for (int i = 0; i < col; i++)
    {
        if (p[i] <= (i + 1)*alpha / col)
        {
            result = true;
            break;
        }
    }
    return(result);
}
