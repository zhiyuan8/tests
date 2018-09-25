#include <iostream>
#include <random>//std::default_random_engine generator
#include <cstdlib>  // std::rand, std::srand
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <chrono>  // for high_resolution_clock

#include "pre_define.h"
#include "utils.h"
using namespace std;

void generate(int current_line, P &p, INPUT &input, int a)
{
    seed_seq seq{ current_line };
    mt19937 generator(seq); // I cannot explain why this works in setting seed
                            //generator.seed(input.seed*i*current_line);
                            //std::default_random_engine generator;
    normal_distribution<double> distribution0(input.mu0, input.sig0); // H0 distribution
    normal_distribution<double> distribution1(input.mu1, input.sig1); // H1 distribution
    if (input.rejected[current_line] == 0) // this line has no rejection
    {
        for (int i = 0; i < input.col; i++)
        {

            p.data[i] = distribution0(generator);
            p.identity[i] = false;
            p.value[i] = cdf_approx(a*(p.data[i] - input.mu0) / input.sig0);//use CDF approximation
            p.value_sort[i] = p.value[i];
        }
    }
    else  // this line has rejection
    {
        random_shuffle(input.index, input.index+input.col);// shuffle index
        for (int i = 0; i < input.rejected[current_line]; i++)
        {
            int x = input.index[i];
            p.data[x] = distribution1(generator);
            p.identity[x] = true;
            p.value[x] = cdf_approx(a*(p.data[x] - input.mu0) / input.sig0);//use CDF approximation
            p.value_sort[x] = p.value[x];
        }
        for (int i = input.rejected[current_line]; i < input.col; i++)
        {
            int x = input.index[i];
            p.data[x] = distribution0(generator);
            p.identity[x] = false;
            p.value[x] = cdf_approx(a*(p.data[x] - input.mu0) / input.sig0);//use CDF approximation
            p.value_sort[x] = p.value[x];
        }
    }
}
