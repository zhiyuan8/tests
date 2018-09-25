#include <iostream>
using namespace std;

// remember the threshold is the biggest p value we make rejection

double mul_bonf(double p[], int col, double alpha)
{
    return(alpha / col);
}

double mul_holm(double p[], int col, double alpha)
{
    int reject = 0;
    for (int j = 0; j<col;j++)
        if (p[j] <= alpha / (col - j)) // remember index starts from 0
            reject = j + 1;
        else break;
        return  (alpha / (col - reject + 1));
}

double mul_hoch(double p[], int col, double alpha)
{
    int reject = 0;
    for (int j = col; j > 0;j--)
    {
        if (p[j - 1] > alpha / (col - j + 1)) // remember maximal index is col-1
            reject = j - 1;
        else break;
    }
    return  (alpha / (col - reject + 1));
}

double mul_homm(double p[], int col, double alpha) {
    int reject = 0; int threshold, decision;
    for (int j = col; j > 0; j--)
    {
        decision = 0; // temp # of times where p[n - j + k] > ka / j
        for (int k = 1; k <= j;k++)
        {
            if (p[col - j + k - 1] > (k*alpha / j))
                decision = decision + 1;
        }
        if (decision == j)
        {
            threshold = j;
            break;
        } //find the largest j, save it in threshold, and break out the loop
    }
    if (threshold > 0) // j exists, reject p-value <= mul_Homm$criteria
        return(alpha / threshold);
    else
        return(1);//if no such j exists, then reject all hypothesis
}

double mul_BH(double p[], int col, double alpha)
{
    int threshold = 0;
    for (int j = col; j > 0; j--)
    {
        if (p[j-1] <= j*alpha / col) // remeber maximal index is input.col-1
        {
            threshold = j;
            break;
        } //find the largest j, save it in threshold, and break out the loop
    }
    if (threshold > 0) // j exists, reject p-value <= mul_Homm$criteria
        return(alpha * threshold / col);
    else
        return(0);//if no such j exists, then all p value are too large to be rejected
}
