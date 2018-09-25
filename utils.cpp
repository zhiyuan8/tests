#include <iostream>
#include <math.h> // use log
#include <algorithm>    // std::random_shuffle
#include <random>//std::default_random_engine generator
#include <vector>

#include "pre_define.h"
#include "global_testing.h"
#include "multiple_testing.h"
using namespace std;

void get_size(INPUT &input, int &n_gl, int &n_mul, int &n_effect, int &n_ratio2)
{
    n_gl = input.gl.size();
    n_mul = input.mul.size();
    n_effect = input.MU1.size();
    n_ratio2 = input.RATIO2.size();

}

void Initialize_size(INPUT &input, int n_gl, int n_mul)
{
    input.rejected.resize(input.row);
    for (int i = 0; i < input.col; i++)
    {
        input.index[i] = i;
    }
}

void Initialize_confusion(double confusion_gl[][GLNEED], double confusion_mul[][MULNEED], int n_gl, int n_mul)
{
    for (int x = 0; x < n_gl; x++)
    {
        for (int y = 0; y < GLNEED; y++)
        {
            confusion_gl[x][y] = 0;
        }
    }
    for (int x = 0; x < n_mul; x++)
    {
        for (int y = 0; y < MULNEED; y++)
        {
            confusion_mul[x][y] = 0;
        }
    }
}

void define_a(double mu0, double mu1, int &a)
{
    if (mu0 <= mu1)
        a = -1; // use area in right tail as p - value
    else
        a = 1;  // use area in left tail as p - value
}

int round_double(double number)
{
    if (number > 0.0)
    {
        return(int(number + 0.5));
    }
    else
    {
        return(int(number - 0.5));
    }
}

void generate_ground_truth(INPUT &input)
{
    int many = round_double(input.ratio2*input.col);    //How many p-value in a row are going to be rejected
    int index = round_double(input.ratio1*input.row);  //The itotal # of row to be rejected
    // If seed<0, then use a random generated seed
    if (input.seed <= 0)
    {
        std::random_device rd;
        input.seed = rd();
    }
    for (int i = 0; i < input.row; i++)
    {
        if (i<index) // from 0 to (index-1)
        {
            input.rejected[i] = many;
        }
        else
        {
            input.rejected[i] = 0;
        }
    }
    std::srand(input.seed); // set.seed here
    std::random_shuffle( input.rejected.begin(), input.rejected.end() ); // finally, I get the shuffled rejected[]
}

double cdf_approx(double x)
{
    double z = fabs(x);
    double c = 0.0;
    if (z <= 37.0)
    {
        const double e = exp(-z*z / 2.0);
        if (z<SPLIT)
        {
            const double n = (((((N6*z + N5)*z + N4)*z + N3)*z + N2)*z + N1)*z + N0;
            const double d = ((((((M7*z + M6)*z + M5)*z + M4)*z + M3)*z + M2)*z + M1)*z + M0;
            c = e*n / d;
        }
        else
        {
            const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
            c = e / (RT2PI*f);
        }
    }
    return x <= 0.0 ? c : 1 - c;
}

void quickSort(double arr[], int left, int right)
{
    int i = left, j = right;
    double tmp;
    double pivot = arr[(left + right) / 2];
    // partition
    while (i <= j) {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    };
    // recursion
    if (left < j)
        quickSort(arr, left, j);
    if (i < right)
        quickSort(arr, i, right);
}

void global_tests(INPUT &input, P &p, bool global[][ROWMAX], int n_gl, int current_row)
{
    for (int i = 0; i < n_gl; i++)
    {
        int x = input.gl[i];//1-Bonf;2-fish;3-chi;4-sims
        if (x == 1)//gl_bonf
        {
            global[0][current_row] = gl_bonf(p.value_sort, input.col, input.alpha);
        }
        else if (x == 2) //gl_chi
        {
            global[1][current_row] = gl_chi(p.data, input.col, input.alpha, input.mu0, input.sig0);
        }
        else if (x == 3) //gl_fish
        {
            global[2][current_row] = gl_fish(p.value_sort, input.col, input.alpha);
        }
        else //gl_sims
        {
            global[3][current_row] = gl_sims(p.value_sort, input.col, input.alpha);
        }
    }
}

void multiple_tests(INPUT &input, P &p, double multiple[][ROWMAX], int n_mul, int current_row)
{
    for (int i = 0; i < n_mul; i++)
    {
        int x = input.mul[i];//1-Bonf;2-holm;3-hoch;4-homm;5-BH
        if (x == 1)//mul_bonf
        {
            multiple[0][current_row] = mul_bonf(p.value_sort, input.col, input.alpha);
        }
        else if (x == 2) //mul_holm
        {
            multiple[1][current_row] = mul_holm(p.value_sort, input.col, input.alpha);
        }
        else if (x == 3) //mul_hoch
        {
            multiple[2][current_row] = mul_hoch(p.value_sort, input.col, input.alpha);
        }
        else if (x == 4) //mul_homm
        {
            multiple[3][current_row] = mul_homm(p.value_sort, input.col, input.alpha);
        }
        else //mul_BH
        {
            multiple[4][current_row] = mul_BH(p.value_sort, input.col, input.alpha);
        }
    }
}

void evaluate_mul(double confusion_mul[][MULNEED],
    P &p, double multiple[][ROWMAX], INPUT &input, int n_mul, int currentRow)
{
    // document the confusion matrix
    for (int i = 0; i < n_mul; i++)
    {
        double temp[4] = { 0 }; //temp confusion matrix
        for (int j = 0; j < input.col; j++) // compare element in one row
            {
            if ( (p.identity[j] == false) && (p.value[j] > multiple[i][currentRow]) )
            {
                temp[0] = temp[0] + 1; // assume = 0 , result=0
            }
            else if ( (p.identity[j] == true) && (p.value[j] > multiple[i][currentRow]) )
            {
                temp[1] = temp[1] + 1; // assume = 1 , result=0
            }
            else if ((p.identity[j] == false) && (p.value[j] <= multiple[i][currentRow]))
            {
                temp[2] = temp[2] + 1; // assume = 0 , result=1
            }
            else  //else if ((p.identity[j] == true) && (p.value[j] <= multiple[i][currentRow]))
            {
                temp[3] = temp[3] + 1; // assume = 1 , result=1
            }
            }
        // get FWER,FDR,power
        if (temp[2] >= 1)
        {
            confusion_mul[i][0] = confusion_mul[i][0] + 1;//V
        }
        if ((temp[2] + temp[3]) >= 1) // denominator cannot be zero
        {
            confusion_mul[i][1] = confusion_mul[i][1] + (temp[2]+0.0) / (temp[2] + temp[3]);//FDP
        }
        if ((temp[1] + temp[3]) >=1 ) // denominator cannot be zero
        {
            confusion_mul[i][2] = confusion_mul[i][3] + (temp[3] + 0.0) / (temp[1] + temp[3]); //power
        }
        else // T+S>=1
        {
            confusion_mul[i][2] = confusion_mul[i][3] + 1; //power
        }
        //cout << "mul matrix " << i+1 << " method " << temp[0] << ' ' << temp[1] << ' ' << temp[2] << ' ' << temp[3] << endl;
    }
}

void evaluate_gl(double confusion_gl[][GLNEED], P &p, bool global[][ROWMAX], INPUT &input, int n_gl)
{
    // document the confusion matrix
    for (int i = 0; i < n_gl; i++)
    {
        double temp[4] = { 0 }; //temp confusion matrix
        for (int j = 0; j < input.row; j++)
        {
            if ((input.rejected[j] == 0) && (global[i][j] == false))
            {
                temp[0] = temp[0] + 1; // assume = 0 , result=0
            }
            else if ((input.rejected[j] > 0) && (global[i][j] == false))
            {
                temp[1] = temp[1] + 1; // assume = 1 , result=0
            }
            else if ((input.rejected[j] == 0) && (global[i][j] == true))
            {
                temp[2] = temp[2] + 1; // assume = 0 , result=1
            }
            //else if ((input.rejected[j] > 0) && (global[i][j] == true))
            else
            {
                temp[3] = temp[3] + 1; // assume = 1 , result=1
            }
        }
        // get alpha, beta, power
        if ( (temp[0] + temp[2]) >=1 ) // denominator cannot be zero
        {
            confusion_gl[i][0] = (temp[2]+0.0) / (temp[0] + temp[2]); //alpha
        }
        else
        {
            confusion_gl[i][0] = 0;
        }
        if ((temp[1] + temp[3]) >= 1) // denominator cannot be zero
        {
            confusion_gl[i][1] = (temp[3] + 0.0) / (temp[1] + temp[3]); //power
        }
        else
        {
            confusion_gl[i][1] = 1;  //power
        }
        //cout << "gl matrix " << i + 1 << " method " << temp[0] << ' ' << temp[1] << ' ' << temp[2] << ' ' << temp[3] << endl;
    }
}



void show_p_value(int col, P &p, int dec1, int dec2, int dec3)// check pvalue
{
    if (dec1 == 1)
    {
        cout << "pdata" << ' ';
        for (int x = 0; x < col; x++)
        {
            cout << p.data[x] << ' ';
        }
        cout << endl;
    }
    if (dec2 == 1)
    {
        cout << "pidentity" << ' ';
        for (int x = 0; x < col; x++)
        {
            cout << p.identity[x] << ' ';
        }
        cout << endl;
    }
    if (dec3 == 1)
    {
        cout << "pvalue" << ' ';
        for (int x = 0; x < col; x++)
        {
            cout << p.value[x] << ' ';
        }
        cout << endl;
    }
}

void show_global_multiple(int n_gl, int n_mul, int row, bool global[][ROWMAX], double multiple[][ROWMAX], int dec1, int dec2)
{
    if (dec1 == 1)
    {
        cout << "global matrix" << endl;
        for (int x = 0; x < n_gl; x++)
        {
            for (int y = 0; y < row; y++)
            {
                cout << global[x][y] << ' ';
            }
            cout << endl;
        }
    }
    if (dec2 == 1)
    {
        cout << "multiple matrix" << endl;
        for (int x = 0; x < n_mul; x++)
        {
            for (int y = 0; y < row; y++)
            {
                cout << multiple[x][y] << ' ';
            }
            cout << endl;
        }
    }
}

void show_evaluation(double confusion_gl[][GLNEED], double confusion_mul[][MULNEED], int n_gl, int n_mul, int row, int dec1, int dec2)
{
    if (dec1 == 1)
    {
        for (int i = 0; i < n_gl; i++)
        {
            cout << "gl method " << i+1 << " alpha " << confusion_gl[i][0]
                << " power " << confusion_gl[i][1] << endl;
        }
    }
    if (dec2 == 1)
    {
        for (int i = 0; i < n_mul; i++)
        {
            cout << "mul method " << i+1 << " FWER " << (confusion_mul[i][0]+0.0)/row << " FDR " << (confusion_mul[i][1]+0.0)/row
                << " power " << (confusion_mul[i][2] + 0.0) / row << endl;
        }
    }
}
