#include <iostream>
#include <chrono>  // for high_resolution_clock
#include <vector>
#include <time.h> //show calculation time

#include "pre_define.h"
#include "main.h"
#include "utils.h"
//#include "io.h" //read in and output
#include "generate.h"
#include "global_testing.h"
#include "multiple_testing.h"

using namespace std;

int main()
{
        cout<<"OK I have sth"<<endl;
        // This part should be user-defined later
        input.gl = {1,4};//1-Bonf;2-fish;3-chi;4-sims
        input.mul = {1,2,3,4,5};//1-Bonf;2-holm;3-hoch;4-homm;5-BH
        input.MU1 = {0.1,0.5,1,3};
        input.RATIO2 = {0.1,0.2,0.4};
        input.rep_time = 2;
        input.row =50;
        input.col = 100;
        input.ratio1 = 0.5;

        startTime = clock(); // begin documenting time
        get_size(input, n_gl, n_mul, n_effect, n_ratio2); // get n_gl, n_mul, n_effect, n_ratio2
        Initialize_size(input, n_gl, n_mul); // initialize input struct
        if ((GLMAX < n_gl) || (MULMAX < n_mul) || (ROWMAX < input.row) || (COLMAX < input.col))
    {
        cout << "Check your GLMAX MULMAX ROWMAX COLMAX size in pre_define.h";
                return 0;//check array size meets our requirement
    }
    for (int r = 0; r < n_effect; r++)
    {
        input.mu1 = input.MU1[r];
        define_a(input.mu0, input.mu1, a); //use area in left/right tail as p - value
        for (int s = 0; s < n_ratio2; s++)
        {
            input.ratio2 = input.RATIO2[s];
            for (int t = 0; t < input.rep_time; t++)
            {
                position++;
                input.seed = input.seed+(t + 1)*1000; // change seed in each loop
                                generate_ground_truth(input); // get input.rejected array
                                Initialize_confusion(confusion_gl, confusion_mul, n_gl, n_mul); // initialize confusion_mul and confusion_gl
                                for (int i = 0; i < input.row; i++)
                                {// A large for loop begins here, deal with data row by row
                                        start = std::chrono::high_resolution_clock::now();
                                        generate(i, p, input, a); // generate a row of pvalue
                                        finish = std::chrono::high_resolution_clock::now();
                                        elapsed = finish - start;
                                        generate_time = generate_time + elapsed.count();
                                        quickSort(p.value_sort, 0, input.col-1); // sort p value by my own function
                                        show_p_value(input.col, p, 1,1,1); // if 1, then show me p data/identity/value
                                        start = std::chrono::high_resolution_clock::now();
                                        global_tests(input, p, global, n_gl,i);//get 1/0 for ith row
                                        multiple_tests(input, p, multiple, n_mul, i);//get criteria for ith row
                                        finish = std::chrono::high_resolution_clock::now();
                                        elapsed = finish - start;
                                        testTime = testTime + elapsed.count();
                                        evaluate_mul(confusion_mul, p, multiple, input, n_mul, i); //get V,FDP,beta,power
                                }// The big for loop ends here
                                evaluate_gl(confusion_gl, p, global, input, n_gl); //get alpha,beta,power
                                show_global_multiple(n_gl, n_mul, input.row, global, multiple,1,1); // if 1, show me global/multiple matrix
                                show_evaluation(confusion_gl, confusion_mul,n_gl,n_mul,input.row,1,1); // if 1, shown me global/multiple confusion matrix
                cout << "finished case" << position << " mu1 "<< input.mu1 << " input.ratio2 " << input.ratio2 << " rep_time " << t+1 << endl;
            }// one case is finished
        }
    }
        endTime = clock();
    // show total time
        cout << "Generate p value time :" << generate_time << "s" << endl;
        cout << "Global/multiple testing time :" << testTime << "s" << endl;
        cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
        return 0;
}
