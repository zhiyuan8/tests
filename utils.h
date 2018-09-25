#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include "pre_define.h"

void get_size(INPUT &input, int &n_gl, int &n_mul, int &n_effect, int &n_ratio2);
void Initialize_size(INPUT &input, int n_gl, int n_mul);
void Initialize_confusion(double confusion_gl[][GLNEED], double confusion_mul[][MULNEED], int n_gl, int n_mul);
int round_double(double number);
void define_a(double mu0, double mu1, int &a);
void generate_ground_truth(INPUT &input);
double cdf_approx(double x);
void quickSort(double arr[], int left, int right);
void global_tests(INPUT &input, P &p, bool global[][ROWMAX], int n_gl, int current_row);
void multiple_tests(INPUT &input, P &p, double multiple[][ROWMAX], int n_mul, int current_row);
void evaluate_mul(double confusion_mul[][MULNEED],
    P &p, double multiple[][ROWMAX], INPUT &input, int n_mul, int currentRow);
void evaluate_gl(double confusion_gl[][GLNEED], P &p, bool global[][ROWMAX], INPUT &input, int n_gl);

void show_p_value(int col, P &p, int dec1, int dec2, int dec3);
void show_global_multiple(int n_gl, int n_mul, int row, bool global[][ROWMAX], double multiple[][ROWMAX], int dec1, int dec2);
void show_evaluation(double confusion_gl[][GLNEED], double confusion_mul[][MULNEED], int n_gl, int n_mul, int row, int dec1, int dec2);
#endif
