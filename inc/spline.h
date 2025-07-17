#ifndef _SPLINE_H_
#define _SPLINE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basic_support.h"


// 样条结构体
typedef struct {
    int dim;            // 数据维度
    int num_breaks;     // 断点数
    int num_segments;   // 分段数 (num_breaks - 1)
    double* breaks;     // 断点数组（自变量）
    double* coefs;      // 系数数组：dim × num_segments × 4
} Spline;

// 三对角方程组求解 (Thomas算法)
void thomas_algorithm(int n, double* a, double* b, double* c, double* d, double* x);
double* compute_h(int n, double* u);
Spline* csape_c(int n, double* u, double* points, int dim, int* bctype, double* bcval);
double spline_deriv(Spline* spline, double t, int dim, int order);
Point spline_deriv_xy(Spline* spline, double t, int order);
int validate_boundary_conditions(Spline* spline, int* bctype, double* bcval);
Point spline_ppval(Spline* spline, double t);
double spline_eval(Spline* spline, double t, int dim);
void spline_free(Spline* spline);
void print_spline(Spline* spline);

#endif  // _SPLINE_H_
