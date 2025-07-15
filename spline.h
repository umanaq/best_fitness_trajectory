#ifndef SPLINE_H
#define SPLINE_H

#include <stdlib.h>
#include <math.h>

// 多维点结构
typedef struct {
    double *coords;   // 坐标数组
    int dim;          // 维度
} PointND;

// 兼容旧代码的2D点结构
typedef struct {
    double x;
    double y;
} Point2D;

typedef struct {
    double *t;       // 参数化值
    double **coefs;  // 系数矩阵
    int pieces;      // 样条段数
    int dim;         // 维度
    int order;       // 阶数 (通常为4，表示三次样条)
} SplinePP;

// 创建ND点
PointND* pointnd_create(int dim);

// 从数组创建ND点
PointND* pointnd_from_array(double *coords, int dim);

// 从2D点数组创建ND点数组
PointND* pointnd_from_point2d(Point2D *points, int num_points);

// 释放ND点
void pointnd_free(PointND *point);

// 使用Clamped边界条件创建多维三次样条曲线
SplinePP* spline_create_nd(PointND *points, int num_points, PointND *end_tangent, double factor_end);

// 兼容性函数：从2D点创建样条
SplinePP* spline_create(Point2D *points, int num_points, Point2D end_tangent, double factor_end);

// 计算样条在参数t处的点 (返回dim维数组)
double* spline_evaluate(SplinePP *pp, double t);

// 计算样条的n阶导数
SplinePP* spline_derivative(SplinePP *pp, int n);

// 弦长参数化
double* chordal_parameterization_nd(PointND *points, int num_points);

// 计算样条曲线在给定参数集合处的曲率 (仅适用于2D/3D)
void spline_compute_curvature(SplinePP *sp, double *t, int n_samples, double *curvature, double *max_curv);

// 计算终点区域的曲率适应度
double spline_calculate_fitness_end(SplinePP *sp, double t_endzone_percent, double *trend_penalty);

// 创建等间距参数数组
double* linspace(double start, double end, int n);

// 采样曲线上的点
void sample_curve_by_arc_length_nd(double **coords, int n_points, int dim, double *s_list, int n_samples, PointND *result);

// 释放样条曲线资源
void spline_free(SplinePP *pp);

#endif /* SPLINE_H */
