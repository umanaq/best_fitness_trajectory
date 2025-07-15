#include "spline.h"
#include <stdio.h>
#include <string.h>
#include <float.h>

// 弦长参数化
double* chordal_parameterization(Point2D *points, int num_points) {
    double *t = (double*)malloc((num_points) * sizeof(double));
    t[0] = 0.0;
    
    double total_length = 0.0;
    for (int i = 1; i < num_points; i++) {
        double dx = points[i].x - points[i-1].x;
        double dy = points[i].y - points[i-1].y;
        double segment_length = sqrt(dx*dx + dy*dy);
        total_length += segment_length;
        t[i] = total_length;
    }
    
    // 归一化参数
    if (total_length > 0) {
        for (int i = 0; i < num_points; i++) {
            t[i] /= total_length;
        }
    }
    
    return t;
}

SplinePP* spline_create(Point2D *points, int num_points, Point2D end_tangent, double factor_end) {
    double path_length = sqrt(pow(points[num_points-1].x - points[0].x, 2) + 
                            pow(points[num_points-1].y - points[0].y, 2));
    double extension_factor_end = factor_end * path_length;
    
    // 归一化终点切线向量
    double et_norm = sqrt(end_tangent.x * end_tangent.x + end_tangent.y * end_tangent.y);
    if (et_norm > 1e-10) {
        end_tangent.x = end_tangent.x / et_norm * extension_factor_end;
        end_tangent.y = end_tangent.y / et_norm * extension_factor_end;
    } else {
        end_tangent.x = extension_factor_end;
        end_tangent.y = 0;
    }
    
    // 创建参数化向量
    double *t = chordal_parameterization(points, num_points);
    
    // 创建并初始化样条结构
    SplinePP *pp = (SplinePP*)malloc(sizeof(SplinePP));
    pp->pieces = num_points - 1;
    pp->dim = 2;  // x和y坐标
    pp->order = 4; // 三次样条
    pp->t = t;
    
    // 分配系数矩阵内存
    pp->coefs = (double**)malloc(pp->dim * sizeof(double*));
    for (int i = 0; i < pp->dim; i++) {
        pp->coefs[i] = (double*)malloc(pp->pieces * pp->order * sizeof(double));
    }
    
    // 使用三对角矩阵算法计算三次样条插值
    for (int dim = 0; dim < pp->dim; dim++) {
        double *y = (double*)malloc(num_points * sizeof(double));
        double *h = (double*)malloc((num_points-1) * sizeof(double));
        double *alpha = (double*)malloc((num_points-1) * sizeof(double));
        double *l = (double*)malloc(num_points * sizeof(double));
        double *mu = (double*)malloc(num_points * sizeof(double));
        double *z = (double*)malloc(num_points * sizeof(double));
        double *c = (double*)malloc(num_points * sizeof(double));
        double *b = (double*)malloc((num_points-1) * sizeof(double));
        double *d = (double*)malloc((num_points-1) * sizeof(double));
        
        // 提取对应维度的坐标
        for (int i = 0; i < num_points; i++) {
            if (dim == 0) y[i] = points[i].x;
            else y[i] = points[i].y;
        }
        
        // 计算步长
        for (int i = 0; i < num_points-1; i++) {
            h[i] = t[i+1] - t[i];
        }
        
        // 计算alpha
        for (int i = 1; i < num_points-1; i++) {
            alpha[i] = 3.0/h[i] * (y[i+1]-y[i]) - 3.0/h[i-1] * (y[i]-y[i-1]);
        }
        
        // 三对角矩阵算法
        l[0] = 1.0;
        mu[0] = 0.0;
        z[0] = 0.0;
        
        for (int i = 1; i < num_points-1; i++) {
            l[i] = 2 * (t[i+1] - t[i-1]) - h[i-1] * mu[i-1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
        }
        
        l[num_points-1] = 1.0;
        z[num_points-1] = 0.0;
        c[num_points-1] = 0.0;
        
        // 在终点应用切线条件
        if (dim == 0) {
            double tangent_dx = end_tangent.x / (t[num_points-1] - t[num_points-2]);
            alpha[num_points-2] = 3.0/h[num_points-2] * tangent_dx - 3.0/h[num_points-2] * (y[num_points-1]-y[num_points-2])/h[num_points-2];
            l[num_points-2] = 2 * h[num_points-2];
            z[num_points-2] = alpha[num_points-2] / l[num_points-2];
        } else {
            double tangent_dy = end_tangent.y / (t[num_points-1] - t[num_points-2]);
            alpha[num_points-2] = 3.0/h[num_points-2] * tangent_dy - 3.0/h[num_points-2] * (y[num_points-1]-y[num_points-2])/h[num_points-2];
            l[num_points-2] = 2 * h[num_points-2];
            z[num_points-2] = alpha[num_points-2] / l[num_points-2];
        }
        
        // 反向替换求解c,b,d系数
        for (int j = num_points-2; j >= 0; j--) {
            c[j] = z[j] - mu[j] * c[j+1];
            b[j] = (y[j+1]-y[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3.0;
            d[j] = (c[j+1]-c[j])/(3.0*h[j]);
        }
        
        // 填充系数矩阵
        for (int i = 0; i < pp->pieces; i++) {
            pp->coefs[dim][i*4 + 0] = y[i];                // 常数项
            pp->coefs[dim][i*4 + 1] = b[i];                // 一次项
            pp->coefs[dim][i*4 + 2] = c[i];                // 二次项
            pp->coefs[dim][i*4 + 3] = d[i];                // 三次项
        }
        
        // 释放临时数组
        free(y);
        free(h);
        free(alpha);
        free(l);
        free(mu);
        free(z);
        free(c);
        free(b);
        free(d);
    }
    
    return pp;
}

double* spline_evaluate(SplinePP *pp, double t) {
    double *result = (double*)malloc(pp->dim * sizeof(double));
    
    // 找到t在哪个区间
    int i = 0;
    while (i < pp->pieces && t > pp->t[i+1]) i++;
    
    // 如果t超出范围，夹到边界
    if (i >= pp->pieces) i = pp->pieces - 1;
    if (t > pp->t[pp->pieces]) t = pp->t[pp->pieces];
    if (t < pp->t[0]) t = pp->t[0];
    
    // 计算相对于区间起点的t值
    double rel_t = (t - pp->t[i]);
    
    // 计算每个维度的值
    for (int dim = 0; dim < pp->dim; dim++) {
        double *c = &pp->coefs[dim][i * pp->order];
        result[dim] = c[0] + rel_t * (c[1] + rel_t * (c[2] + rel_t * c[3]));
    }
    
    return result;
}

SplinePP* spline_derivative(SplinePP *pp, int n) {
    if (n <= 0) return pp; // 零阶导数就是原函数
    
    SplinePP *der_pp = (SplinePP*)malloc(sizeof(SplinePP));
    der_pp->pieces = pp->pieces;
    der_pp->dim = pp->dim;
    der_pp->order = pp->order - n > 0 ? pp->order - n : 1; // 导数后阶数降低
    der_pp->t = pp->t; // 共享参数化
    
    // 分配系数矩阵内存
    der_pp->coefs = (double**)malloc(der_pp->dim * sizeof(double*));
    for (int i = 0; i < der_pp->dim; i++) {
        der_pp->coefs[i] = (double*)malloc(der_pp->pieces * der_pp->order * sizeof(double));
    }
    
    // 计算导数系数
    for (int dim = 0; dim < pp->dim; dim++) {
        for (int i = 0; i < pp->pieces; i++) {
            double h = pp->t[i+1] - pp->t[i];
            double *c = &pp->coefs[dim][i * pp->order];
            double *d = &der_pp->coefs[dim][i * der_pp->order];
            
            if (n == 1) {
                // 一阶导数
                d[0] = c[1];
                d[1] = 2 * c[2];
                d[2] = 3 * c[3];
            } else if (n == 2) {
                // 二阶导数
                d[0] = 2 * c[2];
                d[1] = 6 * c[3];
            } else {
                // 高阶导数
                if (n == 3 && pp->order >= 4) {
                    d[0] = 6 * c[3];
                } else {
                    d[0] = 0; // 三次样条的三阶以上导数为常数或0
                }
            }
        }
    }
    
    return der_pp;
}

void spline_compute_curvature(SplinePP *sp, double *t, int n_samples, double *curvature, double *max_curv) {
    SplinePP *sp_der1 = spline_derivative(sp, 1);
    SplinePP *sp_der2 = spline_derivative(sp, 2);
    
    *max_curv = 0.0;
    
    for (int k = 0; k < n_samples; k++) {
        double *dx_dt = spline_evaluate(sp_der1, t[k]);
        double *ddx_dt = spline_evaluate(sp_der2, t[k]);
        
        double numerator = dx_dt[0] * ddx_dt[1] - dx_dt[1] * ddx_dt[0];
        double denominator = pow(dx_dt[0]*dx_dt[0] + dx_dt[1]*dx_dt[1], 1.5);
        
        if (denominator > 1e-10) {
            curvature[k] = fabs(numerator) / denominator;
        } else {
            curvature[k] = 0.0;
        }
        
        if (curvature[k] > *max_curv) {
            *max_curv = curvature[k];
        }
        
        free(dx_dt);
        free(ddx_dt);
    }
    
    // 释放导数样条
    if (sp != sp_der1) spline_free(sp_der1);
    if (sp != sp_der2) spline_free(sp_der2);
}

double spline_calculate_fitness_end(SplinePP *sp, double t_endzone_percent, double *trend_penalty) {
    int n_samples = 100;
    double *t_endzone = linspace(t_endzone_percent, 1.0, n_samples);
    double *curvatures = (double*)malloc(n_samples * sizeof(double));
    double max_curv = 0.0;
    
    spline_compute_curvature(sp, t_endzone, n_samples, curvatures, &max_curv);
    
    int increasing_count = 0;
    double sum_deriv = 0.0;
    
    for (int i = 0; i < n_samples; i++) {
        if (i > 0) {
            double deriv = curvatures[i] - curvatures[i-1];
            if (deriv > 0) {
                increasing_count++;
            }
            sum_deriv += fabs(deriv);
        }
    }
    
    *trend_penalty = increasing_count + sum_deriv / (n_samples - 1);
    double curvature_end = max_curv + 0.1 * (*trend_penalty);
    
    free(t_endzone);
    free(curvatures);
    
    return curvature_end;
}

double* linspace(double start, double end, int n) {
    double *result = (double*)malloc(n * sizeof(double));
    double step = (end - start) / (n - 1);
    
    for (int i = 0; i < n; i++) {
        result[i] = start + i * step;
    }
    
    return result;
}

// 一维线性插值辅助函数
double interp1(double *x, double *y, int n, double xi) {
    // 简单线性插值
    if (xi <= x[0]) return y[0];
    if (xi >= x[n-1]) return y[n-1];
    
    // 找到xi所在区间
    int i = 0;
    while (i < n - 1 && xi > x[i+1]) i++;
    
    // 线性插值
    double t = (xi - x[i]) / (x[i+1] - x[i]);
    return y[i] + t * (y[i+1] - y[i]);
}

void sample_curve_by_arc_length(double *x, double *y, int n_points, double *s_list, int n_samples, Point2D *result) {
    // 计算累积弧长
    double *dx = (double*)malloc((n_points-1) * sizeof(double));
    double *dy = (double*)malloc((n_points-1) * sizeof(double));
    double *segment_lengths = (double*)malloc((n_points-1) * sizeof(double));
    double *s_values = (double*)malloc(n_points * sizeof(double));
    
    s_values[0] = 0.0;
    for (int i = 0; i < n_points - 1; i++) {
        dx[i] = x[i+1] - x[i];
        dy[i] = y[i+1] - y[i];
        segment_lengths[i] = sqrt(dx[i]*dx[i] + dy[i]*dy[i]);
        s_values[i+1] = s_values[i] + segment_lengths[i];
    }
    
    double total_length = s_values[n_points-1];
    
    // 检查s_list范围
    for (int i = 0; i < n_samples; i++) {
        if (s_list[i] < 0) s_list[i] = 0;
        if (s_list[i] > total_length - 1e-5) s_list[i] = total_length;
    }
    
    // 插值获取对应点
    for (int i = 0; i < n_samples; i++) {
        result[i].x = interp1(s_values, x, n_points, s_list[i]);
        result[i].y = interp1(s_values, y, n_points, s_list[i]);
    }
    
    free(dx);
    free(dy);
    free(segment_lengths);
    free(s_values);
}

void spline_free(SplinePP *pp) {
    if (pp == NULL) return;
    
    if (pp->t) free(pp->t);
    
    if (pp->coefs) {
        for (int i = 0; i < pp->dim; i++) {
            if (pp->coefs[i]) free(pp->coefs[i]);
        }
        free(pp->coefs);
    }
    
    free(pp);
}
