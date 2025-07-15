#include "spline.h"
#include <stdio.h>
#include <string.h>
#include <float.h>

// 创建ND点
PointND* pointnd_create(int dim) {
    PointND *point = (PointND*)malloc(sizeof(PointND));
    point->dim = dim;
    point->coords = (double*)malloc(dim * sizeof(double));
    return point;
}

// 从数组创建ND点
PointND* pointnd_from_array(double *coords, int dim) {
    PointND *point = pointnd_create(dim);
    memcpy(point->coords, coords, dim * sizeof(double));
    return point;
}

// 从2D点数组创建ND点数组
PointND* pointnd_from_point2d(Point2D *points, int num_points) {
    PointND *nd_points = (PointND*)malloc(num_points * sizeof(PointND));
    
    for (int i = 0; i < num_points; i++) {
        nd_points[i].dim = 2;
        nd_points[i].coords = (double*)malloc(2 * sizeof(double));
        nd_points[i].coords[0] = points[i].x;
        nd_points[i].coords[1] = points[i].y;
    }
    
    return nd_points;
}

// 释放ND点
void pointnd_free(PointND *point) {
    if (point) {
        if (point->coords) free(point->coords);
        free(point);
    }
}

// 弦长参数化 (ND版本)
double* chordal_parameterization_nd(PointND *points, int num_points) {
    double *t = (double*)malloc(num_points * sizeof(double));
    t[0] = 0.0;
    
    double total_length = 0.0;
    for (int i = 1; i < num_points; i++) {
        double sum_squared = 0.0;
        for (int d = 0; d < points[i].dim; d++) {
            double diff = points[i].coords[d] - points[i-1].coords[d];
            sum_squared += diff * diff;
        }
        double segment_length = sqrt(sum_squared);
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

// 创建多维样条插值
SplinePP* spline_create_nd(PointND *points, int num_points, PointND *end_tangent, double factor_end) {
    if (num_points < 2) return NULL;
    int dim = points[0].dim;
    
    // 计算路径长度用于比例缩放终点切线
    double path_length = 0.0;
    double sum_squared = 0.0;
    for (int d = 0; d < dim; d++) {
        double diff = points[num_points-1].coords[d] - points[0].coords[d];
        sum_squared += diff * diff;
    }
    path_length = sqrt(sum_squared);
    
    double extension_factor_end = factor_end * path_length;
    
    // 归一化终点切线向量
    double et_norm = 0.0;
    for (int d = 0; d < dim; d++) {
        et_norm += end_tangent->coords[d] * end_tangent->coords[d];
    }
    et_norm = sqrt(et_norm);
    
    if (et_norm > 1e-10) {
        for (int d = 0; d < dim; d++) {
            end_tangent->coords[d] = end_tangent->coords[d] / et_norm * extension_factor_end;
        }
    } else {
        // 默认切线
        for (int d = 0; d < dim; d++) {
            end_tangent->coords[d] = d == 0 ? extension_factor_end : 0;
        }
    }
    
    // 创建参数化向量
    double *t = chordal_parameterization_nd(points, num_points);
    
    // 创建并初始化样条结构
    SplinePP *pp = (SplinePP*)malloc(sizeof(SplinePP));
    pp->pieces = num_points - 1;
    pp->dim = dim;  
    pp->order = 4;  // 三次样条
    pp->t = t;
    
    // 分配系数矩阵内存
    pp->coefs = (double**)malloc(pp->dim * sizeof(double*));
    for (int i = 0; i < pp->dim; i++) {
        pp->coefs[i] = (double*)malloc(pp->pieces * pp->order * sizeof(double));
    }
    
    // 使用三对角矩阵算法计算三次样条插值
    for (int d = 0; d < pp->dim; d++) {
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
            y[i] = points[i].coords[d];
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
        double tangent_d = end_tangent->coords[d] / (t[num_points-1] - t[num_points-2]);
        alpha[num_points-2] = 3.0/h[num_points-2] * tangent_d - 3.0/h[num_points-2] * (y[num_points-1]-y[num_points-2])/h[num_points-2];
        l[num_points-2] = 2 * h[num_points-2];
        z[num_points-2] = alpha[num_points-2] / l[num_points-2];
        
        // 反向替换求解c,b,d系数
        for (int j = num_points-2; j >= 0; j--) {
            c[j] = z[j] - mu[j] * c[j+1];
            b[j] = (y[j+1]-y[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3.0;
            d[j] = (c[j+1]-c[j])/(3.0*h[j]);
        }
        
        // 填充系数矩阵
        for (int i = 0; i < pp->pieces; i++) {
            pp->coefs[d][i*4 + 0] = y[i];                // 常数项
            pp->coefs[d][i*4 + 1] = b[i];                // 一次项
            pp->coefs[d][i*4 + 2] = c[i];                // 二次项
            pp->coefs[d][i*4 + 3] = d[i];                // 三次项
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

// 兼容旧接口的2D样条创建
SplinePP* spline_create(Point2D *points, int num_points, Point2D end_tangent, double factor_end) {
    // 转换为ND点
    PointND *nd_points = pointnd_from_point2d(points, num_points);
    
    // 创建终点切线
    PointND end_tangent_nd;
    end_tangent_nd.dim = 2;
    end_tangent_nd.coords = (double*)malloc(2 * sizeof(double));
    end_tangent_nd.coords[0] = end_tangent.x;
    end_tangent_nd.coords[1] = end_tangent.y;
    
    // 创建样条
    SplinePP *result = spline_create_nd(nd_points, num_points, &end_tangent_nd, factor_end);
    
    // 清理
    free(end_tangent_nd.coords);
    for (int i = 0; i < num_points; i++) {
        free(nd_points[i].coords);
    }
    free(nd_points);
    
    return result;
}

// 计算样条在t处的值
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

// 注意：曲率计算仅适用于2D和3D
void spline_compute_curvature(SplinePP *sp, double *t, int n_samples, double *curvature, double *max_curv) {
    if (sp->dim != 2 && sp->dim != 3) {
        fprintf(stderr, "Warning: Curvature computation only supported for 2D/3D curves\n");
        return;
    }
    
    SplinePP *sp_der1 = spline_derivative(sp, 1);
    SplinePP *sp_der2 = spline_derivative(sp, 2);
    
    *max_curv = 0.0;
    
    for (int k = 0; k < n_samples; k++) {
        double *dx_dt = spline_evaluate(sp_der1, t[k]);
        double *ddx_dt = spline_evaluate(sp_der2, t[k]);
        
        double numerator = 0.0;
        double denominator_term = 0.0;
        
        if (sp->dim == 2) {
            // 2D曲率计算公式
            numerator = dx_dt[0] * ddx_dt[1] - dx_dt[1] * ddx_dt[0];
            denominator_term = dx_dt[0]*dx_dt[0] + dx_dt[1]*dx_dt[1];
        } else if (sp->dim == 3) {
            // 3D曲率计算公式 (Frenet-Serret公式)
            // |r' × r''|/|r'|^3
            double cross_x = dx_dt[1]*ddx_dt[2] - dx_dt[2]*ddx_dt[1];
            double cross_y = dx_dt[2]*ddx_dt[0] - dx_dt[0]*ddx_dt[2];
            double cross_z = dx_dt[0]*ddx_dt[1] - dx_dt[1]*ddx_dt[0];
            numerator = sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
            denominator_term = dx_dt[0]*dx_dt[0] + dx_dt[1]*dx_dt[1] + dx_dt[2]*dx_dt[2];
        }
        
        double denominator = pow(denominator_term, 1.5);
        
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

// 其他函数保持不变...

// 多维版本的曲线弧长采样
void sample_curve_by_arc_length_nd(double **coords, int n_points, int dim, double *s_list, int n_samples, PointND *result) {
    // 计算累积弧长
    double *segment_lengths = (double*)malloc((n_points-1) * sizeof(double));
    double *s_values = (double*)malloc(n_points * sizeof(double));
    
    s_values[0] = 0.0;
    for (int i = 0; i < n_points - 1; i++) {
        double sum_squared = 0.0;
        for (int d = 0; d < dim; d++) {
            double diff = coords[d][i+1] - coords[d][i];
            sum_squared += diff * diff;
        }
        segment_lengths[i] = sqrt(sum_squared);
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
        // 找到s_list[i]所在区间
        int seg = 0;
        while (seg < n_points-1 && s_list[i] > s_values[seg+1]) seg++;
        
        // 计算插值系数
        double alpha = 0;
        if (segment_lengths[seg] > 1e-10) {
            alpha = (s_list[i] - s_values[seg]) / segment_lengths[seg];
        }
        
        // 插值得到对应点
        for (int d = 0; d < dim; d++) {
            result[i].coords[d] = coords[d][seg] + alpha * (coords[d][seg+1] - coords[d][seg]);
        }
    }
    
    free(segment_lengths);
    free(s_values);
}

// 释放样条曲线资源
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
