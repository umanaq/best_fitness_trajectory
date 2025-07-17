#include "spline.h"

// 三对角方程组求解 (Thomas算法)
void thomas_algorithm(int n, double *a, double *b, double *c, double *d, double *x)
{
    if (n < 1)
        return;

    double *cp = (double *)malloc(n * sizeof(double));
    double *dp = (double *)malloc(n * sizeof(double));

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (int i = 1; i < n; i++)
    {
        double m = 1.0 / (b[i] - a[i] * cp[i - 1]);
        cp[i] = c[i] * m;
        dp[i] = (d[i] - a[i] * dp[i - 1]) * m;
    }

    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    free(cp);
    free(dp);
}

// 计算步长h
double *compute_h(int n, double *u)
{
    double *h = (double *)malloc((n - 1) * sizeof(double));
    for (int i = 0; i < n - 1; i++)
    {
        h[i] = u[i + 1] - u[i];
        if (h[i] <= 0) {
            fprintf(stderr, "Breakpoints must be strictly increasing\n");
            free(h);
            return NULL;
        }
    }
    return h;
}

// 构造样条曲线
Spline *csape_c(int n, double *u, double *points, int dim, int *bctype, double *bcval)
{
    if (n < 2) {
        fprintf(stderr, "At least 2 breakpoints required\n");
        return NULL;
    }

    double *h = compute_h(n, u);
    Spline *spline = (Spline *)malloc(sizeof(Spline));

    spline->dim = dim;
    spline->num_breaks = n;
    spline->num_segments = n - 1;
    spline->breaks = (double *)malloc(n * sizeof(double));
    memcpy(spline->breaks, u, n * sizeof(double));

    // 分配系数数组: dim × segments × 4
    spline->coefs = (double *)calloc(dim * (n - 1) * 4, sizeof(double));
    memset(spline->coefs, 0, dim * (n - 1) * 4 * sizeof(double));
    
    for (int d = 0; d < dim; d++) {
        double *y = (double *)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) {
            y[i] = points[i * dim + d];
        }

        double *a = (double *)calloc(n, sizeof(double));     // 下对角线
        double *b = (double *)calloc(n, sizeof(double));     // 主对角线
        double *c = (double *)calloc(n, sizeof(double));     // 上对角线
        double *d_vec = (double *)calloc(n, sizeof(double)); // 右端项
        double *M = (double *)calloc(n, sizeof(double));     // 二阶导数

        // 边界条件
        if (bctype[0] == 1) {
            b[0] = 2.0;
            c[0] = 1.0;
            d_vec[0] = 6 * ((y[1] - y[0]) / h[0] - bcval[d * 2]) / h[0];
        } else if (bctype[0] == 2) {
            b[0] = 1.0;
            c[0] = 0.0;
            d_vec[0] = bcval[d * 2];
        }

        if (bctype[1] == 1) {
            a[n - 1] = 1.0;
            b[n - 1] = 2.0;
            d_vec[n - 1] = 6.0 * (bcval[d * 2 + 1] - (y[n - 1] - y[n - 2]) / h[n - 2]) / h[n - 2];
        } else if (bctype[1] == 2) {
            a[n - 1] = 0.0;
            b[n - 1] = 1.0;
            d_vec[n - 1] = bcval[d * 2 + 1];
        }

        // 内部节点方程
        for (int i = 1; i < n - 1; i++) {
            /*
            double mu = h[i - 1] / (h[i - 1] + h[i]);
            double lambda = h[i] / (h[i - 1] + h[i]);
            double delta_l = (y[i] - y[i - 1]) / h[i - 1];
            double delta_r = (y[i + 1] - y[i]) / h[i];
            double div_diff = 6 * (delta_r - delta_l) / (h[i - 1] + h[i]);

            a[i] = mu;
            b[i] = 2.0;
            c[i] = lambda;
            d_vec[i] = div_diff;
            */
            a[i] = h[i - 1];
            b[i] = 2.0 * (h[i - 1] + h[i]);
            c[i] = h[i];

            // 计算右侧项
            double delta_l = (y[i] - y[i - 1]) / h[i - 1];
            double delta_r = (y[i + 1] - y[i]) / h[i];
            d_vec[i] = 6.0 * (delta_r - delta_l);
        }

        // 求解三对角方程组
        thomas_algorithm(n, a, b, c, d_vec, M);

        // 计算样条系数
        for (int i = 0; i < n - 1; i++) {
            int idx = d * (n - 1) * 4 + i * 4;
            spline->coefs[idx] = y[i];
            spline->coefs[idx + 1] = (y[i + 1] - y[i]) / h[i] - (2 * M[i] + M[i + 1]) * h[i] / 6;
            spline->coefs[idx + 2] = M[i] / 2;
            spline->coefs[idx + 3] = (M[i + 1] - M[i]) / (6 * h[i]);
        }

        free(y);
        free(a);
        free(b);
        free(c);
        free(d_vec);
        free(M);
    }

    free(h);
    return spline;
}

// 计算样条导数
double spline_deriv(Spline* spline, double t, int dim, int order) {
    // 找到包含t的区间
    int segment = 0;
    while (segment < spline->num_segments - 1 && t >= spline->breaks[segment + 1]) {
        segment++;
    }
    if (segment >= spline->num_segments) segment = spline->num_segments - 1;

    // 计算局部参数u
    double u = t - spline->breaks[segment];

    // 获取当前维度的系数
    int idx = dim * spline->num_segments * 4 + segment * 4;
    double a = spline->coefs[idx];
    double b = spline->coefs[idx + 1];
    double c = spline->coefs[idx + 2];
    double d_coef = spline->coefs[idx + 3];

    switch (order) {
    case 0:  // 函数值
        return a + u * (b + u * (c + u * d_coef));

    case 1:  // 一阶导数
        return b + u * (2 * c + u * 3 * d_coef);

    case 2:  // 二阶导数
        return 2 * c + u * 6 * d_coef;

    default:
        fprintf(stderr, "Unsupported derivative order: %d\n", order);
        return 0.0;
    }
}

Point spline_deriv_xy(Spline* spline, double t, int order) {
    Point result = {0.0, 0.0};  // 初始化返回点
    
    // 检查order有效性（所有维度共享错误检查）
    if (order < 0 || order > 2) {
        fprintf(stderr, "Unsupported derivative order: %d\n", order);
        return result;
    }

    // 步骤1：统一查找t的区间（仅一次）
    int segment = 0;
    while (segment < spline->num_segments - 1 && t >= spline->breaks[segment + 1]) {
        segment++;
    }
    if (segment >= spline->num_segments) segment = spline->num_segments - 1;
    
    // 步骤2：计算局部参数（仅一次）
    double u = t - spline->breaks[segment];
    
    // 步骤3：计算两个维度的导数
    for (int dim = 0; dim < 2; dim++) {
        // 获取当前维度的系数起始索引
        int idx = dim * spline->num_segments * 4 + segment * 4;
        double a = spline->coefs[idx];
        double b = spline->coefs[idx + 1];
        double c = spline->coefs[idx + 2];
        double d_coef = spline->coefs[idx + 3];
        
        // 计算该维度的导数值
        double deriv_val;
        switch (order) {
            case 0:  // 函数值
                deriv_val = a + u * (b + u * (c + u * d_coef));
                break;
            case 1:  // 一阶导数
                deriv_val = b + u * (2 * c + u * 3 * d_coef);
                break;
            case 2:  // 二阶导数
                deriv_val = 2 * c + u * 6 * d_coef;
                break;
        }
        
        // 存储结果
        if (dim == 0) result.x = deriv_val;
        else result.y = deriv_val;
    }
    
    return result;
}

// 边界条件自检验
int validate_boundary_conditions(Spline* spline, int* bctype, double* bcval) {
    const double epsilon = 1e-5;
    int valid = 1;

    for (int dim = 0; dim < spline->dim; dim++) {
        // 获取左边界值
        double t_left = spline->breaks[0];
        int bc_index = dim * 2;  // 每维度有两个边界条件

        // 左边界检验
        if (bctype[0] == 2) {  // 二阶导数条件
            double second_deriv = spline_deriv(spline, t_left, dim, 2);
            double diff = fabs(second_deriv - bcval[bc_index]);

            if (diff > epsilon) {
                printf("VALIDATION FAILED [DIM %d, LEFT]: ", dim);
                printf("Second derivative (%.6f) != %.6f\n", second_deriv, bcval[bc_index]);
                valid = 0;
            } else {
                printf("VALID [DIM %d, LEFT]: ", dim);
                printf("Second derivative = %.6f (expected %.6f)\n", second_deriv, bcval[bc_index]);
            }
        }

        // 右边界检验
        double t_right = spline->breaks[spline->num_breaks - 1];

        if (bctype[1] == 1) {  // 一阶导数条件
            // 使用精确的导数计算
            double slope = spline_deriv(spline, t_right, dim, 1);
            double diff = fabs(slope - bcval[bc_index + 1]);

            if (diff > epsilon) {
                printf("VALIDATION FAILED [DIM %d, RIGHT]: ", dim);
                printf("Slope (%.6f) != %.6f\n", slope, bcval[bc_index + 1]);
                valid = 0;
            } else {
                printf("VALID [DIM %d, RIGHT]: ", dim);
                printf("Slope = %.6f (expected %.6f)\n", slope, bcval[bc_index + 1]);
            }
        }
    }

    return valid;
}

// 评估样条曲线
double spline_eval(Spline* spline, double t, int dim) {
    return spline_deriv(spline, t, dim, 0);
}

Point spline_ppval(Spline* spline, double t) {
    return spline_deriv_xy(spline, t, 0);
}

// 释放样条内存
void spline_free(Spline* spline) {
    if (spline) {
        if (spline->breaks) free(spline->breaks);
        if (spline->coefs) free(spline->coefs);
        free(spline);
    }
}

// 打印样条信息（调试用）
void print_spline(Spline* spline) {
    printf("Spline Information:\n");
    printf("  Dimensions: %d\n", spline->dim);
    printf("  Number of breakpoints: %d\n", spline->num_breaks);
    printf("  Number of segments: %d\n", spline->num_segments);

    printf("Breakpoints: ");
    for (int i = 0; i < spline->num_breaks; i++) {
        printf("%.2f ", spline->breaks[i]);
    }
    printf("\n");

    printf("Coefficients per dimension and segment:\n");
    for (int d = 0; d < spline->dim; d++) {
        printf("Dimension %d:\n", d);
        for (int seg = 0; seg < spline->num_segments; seg++) {
            int idx = d * spline->num_segments * 4 + seg * 4;
            printf("  Segment %d: a=%.6f, b=%.6f, c=%.6f, d=%.6f\n", seg,
                spline->coefs[idx], spline->coefs[idx + 1],
                spline->coefs[idx + 2], spline->coefs[idx + 3]);
        }
    }
}
