#include <memory.h>
#include <stdio.h>

#include "pso_frame.h"
//#include "cubic_spline_with_curvature_constraint.h"
#include "basic_support.h"
//
//#include "spline.h"

// 示例适应度函数: Rastrigin函数
double rastrigin(double* x, int dim) {
    double sum = 0.0;
    for (int i = 0; i < dim; i++) {
        sum += x[i] * x[i] - 10 * cos(2 * PI * x[i]);
    }
    return 10 * dim + sum;
}

// 示例适应度函数: Sphere函数
double sphere(double* x, int dim) {
    double sum = 0.0;
    for (int i = 0; i < dim; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

// BestFitnessSpline
double fitness_spline(double* x, int dim) {
    double fitness = 0.0;
    
    return fitness;
}

#if 1
int main() {
    // 配置搜索空间
    int dim = 2;
    double lower_bounds[] = { -5.12, -5.12 }; // Rastrigin函数的搜索范围
    double upper_bounds[] = { 5.12, 5.12 };
    double velocity_range[] = { -100.0, 100.0 };

    // 创建PSO配置
    PSOConfig config = create_pso_config(
        dim,            // 问题维度
        30,             // 粒子数量
        1000,           // 最大迭代次数
        0.5,            // 惯性权重 w
        2.0,            // 个体学习因子 c1
        2.0,            // 群体学习因子 c2
        lower_bounds,   // 搜索下界
        upper_bounds,   // 搜索上界
        velocity_range
    );

    // 运行PSO优化器
    //run_pso(&config, rastrigin); // 传入rastrigin函数指针
    run_pso(&config, sphere); // 可替换为sphere函数

    return 0;
}
#endif

#if 0
// ================= 示例用法 =================
int main() { // cubic_spline_with....h
    // 创建路径点 (x,y)
    Matrix points = matrix_create(2, 3);
    double point_data[6] = { -300.588989, -21.019232,
                           -350.0, -30.0,
                           -453.451843, -49.921623 };
    //points.data = point_data;
    memcpy(points.data, point_data, sizeof(*point_data) * 6);

    // 设置终点切线方向
    double end_tangent[2] = { -1.0, 0.0 };

    // 生成带边界条件的三次样条
    Spline sp = clamped_clamped_spline(&points, end_tangent, 0.35);

    // 在参数空间均匀采样
    const int num_samples = 10;
    for (int i = 0; i <= num_samples; i++) {
        double t = (double)i / num_samples;

        // 计算位置
        double x = spline_eval(&sp, t, 0);
        double y = spline_eval(&sp, t, 0);

        // 计算曲率
        double curvature[2];
        spline_curvature(&sp, t, curvature);

        printf("t=%.2f: (%.3f, %.3f) | Curvature: %.6f, Change: %.6f\n",
            t, x, y, curvature[0], curvature[1]);
    }

    // 清理资源
    free(sp.breaks);
    free(sp.coefs);
    matrix_free(&points);

    return 0;
}
#endif

#if 0
// 主函数 - 测试代码
int main() { // spline.h
    // 参数设置
    int n = 3; // 节点数
    double u[] = { 0, 0.5, 1 }; // 参数节点
    int dim = 2; // 向量维度

    // 输入点：2维点 (0,0), (0.5,0.5), (1,1)
    double points[] = {
        0, 0,    // 点1 (x=0, y=0)
        0.4, 0.8, // 点2 (x=0.5, y=0.5)
        1, 1      // 点3 (x=1, y=1)
    };

    // 边界条件: 
    // [2,1] - 左端二阶导数, 右端一阶导数
    int bctype[] = { 2, 1 };
    // 边界值: 
    // 左端二阶导数(x=0, y=0), 右端一阶导数(x=1, y=0)
    double bcval[] = {
        0, 0, // 左边界 (dim0=x, dim1=y)
        1, 0  // 右边界 (x方向斜率=1, y方向斜率=0)
    };

    // 创建样条
    Spline3D* spline = csape_c(n, u, points, dim, bctype, bcval);

    if (spline == NULL) {
        return EXIT_FAILURE;
    }

    // 在[0,1]区间内采样计算
    int num_samples = 100;
    double result[2]; // 存储二维结果

    printf("参数t   x坐标     y坐标\n");
    printf("------------------------\n");
    for (int i = 0; i <= num_samples; i++) {
        double t = i / (double)num_samples;
        spline_eval(spline, t, result);
        printf("%.2f   %.6f  %.6f\n", t, result[0], result[1]);
    }

    // 释放内存
    free_spline(spline);

    return EXIT_SUCCESS;
}
#endif
