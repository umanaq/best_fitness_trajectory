#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "spline.h"
#include "pso.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "spline.h"
#include "pso.h"

// 测试多维样条功能
void test_nd_spline() {
    printf("\n===== Testing Multi-dimensional Spline =====\n");
    
    // 创建3D点
    PointND points[4];
    for (int i = 0; i < 4; i++) {
        points[i].dim = 3;
        points[i].coords = (double*)malloc(3 * sizeof(double));
    }
    
    // 设置3D坐标
    points[0].coords[0] = 0.0;  points[0].coords[1] = 0.0;  points[0].coords[2] = 0.0;
    points[1].coords[0] = 10.0; points[1].coords[1] = 20.0; points[1].coords[2] = 5.0;
    points[2].coords[0] = 20.0; points[2].coords[1] = 10.0; points[2].coords[2] = 15.0;
    points[3].coords[0] = 30.0; points[3].coords[1] = 0.0;  points[3].coords[2] = 0.0;
    
    // 设置终点切线
    PointND end_tangent;
    end_tangent.dim = 3;
    end_tangent.coords = (double*)malloc(3 * sizeof(double));
    end_tangent.coords[0] = 1.0;
    end_tangent.coords[1] = 0.0;
    end_tangent.coords[2] = -0.5;
    
    // 创建3D样条
    SplinePP *spline3d = spline_create_nd(points, 4, &end_tangent, 0.35);
    
    // 计算并打印样条点
    printf("3D样条点:\n");
    for (int i = 0; i <= 10; i++) {
        double t = i / 10.0;
        double *point = spline_evaluate(spline3d, t);
        printf("t=%.1f: (%.2f, %.2f, %.2f)\n", 
               t, point[0], point[1], point[2]);
        free(point);
    }
    
    // 清理资源
    for (int i = 0; i < 4; i++) {
        free(points[i].coords);
    }
    free(end_tangent.coords);
    spline_free(spline3d);
}

int main() {
    // 设置随机数种子
    srand(2225);
    
    // 定义起点和终点
    Point2D start_pos = {-300.588989, -21.019232};
    Point2D end_pos = {-453.451843, -49.921623};
    Point2D global_end_pos = {-989.307068, -360.078033};
    
    // 计算终点切线方向
    double dx = global_end_pos.x - end_pos.x;
    double dy = global_end_pos.y - end_pos.y;
    double dist = sqrt(dx*dx + dy*dy);
    Point2D desired_end_tangent = {dx/dist, dy/dist};
    
    // 定义障碍物环境
    Environment env;
    env.num_obstacles = 1;
    env.obstacles = (Obstacle*)malloc(env.num_obstacles * sizeof(Obstacle));
    env.obstacles[0].x = 0.0;
    env.obstacles[0].y = 0.0;
    env.obstacles[0].radius = 0.01;
    
    // 配置PSO参数
    PSOConfig config;
    config.x_min = -1200;
    config.x_max = 1200;
    config.y_min = -1200;
    config.y_max = 1200;
    config.num_particles = 100;
    config.max_iter = 100;
    config.num_waypoints = 1;  // 中间路径点数量
    config.w = 0.7;            // 惯性权重
    config.c1 = 0.5;           // 认知因子
    config.c2 = 0.5;           // 社会因子
    config.max_velocity = 100.0;
    config.gamma = 1.0;        // 曲率惩罚系数
    config.beta = 1.0;         // 碰撞惩罚系数
    config.alpha = 1e6;        // 曲率约束系数
    config.tan2 = 1e4;         // 终点曲率惩罚系数
    config.trend_factor = 1e8; // 曲率趋势惩罚系数
    config.path_weight = 1e0;  // 路径长度权重
    config.factor_end = 0.35;  // 终点切线系数
    config.t_endzone_percent = 0.5; // 终点区域百分比
    config.Kmax = 0.015;       // 最大允许曲率
    config.d = 0.0;            // 安全距离
    
    // 初始化PSO
    PSO *pso = pso_init(start_pos, end_pos, desired_end_tangent, &env, config);
    
    // 运行PSO优化
    pso_run(pso, start_pos, end_pos, desired_end_tangent, &env);
    
    // 获取最佳路径
    int path_size;
    double **best_path = pso_get_best_path(pso, &path_size);
    SplinePP *best_spline = pso_get_best_spline(pso);
    
    // 输出最佳路径分析
    printf("\n===== Path Analysis Results =====\n");
    
    // 计算路径长度
    double path_length = 0.0;
    for (int k = 1; k < path_size; k++) {
        double dx = best_path[0][k] - best_path[0][k-1];
        double dy = best_path[1][k] - best_path[1][k-1];
        path_length += sqrt(dx*dx + dy*dy);
    }
    
    // 计算曲率
    double *t_samples = linspace(0, 1, path_size);
    double *curvature = (double*)malloc(path_size * sizeof(double));
    double max_curv = 0.0;
    
    spline_compute_curvature(best_spline, t_samples, path_size, curvature, &max_curv);
    
    // 计算曲率半径
    double min_curvature_radius = DBL_MAX;
    double sum_curvature_radius = 0.0;
    int count_valid = 0;
    
    for (int i = 0; i < path_size; i++) {
        if (curvature[i] > 1e-6) {
            double radius = 1.0 / curvature[i];
            min_curvature_radius = fmin(min_curvature_radius, radius);
            sum_curvature_radius += radius;
            count_valid++;
        }
    }
    
    double mean_curvature_radius = sum_curvature_radius / (count_valid > 0 ? count_valid : 1);
    
    // 输出路径分析结果
    printf("Total path length: %.2f mm\n", path_length);
    printf("Maximum curvature: %.4f\n", max_curv);
    printf("Minimum curvature radius: %.2f mm\n", min_curvature_radius);
    printf("Average curvature radius: %.2f mm\n", mean_curvature_radius);
    
    // 输出一些路径点作为示例
    printf("\n=== Sample Path Points ===\n");
    int num_samples_to_print = path_size > 10 ? 10 : path_size;
    for (int i = 0; i < num_samples_to_print; i++) {
        int idx = i * (path_size / num_samples_to_print);
        if (idx < path_size) {
            printf("Point %d: (%.3f, %.3f)\n", idx, best_path[0][idx], best_path[1][idx]);
        }
    }
    
    // 清理资源
    free(t_samples);
    free(curvature);
    pso_free(pso);
    free(env.obstacles);
    
    return 0;
}
