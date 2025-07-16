#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "spline.h"
#include "pso.h"
#include "speed_profile.h"

int main() {
    // 设置随机数种子
    srand(2225);
    
    // 定义起点和终点
    Point2D start_pos = {-300.588989, -21.019232};
    Point2D end_pos = {-453.451843, -49.921623};
    Point2D global_start_pos = {-300.588989, -21.019232};
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
    printf("\n===== 路径分析结果 =====\n");
    
    // 计算路径长度
    double path_length = 0.0;
    for (int k = 1; k < path_size; k++) {
        double dx = best_path[0][k] - best_path[0][k-1];
        double dy = best_path[1][k] - best_path[1][k-1];
        path_length += sqrt(dx*dx + dy*dy);
    }
    
    // 计算曲率
    int num_samples = 1000;
    double *t_samples = linspace(0, 1, num_samples);
    double *curvature = (double*)malloc(num_samples * sizeof(double));
    double max_curv = 0.0;
    
    spline_compute_curvature(best_spline, t_samples, num_samples, curvature, &max_curv);
    
    // 计算曲率半径
    double min_curvature_radius = DBL_MAX;
    double sum_curvature_radius = 0.0;
    int count_valid = 0;
    
    for (int i = 0; i < num_samples; i++) {
        if (curvature[i] > 1e-6) {
            double radius = 1.0 / curvature[i];
            min_curvature_radius = fmin(min_curvature_radius, radius);
            sum_curvature_radius += radius;
            count_valid++;
        }
    }
    
    double mean_curvature_radius = sum_curvature_radius / (count_valid > 0 ? count_valid : 1);
    
    // 输出路径分析结果
    printf("曲线段路径长度: %.2f mm\n", path_length);
    printf("最大曲率: %.4f\n", max_curv);
    printf("最小曲率半径: %.2f mm\n", min_curvature_radius);
    printf("平均曲率半径: %.2f mm\n", mean_curvature_radius);
    
    // ===== 速度规划 =====
    printf("\n===== 速度规划 =====\n");
    
    // 配置速度规划参数
    SpeedConfig speed_config;
    speed_config.amax_path = 3.0;       // 路径规划中最大允许加速度
    speed_config.amax_speed = 2320.0;   // 第一段直线切向加速度最大值
    speed_config.jerk_factor = 4.5;     // a/J 因子
    speed_config.percent = 21.0;        // 划分二、三段的百分比
    speed_config.time_interval = 0.016; // 16ms采样间隔
    
    // 计算直线段长度
    double dis1 = 0.0; // 第一段直线长度 (这里为0因为全局起点=起点)
    double dis3 = sqrt(pow(global_end_pos.x - end_pos.x, 2) + pow(global_end_pos.y - end_pos.y, 2));
    
    // 进行速度规划
    SpeedProfile *profile = speed_plan_from_curvature(best_spline, dis1, path_length, dis3, speed_config);
    
    // 获取路径采样点
    int num_points;
    PointND **path_points = sample_path_points(best_spline, profile, &num_points);
    
    // 输出速度规划信息
    double max_velocity = 0.0;
    double max_acceleration = 0.0;
    double max_jerk = 0.0;
    
    for (int i = 0; i < profile->num_points; i++) {
        if (profile->velocity[i] > max_velocity) max_velocity = profile->velocity[i];
        if (fabs(profile->acceleration[i]) > max_acceleration) max_acceleration = fabs(profile->acceleration[i]);
        if (fabs(profile->jerk[i]) > max_jerk) max_jerk = fabs(profile->jerk[i]);
    }
    
    double total_time = profile->time[profile->num_points-1];
    
    printf("最大速度: %.2f mm/s\n", max_velocity);
    printf("最大加速度: %.2f mm/s^2\n", max_acceleration);
    printf("最大加加速度: %.2f mm/s^3\n", max_jerk);
    printf("总运动时间: %.3f s\n", total_time);
    printf("曲线段运动时间: %.3f s\n", profile->time[profile->num_points/3*2] - profile->time[profile->num_points/3]);
    
    // 输出一些采样路径点
    printf("\n=== 采样路径点 ===\n");
    int samples_to_print = num_points > 10 ? 10 : num_points;
    for (int i = 0; i < samples_to_print; i++) {
        int idx = i * (num_points / samples_to_print);
        if (idx < num_points && path_points[idx]->dim >= 2) {
            printf("点 %d: (%.3f, %.3f)\n", idx, path_points[idx]->coords[0], path_points[idx]->coords[1]);
        }
    }
    
    // 清理资源
    free(t_samples);
    free(curvature);
    speed_profile_free(profile);
    for (int i = 0; i < num_points; i++) {
        pointnd_free(path_points[i]);
    }
    free(path_points);
    pso
