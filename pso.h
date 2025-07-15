#ifndef PSO_H
#define PSO_H

#include "spline.h"

typedef struct {
    double x;
    double y;
    double radius;
} Obstacle;

typedef struct {
    double x_min;        // 搜索空间x最小值
    double x_max;        // 搜索空间x最大值
    double y_min;        // 搜索空间y最小值
    double y_max;        // 搜索空间y最大值
    int num_particles;   // 粒子数量
    int max_iter;        // 最大迭代次数
    int num_waypoints;   // 路径点数量
    double w;            // 惯性权重
    double c1;           // 认知因子
    double c2;           // 社会因子
    double max_velocity; // 最大粒子速度
    double gamma;        // 曲率惩罚系数
    double beta;         // 碰撞惩罚系数
    double alpha;        // 曲率约束系数
    double tan2;         // 终点曲率惩罚系数
    double trend_factor; // 曲率趋势惩罚系数
    double path_weight;  // 路径长度权重
    double factor_end;   // 终点切线系数
    double t_endzone_percent; // 终点区域百分比
    double Kmax;         // 最大允许曲率
    double d;            // 安全距离
} PSOConfig;

typedef struct {
    int num_obstacles;
    Obstacle *obstacles;
} Environment;

typedef struct {
    double **positions;        // 粒子位置数组 [num_particles][dimension]
    double **velocities;       // 粒子速度数组 [num_particles][dimension]
    double **pbest_positions;  // 个体最佳位置 [num_particles][dimension]
    double *pbest_fitness;     // 个体最佳适应度 [num_particles]
    double *gbest_position;    // 全局最佳位置 [dimension]
    double gbest_fitness;      // 全局最佳适应度
    int dimension;             // 维度 = 2 * num_waypoints
    PSOConfig config;          // PSO配置参数
    SplinePP *best_spline;     // 最佳样条曲线
    double **best_path;        // 最佳路径点集
    int best_path_size;        // 最佳路径点数量
} PSO;

// 初始化PSO算法
PSO* pso_init(Point2D start_pos, Point2D end_pos, Point2D desired_end_tangent, 
              Environment *env, PSOConfig config);

// 运行PSO算法
void pso_run(PSO *pso, Point2D start_pos, Point2D end_pos, 
             Point2D desired_end_tangent, Environment *env);

// 计算路径适应度
double pso_calculate_fitness(double *particle_pos, int dimension, 
                             Point2D start_pos, Point2D end_pos, 
                             Point2D desired_end_tangent, Environment *env, 
                             PSOConfig config, SplinePP **spline, 
                             double ***path_ptr, int *path_size);

// 获取最佳路径
SplinePP* pso_get_best_spline(PSO *pso);

// 获取最佳路径点集
double** pso_get_best_path(PSO *pso, int *size);

// 释放PSO资源
void pso_free(PSO *pso);

#endif /* PSO_H */
