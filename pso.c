#include "pso.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>

// 计算二维向量长度
double norm2d(double x, double y) {
    return sqrt(x*x + y*y);
}

// 生成[0,1]随机数
double rand_01() {
    return (double)rand() / RAND_MAX;
}

// 生成[min,max]随机数
double rand_range(double min, double max) {
    return min + (max - min) * rand_01();
}

// 计算两点间距离
double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

// 创建粒子群优化器
PSO* pso_init(Point2D start_pos, Point2D end_pos, Point2D desired_end_tangent, 
              Environment *env, PSOConfig config) {
    PSO *pso = (PSO*)malloc(sizeof(PSO));
    pso->dimension = 2 * config.num_waypoints;
    pso->config = config;
    pso->best_spline = NULL;
    pso->best_path = NULL;
    pso->best_path_size = 0;
    
    // 分配内存
    pso->positions = (double**)malloc(config.num_particles * sizeof(double*));
    pso->velocities = (double**)malloc(config.num_particles * sizeof(double*));
    pso->pbest_positions = (double**)malloc(config.num_particles * sizeof(double*));
    pso->pbest_fitness = (double*)malloc(config.num_particles * sizeof(double));
    pso->gbest_position = (double*)malloc(pso->dimension * sizeof(double));
    pso->gbest_fitness = DBL_MAX;
    
    // 初始化粒子
    for (int i = 0; i < config.num_particles; i++) {
        pso->positions[i] = (double*)malloc(pso->dimension * sizeof(double));
        pso->velocities[i] = (double*)malloc(pso->dimension * sizeof(double));
        pso->pbest_positions[i] = (double*)malloc(pso->dimension * sizeof(double));
        pso->pbest_fitness[i] = DBL_MAX;
        
        // 随机初始化位置
        for (int j = 0; j < pso->dimension; j += 2) {
            pso->positions[i][j] = rand_range(config.x_min, config.x_max);
            pso->positions[i][j+1] = rand_range(config.y_min, config.y_max);
            
            // 初始化速度
            pso->velocities[i][j] = rand_range(-config.max_velocity, config.max_velocity);
            pso->velocities[i][j+1] = rand_range(-config.max_velocity, config.max_velocity);
            
            // 初始化个体最佳位置
            pso->pbest_positions[i][j] = pso->positions[i][j];
            pso->pbest_positions[i][j+1] = pso->positions[i][j+1];
        }
    }
    
    return pso;
}

// 计算粒子适应度并更新最佳样条和路径
double pso_calculate_fitness(double *particle_pos, int dimension, 
                            Point2D start_pos, Point2D end_pos, 
                            Point2D desired_end_tangent, Environment *env, 
                            PSOConfig config, SplinePP **spline,
                            double ***path_ptr, int *path_size) {
    int num_waypoints = dimension / 2;
    
    // 创建完整路径点（起点+中间点+终点）
    Point2D *waypoints = (Point2D*)malloc((num_waypoints + 2) * sizeof(Point2D));
    waypoints[0] = start_pos;
    
    for (int i = 0; i < num_waypoints; i++) {
        waypoints[i+1].x = particle_pos[2*i];
        waypoints[i+1].y = particle_pos[2*i+1];
    }
    
    waypoints[num_waypoints+1] = end_pos;
    
    // 创建三次样条曲线
    SplinePP *sp = spline_create(waypoints, num_waypoints+2, desired_end_tangent, config.factor_end);
    
    // 采样样条曲线
    int n_samples = 1000;
    double *t_samples = linspace(0, 1, n_samples);
    double **path = (double**)malloc(2 * sizeof(double*));
    path[0] = (double*)malloc(n_samples * sizeof(double));
    path[1] = (double*)malloc(n_samples * sizeof(double));
    
    for (int k = 0; k < n_samples; k++) {
        // 评估样条在t_samples[k]处的值
        double *point = spline_evaluate(sp, t_samples[k]);
        path[0][k] = point[0];
        path[1][k] = point[1];
        free(point);
    }
    
    // 计算路径长度
    double path_length = 0.0;
    for (int k = 1; k < n_samples; k++) {
        double dx = path[0][k] - path[0][k-1];
        double dy = path[1][k] - path[1][k-1];
        path_length += sqrt(dx*dx + dy*dy);
    }
    
    // 计算曲率
    double *curvature = (double*)malloc(n_samples * sizeof(double));
    double max_curv = 0.0;
    spline_compute_curvature(sp, t_samples, n_samples, curvature, &max_curv);
    
    // 计算终点曲率适应度
    double trend_penalty = 0.0;
    double curvature_end = spline_calculate_fitness_end(sp, config.t_endzone_percent, &trend_penalty);
    
    // 计算曲率惩罚
    double curvature_penalty = config.alpha * max_curv;
    if (max_curv > config.Kmax) {
        curvature_penalty += config.gamma * (max_curv - config.Kmax) * (max_curv - config.Kmax);
    }
    curvature_penalty += config.tan2 * curvature_end + config.trend_factor * trend_penalty;
    
    // 碰撞检测
    int collision = 0;
    for (int obs = 0; obs < env->num_obstacles; obs++) {
        double x_center = env->obstacles[obs].x;
        double y_center = env->obstacles[obs].y;
        double radius = env->obstacles[obs].radius;
        
        for (int p = 0; p < n_samples; p++) {
            double dist = distance(path[0][p], path[1][p], x_center, y_center);
            if (dist < (radius + config.d)) {
                collision = 1;
                break;
            }
        }
        if (collision) break;
    }
    double collision_penalty = config.beta * collision;
    
    // 总适应度
    double fitness = config.path_weight * path_length + curvature_penalty + collision_penalty;
    
    // 如果需要返回样条和路径
    if (spline != NULL) {
        *spline = sp;
    } else {
        spline_free(sp);
    }
    
    if (path_ptr != NULL && path_size != NULL) {
        *path_ptr = path;
        *path_size = n_samples;
    } else {
        free(path[0]);
        free(path[1]);
        free(path);
    }
    
    free(waypoints);
    free(t_samples);
    free(curvature);
    
    return fitness;
}

// 运行PSO算法
void pso_run(PSO *pso, Point2D start_pos, Point2D end_pos, 
            Point2D desired_end_tangent, Environment *env) {
    // 初始适应度评估
    for (int i = 0; i < pso->config.num_particles; i++) {
        double fitness = pso_calculate_fitness(
            pso->positions[i], pso->dimension, 
            start_pos, end_pos, desired_end_tangent, env, 
            pso->config, NULL, NULL, NULL
        );
        
        pso->pbest_fitness[i] = fitness;
        
        // 更新全局最佳
        if (fitness < pso->gbest_fitness) {
            pso->gbest_fitness = fitness;
            memcpy(pso->gbest_position, pso->positions[i], pso->dimension * sizeof(double));
            
            // 保存全局最佳样条和路径
            if (pso->best_spline != NULL) {
                spline_free(pso->best_spline);
                pso->best_spline = NULL;
            }
            
            if (pso->best_path != NULL) {
                free(pso->best_path[0]);
                free(pso->best_path[1]);
                free(pso->best_path);
                pso->best_path = NULL;
            }
            
            pso_calculate_fitness(
                pso->gbest_position, pso->dimension, 
                start_pos, end_pos, desired_end_tangent, env, 
                pso->config, &pso->best_spline, &pso->best_path, &pso->best_path_size
            );
        }
    }
    
    // PSO迭代
    for (int iter = 0; iter < pso->config.max_iter; iter++) {
        for (int i = 0; i < pso->config.num_particles; i++) {
            // 更新粒子速度和位置
            for (int j = 0; j < pso->dimension; j++) {
                // 更新速度
                pso->velocities[i][j] = pso->config.w * pso->velocities[i][j] + 
                                       pso->config.c1 * rand_01() * (pso->pbest_positions[i][j] - pso->positions[i][j]) + 
                                       pso->config.c2 * rand_01() * (pso->gbest_position[j] - pso->positions[i][j]);
                
                // 速度限制
                if (pso->velocities[i][j] > pso->config.max_velocity) {
                    pso->velocities[i][j] = pso->config.max_velocity;
                } else if (pso->velocities[i][j] < -pso->config.max_velocity) {
                    pso->velocities[i][j] = -pso->config.max_velocity;
                }
                
                // 更新位置
                pso->positions[i][j] += pso->velocities[i][j];
                
                // 位置边界检查
                if (j % 2 == 0) { // x坐标
                    if (pso->positions[i][j] < pso->config.x_min) {
                        pso->positions[i][j] = pso->config.x_min;
                    } else if (pso->positions[i][j] > pso->config.x_max) {
                        pso->positions[i][j] = pso->config.x_max;
                    }
                } else { // y坐标
                    if (pso->positions[i][j] < pso->config.y_min) {
                        pso->positions[i][j] = pso->config.y_min;
                    } else if (pso->positions[i][j] > pso->config.y_max) {
                        pso->positions[i][j] = pso->config.y_max;
                    }
                }
            }
            
            // 评估适应度
            double fitness = pso_calculate_fitness(
                pso->positions[i], pso->dimension, 
                start_pos, end_pos, desired_end_tangent, env, 
                pso->config, NULL, NULL, NULL
            );
            
            // 更新个体最佳
            if (fitness < pso->pbest_fitness[i]) {
                pso->pbest_fitness[i] = fitness;
                memcpy(pso->pbest_positions[i], pso->positions[i], pso->dimension * sizeof(double));
                
                    // 更新全局最佳
                    if (fitness < pso->gbest_fitness) {
                        pso->gbest_fitness = fitness;
                        memcpy(pso->gbest_position, pso->positions[i], pso->dimension * sizeof(double));
                        
                        // 保存全局最佳样条和路径
                        if (pso->best_spline != NULL) {
                            spline_free(pso->best_spline);
                            pso->best_spline = NULL;
                        }
                        
                        if (pso->best_path != NULL) {
                            free(pso->best_path[0]);
                            free(pso->best_path[1]);
                            free(pso->best_path);
                            pso->best_path = NULL;
                        }
                        
                        pso_calculate_fitness(
                            pso->gbest_position, pso->dimension, 
                            start_pos, end_pos, desired_end_tangent, env, 
                            pso->config, &pso->best_spline, &pso->best_path, &pso->best_path_size
                        );
                    }
                }
            }
        }
        
        printf("Iteration %d, Best Fitness: %.4f\n", iter+1, pso->gbest_fitness);
    }
}

// 获取最佳样条曲线
SplinePP* pso_get_best_spline(PSO *pso) {
    return pso->best_spline;
}

// 获取最佳路径点集
double** pso_get_best_path(PSO *pso, int *size) {
    if (size != NULL) {
        *size = pso->best_path_size;
    }
    return pso->best_path;
}

// 释放PSO资源
void pso_free(PSO *pso) {
    if (pso == NULL) return;
    
    // 释放粒子内存
    for (int i = 0; i < pso->config.num_particles; i++) {
        if (pso->positions[i]) free(pso->positions[i]);
        if (pso->velocities[i]) free(pso->velocities[i]);
        if (pso->pbest_positions[i]) free(pso->pbest_positions[i]);
    }
    
    if (pso->positions) free(pso->positions);
    if (pso->velocities) free(pso->velocities);
    if (pso->pbest_positions) free(pso->pbest_positions);
    if (pso->pbest_fitness) free(pso->pbest_fitness);
    if (pso->gbest_position) free(pso->gbest_position);
    
    // 释放最佳样条和路径
    if (pso->best_spline) spline_free(pso->best_spline);
    
    if (pso->best_path) {
        free(pso->best_path[0]);
        free(pso->best_path[1]);
        free(pso->best_path);
    }
    
    free(pso);
}
