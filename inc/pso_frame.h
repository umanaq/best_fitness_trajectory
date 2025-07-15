#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>

// 函数指针类型定义 - 适应度函数
typedef double (*FitnessFunc)(double*, int);

// PSO配置结构体
typedef struct {
    int dim;           // 问题维度
    int pop_size;      // 粒子数量
    int max_iter;      // 最大迭代次数
    double w;          // 惯性权重
    double c1;         // 个体学习因子
    double c2;         // 群体学习因子
    double* lower_bounds; // 搜索空间下界
    double* upper_bounds; // 搜索空间上界
    double* velocity_range; // 搜索速度范围
} PSOConfig;

// 粒子结构体
typedef struct {
    double* position;      // 当前位置
    double* velocity;      // 当前速度
    double* best_position; // 个体历史最优位置
    double best_fitness;   // 个体历史最优适应度
} Particle;

// PSO状态结构体
typedef struct {
    Particle* particles;       // 粒子群
    double* global_best_position; // 全局最优位置
    double global_best_fitness;   // 全局最优适应度
} PSOState;

// 初始化PSO配置
PSOConfig create_pso_config(int dim, int pop_size, int max_iter,
    double w, double c1, double c2,
    double* lower_bounds, double* upper_bounds, double* velocity_range);

// 初始化粒子
void initialize_particle(Particle* particle, PSOConfig* config);

// 初始化PSO
PSOState initialize_pso(PSOConfig* config, FitnessFunc fitness_function);

// 更新粒子
void update_particle(Particle* particle, PSOConfig* config, double* global_best_position, FitnessFunc fitness_function);

// 执行PSO优化
void run_pso(PSOConfig* config, FitnessFunc fitness_function);
