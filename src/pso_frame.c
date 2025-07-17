#include "pso_frame.h"
#include "basic_support.h"

// 初始化PSO配置
PSOConfig create_pso_config(int dim, int pop_size, int max_iter,
    double w, double c1, double c2, double* init_min, double* init_max,
    double* search_min, double* search_max, double velocity_limit) {
    PSOConfig config = { 0 };
    config.dim = dim;
    config.pop_size = pop_size;
    config.max_iter = max_iter;
    config.w = w;
    config.c1 = c1;
    config.c2 = c2;
    config.init_min = init_min;
    config.init_max = init_max;
    config.search_min = search_min;
    config.search_max = search_max;
    config.velocity_limit = velocity_limit;
    return config;
}

// 初始化粒子
void initialize_particle(Particle* particle, PSOConfig* config) {
    particle->position = (double*)malloc(config->dim * sizeof(double));
    particle->velocity = (double*)malloc(config->dim * sizeof(double));
    particle->best_position = (double*)malloc(config->dim * sizeof(double));

    for (int j = 0; j < config->dim; j++) {
        // 在指定范围内初始化位置
        particle->position[j] = random_in_range(config->init_min[j], config->init_max[j]);

        // 初始化速度
        particle->velocity[j] = random_in_range(-config->velocity_limit, config->velocity_limit);

        particle->best_position[j] = particle->position[j];
    }

    particle->best_fitness = DBL_MAX;
}

// 初始化PSO
PSOState initialize_pso(PSOConfig* config, FitnessFunc fitness_function) {
    srand(time(NULL));

    PSOState state;
    state.particles = (Particle*)malloc(config->pop_size * sizeof(Particle));
    state.global_best_position = (double*)malloc(config->dim * sizeof(double));
    state.global_best_fitness = DBL_MAX;

    // 初始化所有粒子
    for (int i = 0; i < config->pop_size; i++) {
        initialize_particle(&state.particles[i], config);

        // 计算初始适应度
        double fitness = fitness_function(state.particles[i].position, config->dim);
        state.particles[i].best_fitness = fitness;

        // 更新全局最优
        if (fitness < state.global_best_fitness) {
            state.global_best_fitness = fitness;
            for (int j = 0; j < config->dim; j++) {
                state.global_best_position[j] = state.particles[i].position[j];
            }
        }
    }

    return state;
}

// 更新粒子
void update_particle(Particle* particle, PSOConfig* config, double* global_best_position,
    FitnessFunc fitness_function) {
    for (int j = 0; j < config->dim; j++) {
        // 生成[0,1]随机数
        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;

        // 速度更新公式
        particle->velocity[j] = config->w * particle->velocity[j]
            + config->c1 * r1 * (particle->best_position[j] - particle->position[j])
            + config->c2 * r2 * (global_best_position[j] - particle->position[j]);
    
        if (particle->velocity[j] + config->velocity_limit < 0)
            particle->velocity[j] = - config->velocity_limit;
        else if (particle->velocity[j] > config->velocity_limit)
            particle->velocity[j] = config->velocity_limit;

        // 位置更新
        particle->position[j] += particle->velocity[j];

        // 边界处理
        if (particle->position[j] < config->search_min[j]) {
            particle->position[j] = config->search_min[j];
            // particle->velocity[j] *= -0.5;  // 反弹
        } else if (particle->position[j] > config->search_max[j]) {
            particle->position[j] = config->search_max[j];
            // particle->velocity[j] *= -0.5;  // 反弹
        }
    }

    // 计算新适应度
    double current_fitness = fitness_function(particle->position, config->dim);

    // 更新个体最优
    if (current_fitness < particle->best_fitness) {
        particle->best_fitness = current_fitness;
        for (int j = 0; j < config->dim; j++) {
            particle->best_position[j] = particle->position[j];
        }
    }
}

// 执行PSO优化
void run_pso(PSOConfig* config, FitnessFunc fitness_function) {
    PSOState state = initialize_pso(config, fitness_function);

    printf("PSO优化开始\n");
    printf("参数: 维度=%d, 粒子数=%d, 迭代=%d, w=%.2f, c1=%.2f, c2=%.2f\n",
        config->dim, config->pop_size, config->max_iter, config->w, config->c1, config->c2);

    // 主迭代循环
    for (int iter = 0; iter < config->max_iter; iter++) {
        for (int i = 0; i < config->pop_size; i++) {
            update_particle(&state.particles[i], config, state.global_best_position, fitness_function);

            // 更新全局最优
            if (state.particles[i].best_fitness < state.global_best_fitness) {
                state.global_best_fitness = state.particles[i].best_fitness;
                for (int j = 0; j < config->dim; j++) {
                    state.global_best_position[j] = state.particles[i].best_position[j];
                }
            }
        }

        // 每100次迭代输出进度
        if (iter % 100 == 0) {
            printf("迭代 %4d: 当前最优值 = %.6f\n", iter, state.global_best_fitness);
        }
    }

    // 输出最终结果
    printf("\n===== 优化完成 =====\n");
    printf("全局最优值: %.10f\n", state.global_best_fitness);
    printf("最优解位置: (");
    for (int j = 0; j < config->dim; j++) {
        printf("%.6f", state.global_best_position[j]);
        if (j < config->dim - 1) printf(", ");
    }
    printf(")\n");

    // 释放内存
    for (int i = 0; i < config->pop_size; i++) {
        free(state.particles[i].position);
        free(state.particles[i].velocity);
        free(state.particles[i].best_position);
    }
    free(state.particles);
    free(state.global_best_position);
}
