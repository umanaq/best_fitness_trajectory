#include "pso_frame.h"


// ��ʼ��PSO����
PSOConfig create_pso_config(int dim, int pop_size, int max_iter,
    double w, double c1, double c2,
    double* lower_bounds, double* upper_bounds, double* velocity_range) {
    PSOConfig config = {};
    config.dim = dim;
    config.pop_size = pop_size;
    config.max_iter = max_iter;
    config.w = w;
    config.c1 = c1;
    config.c2 = c2;
    config.lower_bounds = lower_bounds;
    config.upper_bounds = upper_bounds;
    config.velocity_range = velocity_range;
    return config;
}

// ��ʼ������
void initialize_particle(Particle* particle, PSOConfig* config) {
    particle->position = (double*)malloc(config->dim * sizeof(double));
    particle->velocity = (double*)malloc(config->dim * sizeof(double));
    particle->best_position = (double*)malloc(config->dim * sizeof(double));

    for (int j = 0; j < config->dim; j++) {
        // ��ָ����Χ�ڳ�ʼ��λ��
        double range = config->upper_bounds[j] - config->lower_bounds[j];
        particle->position[j] = config->lower_bounds[j] + ((double)rand() / RAND_MAX) * range;

        // ��ʼ���ٶȣ��ٶȷ�ΧΪλ�÷�Χ��10%��
        particle->velocity[j] = range * 0.1 * ((double)rand() / RAND_MAX - 0.5);

        particle->best_position[j] = particle->position[j];
    }

    particle->best_fitness = DBL_MAX;
}

// ��ʼ��PSO
PSOState initialize_pso(PSOConfig* config, FitnessFunc fitness_function) {
    srand(time(NULL));

    PSOState state;
    state.particles = (Particle*)malloc(config->pop_size * sizeof(Particle));
    state.global_best_position = (double*)malloc(config->dim * sizeof(double));
    state.global_best_fitness = DBL_MAX;

    // ��ʼ����������
    for (int i = 0; i < config->pop_size; i++) {
        initialize_particle(&state.particles[i], config);

        // �����ʼ��Ӧ��
        double fitness = fitness_function(state.particles[i].position, config->dim);
        state.particles[i].best_fitness = fitness;

        // ����ȫ������
        if (fitness < state.global_best_fitness) {
            state.global_best_fitness = fitness;
            for (int j = 0; j < config->dim; j++) {
                state.global_best_position[j] = state.particles[i].position[j];
            }
        }
    }

    return state;
}

// ��������
void update_particle(Particle* particle, PSOConfig* config, double* global_best_position, FitnessFunc fitness_function) {
    for (int j = 0; j < config->dim; j++) {
        // ����[0,1]�����
        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;

        // �ٶȸ��¹�ʽ
        particle->velocity[j] = config->w * particle->velocity[j]
            + config->c1 * r1 * (particle->best_position[j] - particle->position[j])
            + config->c2 * r2 * (global_best_position[j] - particle->position[j]);

        // λ�ø���
        particle->position[j] += particle->velocity[j];

        // �߽紦��
        if (particle->position[j] < config->lower_bounds[j])
            particle->position[j] = config->lower_bounds[j];
        if (particle->position[j] > config->upper_bounds[j])
            particle->position[j] = config->upper_bounds[j];
    }

    // ��������Ӧ��
    double current_fitness = fitness_function(particle->position, config->dim);

    // ���¸�������
    if (current_fitness < particle->best_fitness) {
        particle->best_fitness = current_fitness;
        for (int j = 0; j < config->dim; j++) {
            particle->best_position[j] = particle->position[j];
        }
    }
}

// ִ��PSO�Ż�
void run_pso(PSOConfig* config, FitnessFunc fitness_function) {
    PSOState state = initialize_pso(config, fitness_function);

    printf("PSO�Ż���ʼ\n");
    printf("����: ά��=%d, ������=%d, ����=%d, w=%.2f, c1=%.2f, c2=%.2f\n",
        config->dim, config->pop_size, config->max_iter, config->w, config->c1, config->c2);

    // ������ѭ��
    for (int iter = 0; iter < config->max_iter; iter++) {
        for (int i = 0; i < config->pop_size; i++) {
            update_particle(&state.particles[i], config, state.global_best_position, fitness_function);

            // ����ȫ������
            if (state.particles[i].best_fitness < state.global_best_fitness) {
                state.global_best_fitness = state.particles[i].best_fitness;
                for (int j = 0; j < config->dim; j++) {
                    state.global_best_position[j] = state.particles[i].best_position[j];
                }
            }
        }

        // ÿ100�ε����������
        if (iter % 100 == 0) {
            printf("���� %4d: ��ǰ����ֵ = %.6f\n", iter, state.global_best_fitness);
        }
    }

    // ������ս��
    printf("\n===== �Ż���� =====\n");
    printf("ȫ������ֵ: %.10f\n", state.global_best_fitness);
    printf("���Ž�λ��: (");
    for (int j = 0; j < config->dim; j++) {
        printf("%.6f", state.global_best_position[j]);
        if (j < config->dim - 1) printf(", ");
    }
    printf(")\n");

    // �ͷ��ڴ�
    for (int i = 0; i < config->pop_size; i++) {
        free(state.particles[i].position);
        free(state.particles[i].velocity);
        free(state.particles[i].best_position);
    }
    free(state.particles);
    free(state.global_best_position);
}
