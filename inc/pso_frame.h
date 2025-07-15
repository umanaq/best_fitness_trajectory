#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>

// ����ָ�����Ͷ��� - ��Ӧ�Ⱥ���
typedef double (*FitnessFunc)(double*, int);

// PSO���ýṹ��
typedef struct {
    int dim;           // ����ά��
    int pop_size;      // ��������
    int max_iter;      // ����������
    double w;          // ����Ȩ��
    double c1;         // ����ѧϰ����
    double c2;         // Ⱥ��ѧϰ����
    double* lower_bounds; // �����ռ��½�
    double* upper_bounds; // �����ռ��Ͻ�
    double* velocity_range; // �����ٶȷ�Χ
} PSOConfig;

// ���ӽṹ��
typedef struct {
    double* position;      // ��ǰλ��
    double* velocity;      // ��ǰ�ٶ�
    double* best_position; // ������ʷ����λ��
    double best_fitness;   // ������ʷ������Ӧ��
} Particle;

// PSO״̬�ṹ��
typedef struct {
    Particle* particles;       // ����Ⱥ
    double* global_best_position; // ȫ������λ��
    double global_best_fitness;   // ȫ��������Ӧ��
} PSOState;

// ��ʼ��PSO����
PSOConfig create_pso_config(int dim, int pop_size, int max_iter,
    double w, double c1, double c2,
    double* lower_bounds, double* upper_bounds, double* velocity_range);

// ��ʼ������
void initialize_particle(Particle* particle, PSOConfig* config);

// ��ʼ��PSO
PSOState initialize_pso(PSOConfig* config, FitnessFunc fitness_function);

// ��������
void update_particle(Particle* particle, PSOConfig* config, double* global_best_position, FitnessFunc fitness_function);

// ִ��PSO�Ż�
void run_pso(PSOConfig* config, FitnessFunc fitness_function);
