#include "speed_profile.h"
#include <string.h>
#include <float.h>

// 基于曲率的速度规划
SpeedProfile* speed_plan_from_curvature(SplinePP *path_spline, 
                                       double s1, double s2, double s3, 
                                       SpeedConfig config) {
    int N = 1000; // 采样点数量
    
    // 创建速度规划结果结构
    SpeedProfile *profile = (SpeedProfile*)malloc(sizeof(SpeedProfile));
    profile->num_points = N * 3; // 三个段的总点数
    profile->time = (double*)malloc(profile->num_points * sizeof(double));
    profile->position = (double*)malloc(profile->num_points * sizeof(double));
    profile->velocity = (double*)malloc(profile->num_points * sizeof(double));
    profile->acceleration = (double*)malloc(profile->num_points * sizeof(double));
    profile->jerk = (double*)malloc(profile->num_points * sizeof(double));
    profile->curvature = (double*)malloc(profile->num_points * sizeof(double));
    profile->normal_acc = (double*)malloc(profile->num_points * sizeof(double));
    profile->total_acc = (double*)malloc(profile->num_points * sizeof(double));
    
    // 计算采样点处的曲率
    double *t_samples = linspace(0, 1, N);
    double *curvature = (double*)malloc(N * sizeof(double));
    double max_curv = 0.0;
    
    spline_compute_curvature(path_spline, t_samples, N, curvature, &max_curv);
    
    // 计算曲率半径
    double *curvature_radius = (double*)malloc(N * sizeof(double));
    double min_radius = DBL_MAX;
    
    for (int i = 0; i < N; i++) {
        if (curvature[i] > 1e-6) {
            curvature_radius[i] = 1.0 / curvature[i];
            if (curvature_radius[i] < min_radius) {
                min_radius = curvature_radius[i];
            }
        } else {
            curvature_radius[i] = DBL_MAX;
        }
    }
    
    // 基于最小曲率半径计算最大速度
    double vmax = sqrt(config.amax_path * min_radius * 0.001) * 0.999;
    
    // 第一段速度规划 (加速)
    double Jmax = config.amax_speed / config.jerk_factor / 0.016;
    
    double t11 = config.amax_speed / Jmax;
    double v1 = 0.5 * Jmax * t11 * t11;
    double a1 = Jmax * t11;
    double d1 = Jmax * t11 * t11 * t11 / 6;
    double d2 = s2 - d1;
    
    // 二分法求解t12
    double t12 = 0.1;
    double t12_min = 0.0;
    double t12_max = 10.0;
    double tolerance = 1e-6;
    double error = 1.0;
    
    while (fabs(error) > tolerance) {
        double F = Jmax * t11 * t12 * t12 / 2 + 0.5 * Jmax * t11 * t11 * t12 - d2;
        error = F;
        
        if (error > 0) {
            t12_max = t12;
        } else {
            t12_min = t12;
        }
        
        t12 = (t12_min + t12_max) / 2;
        
        if (t12_max - t12_min < tolerance) {
            break;
        }
    }
    
    double v2 = Jmax * t11 * t12 + 0.5 * Jmax * t11 * t11;
    double a2 = Jmax * t11;
    
    // 填充第一段速度规划数据
    int n_part1 = N / 2;
    int n_part2 = N - n_part1;
    
    double *t1_part1 = linspace(0, t11, n_part1);
    double *t1_part2 = linspace(t11, t11 + t12, n_part2);
    
    for (int i = 0; i < n_part1; i++) {
        double t = t1_part1[i];
        profile->time[i] = t;
        profile->position[i] = (1.0/6.0) * Jmax * t * t * t;
        profile->velocity[i] = 0.5 * Jmax * t * t;
        profile->acceleration[i] = Jmax * t;
        profile->jerk[i] = Jmax;
    }
    
    for (int i = 0; i < n_part2; i++) {
        double t = t1_part2[i];
        double dt = t - t11;
        int idx = n_part1 + i;
        
        profile->time[idx] = t;
        profile->position[idx] = d1 + v1 * dt + 0.5 * a1 * dt * dt;
        profile->velocity[idx] = v1 + a1 * dt;
        profile->acceleration[idx] = a1;
        profile->jerk[idx] = 0;
    }
    
    // 第二段速度规划 (匀速或部分匀速)
    double main_s2_segment = s3 * config.percent / 100.0;
    
    // 二分法求解t13
    double t13 = t12;
    double t13_min = t12;
    double t13_max = 10.0;
    error = 1.0;
    
    while (fabs(error) > tolerance) {
        double F = Jmax * t11 * t13 * t13 / 2 + 0.5 * Jmax * t11 * t11 * t13 - main_s2_segment - d2;
        error = F;
        
        if (error > 0) {
            t13_max = t13;
        } else {
            t13_min = t13;
        }
        
        t13 = (t13_min + t13_max) / 2;
        
        if (t13_max - t13_min < tolerance) {
            break;
        }
    }
    
    double t3 = t13 - t12;
    double t2_start = t11 + t12;
    double t2_end = t2_start + t3;
    
    // 填充第二段速度规划数据
    double *t2_all = linspace(t2_start, t2_end, N);
    
    for (int i = 0; i < N; i++) {
        double t = t2_all[i];
        double dt = t - t2_start;
        int idx = N + i;
        
        profile->time[idx] = t;
        profile->position[idx] = s2 + v2 * dt + 0.5 * a2 * dt * dt;
        profile->velocity[idx] = v2 + a2 * dt;
        profile->acceleration[idx] = a2;
        profile->jerk[idx] = 0;
    }
    
    // 第三段速度规划 (减速)
    double a0 = config.amax_speed;
    double v0 = profile->velocity[N*2-1]; // 使用第二段末尾的速度
    double J_k = a0 / 0.016 / 5;
    
    // 寻找合适的J_max使得最终位移满足要求
    double J_max = J_k;
    double remaining_s3 = s3 * (100 - config.percent) / 100.0;
    
    // 二分法搜索合适的J_max
    double J_min = J_k / 20;
    double J_max_val = J_k * 20;
    double best_error = DBL_MAX;
    double best_J = J_k;
    
    for (int iter = 0; iter < 50; iter++) {
        J_max = (J_min + J_max_val) / 2;
        
        MotionStage stage1, stage2, stage3, stage4;
        double total_time, total_displacement;
        
        calculate_smooth_profile(v0, a0, J_max, &stage1, &stage2, &stage3, &stage4,
                              &total_time, &total_displacement);
                              
        double current_error = fabs(total_displacement - remaining_s3);
        
        if (current_error < best_error) {
            best_error = current_error;
            best_J = J_max;
        }
        
        if (total_displacement > remaining_s3) {
            J_min = J_max;
        } else {
            J_max_val = J_max;
        }
        
        if (fabs(total_displacement - remaining_s3) < 0.001 || 
            fabs(J_max_val - J_min) < 0.001) {
            break;
        }
    }
    
    // 使用找到的最佳J_max生成第三段轨迹
    SpeedProfile *profile3 = generate_speed_profile(v0, a0, best_J, N);
    
    // 将第三段轨迹加入总体结果
    double t3_offset = profile->time[N*2-1];
    
    for (int i = 0; i < N; i++) {
        int idx = N*2 + i;
        
        profile->time[idx] = t3_offset + profile3->time[i];
        profile->position[idx] = profile->position[N*2-1] + profile3->position[i];
        profile->velocity[idx] = profile3->velocity[i];
        profile->acceleration[idx] = profile3->acceleration[i];
        profile->jerk[idx] = profile3->jerk[i];
    }
    
    // 计算曲率和法向加速度
    for (int i = 0; i < profile->num_points; i++) {
        double t_normalized;
        
        if (i < N) {
            // 第一段
            t_normalized = (double)i / N;
        } else if (i < N*2) {
            // 第二段
            t_normalized = (double)(i-N) / N;
        } else {
            // 第三段
            t_normalized = (double)(i-N*2) / N;
        }
        
        // 找到最近的曲率值
        int curv_idx = (int)(t_normalized * (N-1));
        if (curv_idx >= N) curv_idx = N-1;
        
        profile->curvature[i] = curvature[curv_idx];
        
        double v_squared = profile->velocity[i] * profile->velocity[i];
        double normal_acc = v_squared * profile->curvature[i];
        profile->normal_acc[i] = normal_acc;
        
        // 计算合成加速度
        profile->total_acc[i] = sqrt(normal_acc*normal_acc + 
                                   profile->acceleration[i]*profile->acceleration[i]);
    }
    
    // 清理
    free(t_samples);
    free(curvature);
    free(curvature_radius);
    free(t1_part1);
    free(t1_part2);
    free(t2_all);
    speed_profile_free(profile3);
    
    return profile;
}

// 计算平滑加减速方案
void calculate_smooth_profile(double v0, double a0, double jerk, 
                             MotionStage *stage1, MotionStage *stage2,
                             MotionStage *stage3, MotionStage *stage4,
                             double *total_time, double *total_displacement) {
    // Stage 1: 减小加速度
    double t1 = a0 / jerk;
    double s1 = v0 * t1 + 0.5 * a0 * t1 * t1 - (1.0/6.0) * jerk * t1 * t1 * t1;
    double v1 = v0 + a0 * t1 - 0.5 * jerk * t1 * t1;
    
    // Stage 2: 增加负加速度
    double t2 = a0 / jerk;
    double s2 = s1 + v1 * t2 - (1.0/6.0) * jerk * t2 * t2 * t2;
    double v2 = v1 - 0.5 * jerk * t2 * t2;
    
    // Stage 3: 均匀减速
    double t3 = (v2 - 0.5 * a0 * a0 / jerk) / a0;
    if (t3 < 0) t3 = 0;
    double s3 = s2 + v2 * t3 - 0.5 * a0 * t3 * t3;
    double v3 = v2 - a0 * t3;
    
    // Stage 4: 降低减速度到零
    double t4 = a0 / jerk;
    double s4 = s3 + v3 * t4 - 0.5 * a0 * t4 * t4 + (1.0/6.0) * jerk * t4 * t4 * t4;
    
    // 填充返回值
    stage1->time = t1;
    stage1->displacement = s1;
    stage1->velocity = v1;
    
    stage2->time = t2;
    stage2->displacement = s2;
    stage2->velocity = v2;
    
    stage3->time = t3;
    stage3->displacement = s3;
    stage3->velocity = v3;
    
    stage4->time = t4;
    stage4->displacement = s4;
    stage4->velocity = 0; // 终点速度为0
    
    *total_time = t1 + t2 + t3 + t4;
    *total_displacement = s4;
}

// 生成时间点采样的轨迹
SpeedProfile* generate_speed_profile(double v0, double a0, double jerk, int num_samples) {
    MotionStage stage1, stage2, stage3, stage4;
    double total_time, total_displacement;
    
    calculate_smooth_profile(v0, a0, jerk, &stage1, &stage2, &stage3, &stage4,
                           &total_time, &total_displacement);
    
    SpeedProfile *profile = (SpeedProfile*)malloc(sizeof(SpeedProfile));
    profile->num_points = num_samples;
    profile->time = (double*)malloc(num_samples * sizeof(double));
    profile->position = (double*)malloc(num_samples * sizeof(double));
    profile->velocity = (double*)malloc(num_samples * sizeof(double));
    profile->acceleration = (double*)malloc(num_samples * sizeof(double));
    profile->jerk = (double*)malloc(num_samples * sizeof(double));
    
    // 分配其他数组但不使用
    profile->curvature = (double*)calloc(num_samples, sizeof(double));
    profile->normal_acc = (double*)calloc(num_samples, sizeof(double));
    profile->total_acc = (double*)calloc(num_samples, sizeof(double));
    
    // 采样点
    double *t_all = linspace(0, total_time, num_samples);
    
    for (int i = 0; i < num_samples; i++) {
        double t = t_all[i];
        profile->time[i] = t;
        
        // 确定当前时刻所在的阶段
        if (t < stage1.time) {
            // Stage 1: 减小加速度
            double dt = t;
            profile->position[i] = v0 * dt + 0.5 * a0 * dt * dt - (1.0/6.0) * jerk * dt * dt * dt;
            profile->velocity[i] = v0 + a0 * dt - 0.5 * jerk * dt * dt;
            profile->acceleration[i] = a0 - jerk * dt;
            profile->jerk[i] = -jerk;
        }
        else if (t < (stage1.time + stage2.time)) {
            // Stage 2: 增加负加速度
            double dt = t - stage1.time;
            profile->position[i] = stage1.displacement + stage1.velocity * dt - (1.0/6.0) * jerk * dt * dt * dt;
            profile->velocity[i] = stage1.velocity - 0.5 * jerk * dt * dt;
            profile->acceleration[i] = -jerk * dt;
            profile->jerk[i] = -jerk;
        }
        else if (t < (stage1.time + stage2.time + stage3.time)) {
            // Stage 3: 均匀减速
            double dt = t - (stage1.time + stage2.time);
            profile->position[i] = stage2.displacement + stage2.velocity * dt - 0.5 * a0 * dt * dt;
            profile->velocity[i] = stage2.velocity - a0 * dt;
            profile->acceleration[i] = -a0;
            profile->jerk[i] = 0;
        }
        else {
            // Stage 4: 降低减速度到零
            double dt = t - (stage1.time + stage2.time + stage3.time);
            profile->position[i] = stage3.displacement + stage3.velocity * dt - 0.5 * a0 * dt * dt + (1.0/6.0) * jerk * dt * dt * dt;
            profile->velocity[i] = stage3.velocity - a0 * dt + 0.5 * jerk * dt * dt;
            profile->acceleration[i] = -a0 + jerk * dt;
            profile->jerk[i] = jerk;
        }
    }
    
    free(t_all);
    return profile;
}

// 生成等间隔时间的位移数组
double* choose_time_samples(double *time_array, double *disp_array, int num_points, 
                          double interval, int *output_size) {
    double t_length = time_array[num_points-1] - time_array[0];
    int n1 = (int)ceil(t_length / interval);
    double t_stop = t_length / n1;
    
    // 创建目标时间点
    double *t_targets = (double*)malloc(n1 * sizeof(double));
    for (int i = 0; i < n1; i++) {
        t_targets[i] = time_array[0] + i * t_stop;
    }
    
    // 找到最接近目标时间点的索引
    int *indices = (int*)malloc(n1 * sizeof(int));
    for (int i = 0; i < n1; i++) {
        double min_diff = DBL_MAX;
        indices[i] = 0;
        
        for (int j = 0; j < num_points; j++) {
            double diff = fabs(time_array[j] - t_targets[i]);
            if (diff < min_diff) {
                min_diff = diff;
                indices[i] = j;
            }
        }
    }
    
    // 去除重复索引
    int unique_count = 0;
    int *unique_indices = (int*)malloc(n1 * sizeof(int));
    
    for (int i = 0; i < n1; i++) {
        int is_unique = 1;
        for (int j = 0; j < unique_count; j++) {
            if (indices[i] == unique_indices[j]) {
                is_unique = 0;
                break;
            }
        }
        
        if (is_unique) {
            unique_indices[unique_count++] = indices[i];
        }
    }
    
    // 提取对应的位移值
    double *s_targets = (double*)malloc(unique_count * sizeof(double));
    for (int i = 0; i < unique_count; i++) {
        s_targets[i] = disp_array[unique_indices[i]];
    }
    
    *output_size = unique_count;
    
    // 清理
    free(t_targets);
    free(indices);
    free(unique_indices);
    
    return s_targets;
}

// 生成直线段坐标点
PointND** generate_line_points(PointND *start, PointND *end, double *dists, int num_points) {
    int dim = start->dim;
    
    // 计算方向向量
    double *direction = (double*)malloc(dim * sizeof(double));
    double length = 0.0;
    
    for (int i = 0; i < dim; i++) {
        direction[i] = end->coords[i] - start->coords[i];
        length += direction[i] * direction[i];
    }
    length = sqrt(length);
    
    // 归一化方向向量
    if (length > 1e-10) {
        for (int i = 0; i < dim; i++) {
            direction[i] /= length;
        }
    }
    
    // 生成点
    PointND **points = (PointND**)malloc(num_points * sizeof(PointND*));
    
    for (int i = 0; i < num_points; i++) {
        points[i] = pointnd_create(dim);
        
        for (int j = 0; j < dim; j++) {
            points[i]->coords[j] = start->coords[j] + dists[i] * direction[j];
        }
    }
    
    free(direction);
    return points;
}

// 时间点采样获取路径点
PointND** sample_path_points(SplinePP *path_spline, SpeedProfile *profile, int *num_points) {
    // 生成等间隔时间的位移采样
    int output_size;
    double *s_targets = choose_time_samples(profile->time, profile->position, 
                                          profile->num_points, 0.016, &output_size);
    
    // 采样曲线上的点
    PointND **points = (PointND**)malloc(output_size * sizeof(PointND*));
    
    for (int i = 0; i < output_size; i++) {
        // 找到对应的归一化参数t
        double t = s_targets[i] / profile->position[profile->num_points-1];
        if (t > 1.0) t = 1.0;
        if (t < 0.0) t = 0.0;
        
        // 计算样条点
        double *point_coords = spline_evaluate(path_spline, t);
        
        // 创建点
        points[i] = pointnd_create(path_spline->dim);
        memcpy(points[i]->coords, point_coords, path_spline->dim * sizeof(double));
        
        free(point_coords);
    }
    
    *num_points = output_size;
    free(s_targets);
    return points;
}

// 清理速度规划资源
void speed_profile_free(SpeedProfile *profile) {
    if (profile) {
        if (profile->time) free(profile->time);
        if (profile->position) free(profile->position);
        if (profile->velocity) free(profile->velocity);
        if (profile->acceleration) free(profile->acceleration);
        if (profile->jerk) free(profile->jerk);
        if (profile->curvature) free(profile->curvature);
        if (profile->normal_acc) free(profile->normal_acc);
        if (profile->total_acc) free(profile->total_acc);
        free(profile);
    }
}
