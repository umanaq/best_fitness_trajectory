#ifndef SPEED_PROFILE_H
#define SPEED_PROFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spline.h"

// 速度规划参数
typedef struct {
    double amax_path;         // 路径规划中最大允许加速度
    double amax_speed;        // 第一段直线切向加速度最大值
    double jerk_factor;       // a/J 因子
    double percent;           // 划分二、三段的百分比
    double time_interval;     // 采样时间间隔(s)
} SpeedConfig;

// 速度规划结果
typedef struct {
    double *time;             // 时间点数组
    double *position;         // 位置数组
    double *velocity;         // 速度数组
    double *acceleration;     // 加速度数组
    double *jerk;             // 加加速度数组
    double *curvature;        // 曲率数组
    double *normal_acc;       // 法向加速度数组
    double *total_acc;        // 合成加速度数组
    int num_points;           // 采样点数量
} SpeedProfile;

// 三段式运动阶段
typedef struct {
    double time;              // 阶段时间
    double displacement;      // 阶段位移
    double velocity;          // 阶段末速度
} MotionStage;

// 基于曲率的速度规划
SpeedProfile* speed_plan_from_curvature(SplinePP *path_spline, 
                                       double s1, double s2, double s3, 
                                       SpeedConfig config);

// 时间点采样获取路径点
PointND** sample_path_points(SplinePP *path_spline, SpeedProfile *profile, int *num_points);

// 计算平滑加减速方案
void calculate_smooth_profile(double v0, double a0, double jerk, 
                             MotionStage *stage1, MotionStage *stage2,
                             MotionStage *stage3, MotionStage *stage4,
                             double *total_time, double *total_displacement);

// 生成时间点采样的轨迹
SpeedProfile* generate_speed_profile(double v0, double a0, double jerk, int num_samples);

// 生成等间隔时间的位移数组
double* choose_time_samples(double *time_array, double *disp_array, int num_points, 
                          double interval, int *output_size);

// 生成直线段坐标点
PointND** generate_line_points(PointND *start, PointND *end, double *dists, int num_points);

// 清理速度规划资源
void speed_profile_free(SpeedProfile *profile);

#endif /* SPEED_PROFILE_H */
