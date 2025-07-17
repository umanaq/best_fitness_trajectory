#ifndef _FITNESS_H_
#define _FITNESS_H_

#include "spline.h"


typedef struct {
    double x;
    double y;
    double radius;
} Obstacle;


// 碰撞检测
int check_collision(Spline* sp, Obstacle* obstacles,
    int num_obstacles, double d, int samples);

// 计算两点间距离
double distance(Point a, Point b);

// 计算曲率
double curvature(Spline* sp, double t);

// 计算路径长度
double path_length(Spline* sp, int samples);

// 端点区域曲率和趋势惩罚
void end_curvature_penalty(Spline* sp, double t_start,
    double* curv_end, double* trend_penalty);


// 主适应度计算函数
double calculate_fitness(
    Point* points,        // 路径点数组
    int num_points,             // 路径点数量
    Obstacle* obstacles,  // 障碍物数组
    int num_obstacles,          // 障碍物数量
    double d,                   // 安全距离
    double desired_dir,         // 终点切线方向
    double factor_end,          // 终点延长系数
    double Kmax,                // 最大允许曲率
    double t_endzone_percent,   // 终点区域百分比
    double alpha,               // 曲率基础惩罚系数
    double gamma,               // 曲率约束惩罚系数
    double tan2,                // 端点曲率惩罚系数
    double trend_factor,        // 曲率趋势惩罚系数
    double path_weight,         // 路径长度权重
    double beta                 // 碰撞惩罚系数
);

#endif  // _FITNESS_H_
