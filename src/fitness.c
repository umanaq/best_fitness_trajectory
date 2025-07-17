#include "fitness.h"


// 碰撞检测
int check_collision(Spline* sp, Obstacle* obstacles,
    int num_obstacles, double d, int samples) {
    for (int i = 0; i <= samples; i++) {
        double t = (1.0 * i) / (1.0 * samples);
        Point pt = spline_ppval(sp, t);

        for (int j = 0; j < num_obstacles; j++) {
            Obstacle obs = obstacles[j];
            double dist = sqrt(pow(pt.x - obs.x, 2) + pow(pt.y - obs.y, 2));
            if (dist < (obs.radius + d)) {
                return 1;
            }
        }
    }
    return 0;
}

// 计算两点间距离
double distance(Point a, Point b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

// 计算曲率
double curvature(Spline* sp, double t) {
    Point d1 = spline_deriv_xy(sp, t, 1);
    Point d2 = spline_deriv_xy(sp, t, 2);

    double dx = d1.x, dy = d1.y;
    double ddx = d2.x, ddy = d2.y;

    double numerator = fabs(dx * ddy - dy * ddx);
    double denominator = pow(dx * dx + dy * dy, 1.5);

    return denominator > 1e-10 ? numerator / denominator : 0.0;
}

// 计算路径长度
double path_length(Spline* sp, int samples) {
    double length = 0.0;
    Point prev = spline_ppval(sp, 0.0);

    for (int i = 1; i <= samples; i++) {
        double t = (1.0 * i) / (1.0 * samples);
        Point curr = spline_ppval(sp, t);
        length += distance(prev, curr);
        prev = curr;
    }
    return length;
}

// 端点区域曲率和趋势惩罚
void end_curvature_penalty(Spline* sp, double t_start,
    double* curv_end, double* trend_penalty) {
    *curv_end = 0.0;
    *trend_penalty = 0.0;

    const int samples = 100;
    double prev_curv = curvature(sp, t_start);

    for (int i = 1; i <= samples; i++) {
        double t = t_start + (1.0 - t_start) * i / samples;
        double curv = curvature(sp, t);

        // 端点曲率（最后10%）
        if (i > samples * 0.9) {
            *curv_end += curv;
        }

        // 曲率趋势惩罚（曲率增加时惩罚）
        if (curv > prev_curv) {
            *trend_penalty += curv - prev_curv;
        }
        prev_curv = curv;
    }
}

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
) {
    // 1. 构造三次样条

    double sl = distance(points[0], points[1]);
    double lv = distance(points[1], points[2]);
    double sv = distance(points[0], points[2]);
    double u[] = { 0, sl / (sl + lv), 1 };  // 参数节点

    double pos_seq[] = {
        points[0].x, points[0].y,
        points[1].x, points[1].y,
        points[2].x, points[2].y
    };
    int dim = 2;  // 向量维度
    // 边界条件: 
    // [2,1] - 左端二阶导数, 右端一阶导数
    int bctype[] = { 2, 1 };
    // 边界值: 
    // 左端二阶导数(x=0, y=0), 右端一阶导数(x=1, y=0)
    double bcval[] = {
        0, cos(desired_dir) * factor_end * sv,  // x(t)左边界二阶导数为0, 二阶导数为
        // 1, 0  // 右边界 (x方向斜率=1, y方向斜率=0)
        // cos(desired_dir), sin(desired_dir)  // 右边界 (x方向斜率=1, y方向斜率=0)
        0, sin(desired_dir) * factor_end * sv  // 右边界 (x方向斜率=1, y方向斜率=0)
    };

    Spline* sp = csape_c(num_points, u, pos_seq, dim, bctype, bcval);
    
    // 2. 计算路径长度
    int samples = 1000;
    double length = path_length(sp, samples);

    // 3. 计算曲率相关指标
    double max_curv = 0.0;
    for (int i = 0; i <= samples; i++) {
        double t = (1.0 * i) / (1.0 * samples);
        double curv = curvature(sp, t);
        if (curv > max_curv) max_curv = curv;
    }

    // 4. 计算端点区域惩罚
    double end_curv, trend_penalty;
    end_curvature_penalty(sp, t_endzone_percent, &end_curv, &trend_penalty);

    // 5. 曲率惩罚项
    double curvature_penalty = alpha * max_curv;
    if (max_curv > Kmax) {
        curvature_penalty += gamma * pow(max_curv - Kmax, 2);
    }
    curvature_penalty += tan2 * end_curv + trend_factor * trend_penalty;

    // 6. 碰撞检测
    double collision_penalty = 0.0;
    if (obstacles && num_obstacles > 0) {
        if (check_collision(sp, obstacles, num_obstacles, d, samples)) {
            collision_penalty = beta;
        }
    }

    // 7. 总适应度
    return path_weight * length + curvature_penalty + collision_penalty;
}
