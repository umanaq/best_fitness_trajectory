#ifndef _BASIC_SUPPORT_H_  
#define _BASIC_SUPPORT_H_  

#define PI (3.14159265358979323846)

// ================= 数据结构定义 =================
typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    double* data;
    int rows;
    int cols;
} Matrix;

// ================= 辅助函数 =================
Matrix matrix_create(int rows, int cols);

void matrix_free(Matrix* mat);

double random_in_range(double min, double max);

#endif  // _BASIC_SUPPORT_H_  
