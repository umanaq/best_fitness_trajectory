#pragma once

#define PI (3.14159265358979323846)

// ================= 数据结构定义 =================
typedef struct {
    double* data;
    int rows;
    int cols;
} Matrix;

// ================= 辅助函数 =================
Matrix matrix_create(int rows, int cols);

void matrix_free(Matrix* mat);