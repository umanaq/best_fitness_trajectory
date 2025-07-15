#pragma once

#define PI (3.14159265358979323846)

// ================= ���ݽṹ���� =================
typedef struct {
    double* data;
    int rows;
    int cols;
} Matrix;

// ================= �������� =================
Matrix matrix_create(int rows, int cols);

void matrix_free(Matrix* mat);