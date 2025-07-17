#include "basic_support.h"
#include <stdlib.h>

Matrix matrix_create(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double*)malloc(rows * cols * sizeof(double));
    return mat;
}

void matrix_free(Matrix* mat) {
    if (mat && mat->data) {
        free(mat->data);
        mat->data = NULL;
    }
}

// 随机数生成器: [min, max]范围内的随机浮点数
double random_in_range(double min, double max) {
    return min + ((double)rand() / RAND_MAX) * (max - min);
}
