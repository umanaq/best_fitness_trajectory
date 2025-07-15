#include <memory.h>
#include <stdio.h>

#include "pso_frame.h"
//#include "cubic_spline_with_curvature_constraint.h"
#include "basic_support.h"
//
//#include "spline.h"

// ʾ����Ӧ�Ⱥ���: Rastrigin����
double rastrigin(double* x, int dim) {
    double sum = 0.0;
    for (int i = 0; i < dim; i++) {
        sum += x[i] * x[i] - 10 * cos(2 * PI * x[i]);
    }
    return 10 * dim + sum;
}

// ʾ����Ӧ�Ⱥ���: Sphere����
double sphere(double* x, int dim) {
    double sum = 0.0;
    for (int i = 0; i < dim; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

// BestFitnessSpline
double fitness_spline() {
    double fitness = 0.0;

    return fitness;
}

#if 1
int main() {
    // ���������ռ�
    int dim = 2;
    double lower_bounds[] = { -5.12, -5.12 }; // Rastrigin������������Χ
    double upper_bounds[] = { 5.12, 5.12 };
    double velocity_range[] = { -100.0, 100.0 };

    // ����PSO����
    PSOConfig config = create_pso_config(
        dim,            // ����ά��
        30,             // ��������
        1000,           // ����������
        0.5,            // ����Ȩ�� w
        2.0,            // ����ѧϰ���� c1
        2.0,            // Ⱥ��ѧϰ���� c2
        lower_bounds,   // �����½�
        upper_bounds,   // �����Ͻ�
        velocity_range
    );

    // ����PSO�Ż���
    //run_pso(&config, rastrigin); // ����rastrigin����ָ��
    run_pso(&config, sphere); // ���滻Ϊsphere����

    return 0;
}
#endif

#if 0
// ================= ʾ���÷� =================
int main() { // cubic_spline_with....h
    // ����·���� (x,y)
    Matrix points = matrix_create(2, 3);
    double point_data[6] = { -300.588989, -21.019232,
                           -350.0, -30.0,
                           -453.451843, -49.921623 };
    //points.data = point_data;
    memcpy(points.data, point_data, sizeof(*point_data) * 6);

    // �����յ����߷���
    double end_tangent[2] = { -1.0, 0.0 };

    // ���ɴ��߽���������������
    Spline sp = clamped_clamped_spline(&points, end_tangent, 0.35);

    // �ڲ����ռ���Ȳ���
    const int num_samples = 10;
    for (int i = 0; i <= num_samples; i++) {
        double t = (double)i / num_samples;

        // ����λ��
        double x = spline_eval(&sp, t, 0);
        double y = spline_eval(&sp, t, 0);

        // ��������
        double curvature[2];
        spline_curvature(&sp, t, curvature);

        printf("t=%.2f: (%.3f, %.3f) | Curvature: %.6f, Change: %.6f\n",
            t, x, y, curvature[0], curvature[1]);
    }

    // ������Դ
    free(sp.breaks);
    free(sp.coefs);
    matrix_free(&points);

    return 0;
}
#endif

#if 0
// ������ - ���Դ���
int main() { // spline.h
    // ��������
    int n = 3; // �ڵ���
    double u[] = { 0, 0.5, 1 }; // �����ڵ�
    int dim = 2; // ����ά��

    // ����㣺2ά�� (0,0), (0.5,0.5), (1,1)
    double points[] = {
        0, 0,    // ��1 (x=0, y=0)
        0.4, 0.8, // ��2 (x=0.5, y=0.5)
        1, 1      // ��3 (x=1, y=1)
    };

    // �߽�����: 
    // [2,1] - ��˶��׵���, �Ҷ�һ�׵���
    int bctype[] = { 2, 1 };
    // �߽�ֵ: 
    // ��˶��׵���(x=0, y=0), �Ҷ�һ�׵���(x=1, y=0)
    double bcval[] = {
        0, 0, // ��߽� (dim0=x, dim1=y)
        1, 0  // �ұ߽� (x����б��=1, y����б��=0)
    };

    // ��������
    Spline3D* spline = csape_c(n, u, points, dim, bctype, bcval);

    if (spline == NULL) {
        return EXIT_FAILURE;
    }

    // ��[0,1]�����ڲ�������
    int num_samples = 100;
    double result[2]; // �洢��ά���

    printf("����t   x����     y����\n");
    printf("------------------------\n");
    for (int i = 0; i <= num_samples; i++) {
        double t = i / (double)num_samples;
        spline_eval(spline, t, result);
        printf("%.2f   %.6f  %.6f\n", t, result[0], result[1]);
    }

    // �ͷ��ڴ�
    free_spline(spline);

    return EXIT_SUCCESS;
}
#endif
