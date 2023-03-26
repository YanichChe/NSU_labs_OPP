#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 6000
#define EPSILON 1e-6
#define MAX_ITERATION_COUNT 50000
#define TAU 1e-5
double calcSquareNorm(const double* vector)
{
    double norm = 0;
    for (int i = 0; i < N; i++)
    {
        norm += vector[i] * vector[i];
    }
    return norm;
}
void generate_vector(double* vector)
{
    for (int i = 0; i < N; i++)
    {
        vector[i] = (double)rand() / RAND_MAX * 10.0 - 5.0;
    }
}

void generate_matrix(double* matrix)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < i; j++)
        {
            matrix[i * N + j] = matrix[j * N + i];
        }

        for(int j = i; j < N; j++)
        {
            matrix[i * N + j] = (double)rand() / RAND_MAX * 2.0 - 1.0; //float in range -1 to 1
            if(i == j) matrix[i * N + j] = matrix[i * N + j] + N;

        }
    }
}

void check(double* A, double*x , double* b){
    double * result = malloc(sizeof(double) * N);
    for(int i = 0; i < N; i++) {
        result[i] = 0;
        for(int j = 0; j < N; j++) {
            result[i] += A[i * N + j] * x[j];
        }
    }
    for (int i = 0; i < N; i++){
        if (result[i] - b[i] < EPSILON || b[i] - result[i] < EPSILON){}
        else{
            printf("not okay\n");
            return;
        }
    }
    printf("correct\n");
}

int main(int argc, char *argv[]) {

    double* A = malloc(sizeof(double) * N * N);
    double* x = malloc(sizeof(double) * N);
    double* b = malloc(sizeof(double) * N);

    double *buffer = malloc(sizeof(double) * N);

    generate_matrix(A);
    generate_vector(x);
    generate_vector(b);

    double norm_b = sqrt(calcSquareNorm(b));
    double res = 1;
    int iterationCount = 1;

    int start, end;
    start = omp_get_wtime();
    unsigned int i, j;
    double norm;

#pragma omg parallel private(i, j)
    {
        while (res > EPSILON && iterationCount < MAX_ITERATION_COUNT) {

#pragma omp parallel for private (j) schedule(runtime)
            for(i = 0; i < N; i++) {
                buffer[i] = 0;
                for(j = 0; j < N; j++) {
                    buffer[i] += A[i * N + j] * x[j];
                }
            }

#pragma omp parallel for schedule(runtime)
            for (i = 0; i < N; i++) {
                buffer[i] = buffer[i] - b[i];
            }

#pragma omp single
            norm = 0;

#pragma omp parallel for reduction (+:norm) schedule(runtime)
            for (i = 0; i < N; i++){
                norm += buffer[i] * buffer[i];
            }

#pragma omp parallel for schedule(runtime)
            for (i = 0; i < N; ++i) {
                x[i] = x[i] - buffer[i] * TAU;
            }

#pragma omp single
            {
                res = sqrt(norm / norm_b);
                printf("%f\n", res);
                iterationCount++;
            }
        }
    }
    end = omp_get_wtime();
    printf("Total time is %d seconds\n", (end - start));
    check(A, x, b);

    free(A);
    free(x);
    free(b);
    free(buffer);

    return 0;
}



