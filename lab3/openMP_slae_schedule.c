#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 5000
#define EPSILON 1e-6
#define MAX_ITERATION_COUNT 50000
#define TAU  1e-5

void generate_matrix(double* matrix);
void generate_vector(double* vector);

void print_matrix(const double* matrix);
void print_vector(const double* vector);

double count_square_norm(const double *vector, int size);
void set_matrix_part(int* line_counts, int* offsets, int size, int thread_num);

void mul(const double* matrix, const double* vector, double* result, int lines);
void sub_vectors(const double* vector1, const double* vector2, double* result, int size);
void count_new_x(double* x, const double* vector, int size);
void check(double* A, double*x , double* b);

int main(int argc, char **argv)
{
    double* A = malloc(sizeof(double) * N * N);
    double* x = malloc(sizeof(double) * N);
    double* b = malloc(sizeof(double) * N);

    int  num_threads = omp_get_max_threads();
    int* line_counts = malloc(sizeof(int) * num_threads);
    int* offsets = malloc(sizeof(int) * num_threads);
    double *buffer = malloc(sizeof(double) * N);

    set_matrix_part(line_counts, offsets, N, num_threads);

    generate_matrix(A);
    generate_vector(x);
    generate_vector(b);

    double b_norm = sqrt(count_square_norm(b, N));
    int count_iterations = 0;
    double res = 1;
    double sum_norm = 0;

    double begin = omp_get_wtime();
    //print_matrix(A);
    //print_vector(b);

    int iter_count = 0;

    int thread_id = omp_get_thread_num();
    for (iter_count = 0; res > EPSILON && iter_count < MAX_ITERATION_COUNT; ++iter_count)
    {
        mul(A + offsets[thread_id] * N, x, buffer + offsets[thread_id], line_counts[thread_id]);
        sub_vectors(buffer + offsets[thread_id], b + offsets[thread_id], buffer + offsets[thread_id], line_counts[thread_id]);

        count_new_x(x + offsets[thread_id], buffer + offsets[thread_id], line_counts[thread_id]);


        sum_norm = 0;

        sum_norm += count_square_norm(buffer + offsets[thread_id], line_counts[thread_id]);

        res = sqrt(sum_norm) / b_norm;
    }

    double end = omp_get_wtime();

    if (count_iterations == MAX_ITERATION_COUNT){
        printf("Wrong tau\n");
    }

    else{
        printf("%f \n", (end - begin));
        //check(A, x, b);
        //print_vector(x);
    }

    free(A);
    free(x);
    free(b);
    free(buffer);

    return EXIT_SUCCESS;
}

void check(double* A, double*x , double* b){
    double * result = malloc(sizeof(double) * N);
    mul(A,x, result, N);
    for (int i = 0; i < N; i++){
        if (result[i] - b[i] < EPSILON || b[i] - result[i] < EPSILON){}
        else{
            printf("not okay\n");
            return;
        }
    }
    printf("correct\n");
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


void print_matrix(const double* matrix)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%f ", matrix[i * N + j]);
        }

        printf("\n");
    }

    printf("\n");
}

void print_vector(const double* vector)
{
    for(int i = 0; i < N; i++)
    {
        printf("%f ", vector[i]);
    }

    printf("\n");
}


double count_square_norm(const double* vector, int size)
{
    double norm_value = 0;

    #pragma omp parallel for schedule(runtime) reduction(+: norm_value)
    for (int i = 0; i < size; i++)
    {
        norm_value += vector[i] * vector[i];
    }

    return norm_value;
}

void set_matrix_part(int* line_counts, int* offsets, int size, int thread_num)
{
    int offset = 0;
    for (int i = 0; i < thread_num; ++i)
    {
        line_counts[i] = size / thread_num;

        if (i < size % thread_num)
        {
            ++line_counts[i];
        }

        offsets[i] = offset;
        offset += line_counts[i];

    }
}

void mul(const double* matrix, const double* vector, double* result, int lines)
{
    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < lines; i++)
    {
        result[i] = 0;

        for (int j = 0; j < N; j++)
        {
            result[i] += matrix[i * N + j] * vector[j];
        }
    }

}

void sub_vectors(const double* vector1, const double* vector2, double* result, int size)
{
    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < size; i++)
    {
        result[i] = vector1[i] - vector2[i];
    }
}

void count_new_x(double* x, const double* vector, int size)
{
    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < size; i++)
    {
        x[i] = x[i] - TAU * vector[i];
    }
}