#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define N 6000
#define EPSILON pow(10, -5)
#define MAX_ITERATION_COUNT 50000

void generateMatrix(double* matrix);
void generateVector(double* vector);

void printMatrix(const double* matrix);
void printVector(const double* vector);

double countNorm(const double* vector);
double* subVectors(const double* vector1, const double* vector2,double* result);
void mul(const double* matrix, const double* vector, double* result, int shift);
void countNewX(double* x, double* vector);

double tau =  0.01;
int main(int argc, char **argv)
{
    srand (time (NULL));

    double* A = (double*)malloc(sizeof(double) * N * N);
    generateMatrix(A);

    double* x = malloc(sizeof(double) * N);
    generateVector(x);

    double* b = malloc(sizeof(double) * N);
    generateVector(b);

    double* tmpX = malloc(sizeof(double) * N);
    int countIterations = 0;
    double res = 1;
    double prevRes = 0;

    time_t begin = time(NULL);
    while(res > EPSILON)
    {
        mul(A, x, tmpX, N);
        subVectors(tmpX, b, tmpX);//Ax - b

        res = countNorm(tmpX) / countNorm(b);

        countNewX(x, tmpX);
        countIterations++;

        if (countIterations > MAX_ITERATION_COUNT && prevRes > res)
        {
            printf("Wrong tau...\n");
            free(A);
            free(x);
            free(b);
            return EXIT_SUCCESS;
        }
        prevRes = res;
    }

    time_t end = time(NULL);

    printf("Total time is %ld seconds\n", (end - begin));
    free(A);
    free(x);
    free(b);
    free(tmpX);

    return EXIT_SUCCESS;
}

void countNewX(double* x, double* vector)
{
    for (int i = 0; i < N; i++)
    {
        x[i] = x[i] - tau * vector[i];
    }
}

void generateVector(double* vector)
{
    for (int i = 0; i < N; i++)
    {
        vector[i] = (double)rand()/RAND_MAX*10.0-5.0;
    }
}

void generateMatrix(double* matrix)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < i; j++)
        {
            matrix[i*N + j] = matrix[j*N + i];
        }

        for(int j = i; j < N; j++)
        {
            matrix[i*N + j] = (double)rand()/RAND_MAX*2.0-1.0; //float in range -1 to 1
            if(i == j) matrix[i*N + j] = matrix[i*N + j] + N;

        }
    }
}

void printMatrix(const double* matrix)
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

void printVector(const double* vector)
{
    for(int i = 0; i < N; i++)
    {
        printf("%f ", vector[i]);
    }

    printf("\n");
}

void mul(const double* matrix, const double* vector, double* result, int shift)
{
    for (int i = 0; i < shift; i++)
    {
        result[i] = 0;

        for (int j = 0; j < shift; j++)
        {
            result[i] += matrix[i * shift + j] * vector[j];
        }
    }

}

double* subVectors(const double* vector1, const double* vector2, double* result)
{
    for (int i = 0; i < N; i++)
    {
        result[i] = vector1[i] - vector2[i];
    }
}

double countNorm(const double* vector)
{
    double normValue = 0;
    for (int i = 0; i < N; i++)
    {
        normValue += vector[i] * vector[i];
    }

    return sqrt(normValue);
}