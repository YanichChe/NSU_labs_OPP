#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 6000
#define EPSILON pow(10, -5)
#define MAX_ITERATION_COUNT 50000
#define TAU  0.01

void setMatrixPart(int* lineCounts, int* offsets, int* sendCounts, int* displs, int size, int numProc);
double countSquareNorm(const double *vector, int size);

void generateVector(double* vector);
void generateMatrix(double* matrix);
void printMatrix(const double* matrix);
void printVector(const double* vector, int size);

int main(int argc, char **argv)
{
    int rank, numProc;
    double startTime, endTime;

    MPI_Init(&argc, &argv);
    startTime = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    int* linesNumber = malloc(sizeof(int) * numProc);
    int* offsets = malloc(sizeof(int) * numProc);
    int* sendCounts = malloc(sizeof(int) * N);
    int* displs = malloc(sizeof(int) * N);
    setMatrixPart(linesNumber, offsets, sendCounts, displs, N, numProc);

    double* partA = malloc(sizeof(double) * linesNumber[rank] * N);
    double* A = malloc(sizeof(double) * N * N);
    double* x = malloc(sizeof(double) * N);
    double* b = malloc(sizeof(double) * N);

    double result;
    double normB;
    if (rank == 0)
    {
        generateMatrix(A);
        generateVector(x);
        generateVector(b);
        result = 1;
        normB = sqrt(countSquareNorm(b, N));
    }

    MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(A, sendCounts, displs, MPI_DOUBLE, partA,
                 sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* buffer = malloc(sizeof(double) * linesNumber[rank]);
    double* partX = malloc(sizeof(double) * linesNumber[rank]);

    double sumNorm = 0;
    int countIteration = 0;

    while (result > EPSILON && countIteration < MAX_ITERATION_COUNT)
    {
        for (int i = 0; i < linesNumber[rank]; ++i)
        {
            for (int j = 0; j < N; ++j)
                buffer[i] += partA[i * N + j] * x[j];
            buffer[i] = buffer[i] - b[offsets[rank] + i];
        }

        for (int i = 0; i < linesNumber[rank]; ++i)
            partX[i] = x[offsets[rank] + i] - TAU * buffer[i];

        MPI_Allgatherv(partX, linesNumber[rank], MPI_DOUBLE,
                       x, linesNumber, offsets, MPI_DOUBLE, MPI_COMM_WORLD);

        double partNorm = countSquareNorm(buffer, linesNumber[rank]);
        MPI_Reduce(&partNorm, &sumNorm, 1,
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            result  = sqrt(sumNorm) / normB;
            countIteration++;
        }
        MPI_Bcast(&countIteration, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (rank == 0)
    {
        if (countIteration < MAX_ITERATION_COUNT)
        {
            //printVector(x, N);
            endTime = MPI_Wtime();
            printf("Total time is %f\n", endTime - startTime);
        }

        else
        {
            printf("WRONG TAU\n");
        }
    }

    free(linesNumber);
    free(offsets);
    free(x);
    free(b);
    free(A);
    free(partA);
    free(buffer);
    free(partX);

    MPI_Finalize();

    return 0;
}


void setMatrixPart(int* lineCounts, int* offsets, int* sendCounts, int* displs, int size, int numProc)
{
    int offset = 0;
    for (int i = 0; i < numProc; ++i)
    {
        lineCounts[i] = size / numProc;

        if (i < size % numProc)
        {
            ++lineCounts[i];
        }

        offsets[i] = offset;
        offset += lineCounts[i];
        sendCounts[i] = lineCounts[i] * N;
        displs[i] = offsets[i] * N;
    }
}

double countSquareNorm(const double *vector, int size)
{
    double norm_square = 0;
    for (int i = 0; i < size; ++i)
        norm_square += vector[i] * vector[i];

    return norm_square;
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
            matrix[i * N + j] = (double)rand() / RAND_MAX*2.0-1.0; //float in range -1 to 1
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

void printVector(const double* vector, int size)
{
    for(int i = 0; i < size; i++)
    {
        printf("%f ", vector[i]);
    }

    printf("\n");
}

