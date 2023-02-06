#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "limits.h"

const int SIZE = INT_MAX / 2;

void createVector(int* vector, int size)
{
    for (int i = 0; i < size; i++) {
        vector[i] = (rand() % 100);
    }
}

long mult(const int* vector1, const int* vector2, int size)
{
    long result = 0;
    for (int i = 0; i < size; i++) {
        result += vector1[i] * vector2[i];
    }
    return result;
}

int main(int argc, char *argv[])
{
    int rank, numProc = 0;

    int* vector1;
    int* vector2;

    int* currentVector1;
    int* currentVector2;

    long totalResult, currentResult;
    double startTime, endTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &numProc);

    int currentSize = SIZE / numProc;

    if (rank == 0)
    {
        totalResult = 0;

        vector1 = malloc(SIZE * sizeof(int));
        vector2 = malloc(SIZE * sizeof(int));

        startTime = MPI_Wtime();
        createVector(vector1, SIZE);
        createVector(vector2, SIZE);
        //printf("Result %ld\n", mult(vector1, vector2, SIZE));
    }

    currentVector1 = malloc(currentSize * sizeof(int));
    currentVector2 = malloc(currentSize * sizeof(int));

    MPI_Scatter(vector1, currentSize, MPI_INT,
                currentVector1, currentSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(vector2, currentSize, MPI_INT,
                currentVector2, currentSize, MPI_INT, 0, MPI_COMM_WORLD);

    currentResult = mult(currentVector1, currentVector2, currentSize);

    MPI_Reduce(&currentResult, &totalResult, 1,
               MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("Total result = %ld\n", totalResult);
        endTime = MPI_Wtime();
        printf("Total time is %f\n", endTime - startTime);
    }


    MPI_Finalize();

    return (0);
}
