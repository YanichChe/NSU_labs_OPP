#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "limits.h"

const unsigned int SIZE = 99000;

void createVector(int* vector1, int* vector2, int size)
{
    for (unsigned int i = 0; i < size; i++) {
        vector1[i] = (rand() % 100);
        vector2[i] = (rand() % 100);
    }
}

long long mult(const int* vector1, const int* vector2, int size)
{
    long long result = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < SIZE; j++)
        {
            result += vector1[i] * vector2[j];
        }
    }
    return result;
}

int main(int argc, char *argv[])
{
    int rank, numProc = 0;

    int* vector1;
    int* vector2;

    int* currentVector;

    long long totalResult, currentResult;
    double startTime, endTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &numProc);

    int currentSize = SIZE / numProc;
    vector2 = malloc(SIZE * sizeof(int));

    if (rank == 0)
    {
        totalResult = 0;

        vector1 = malloc(SIZE * sizeof(int));
        createVector(vector1, vector2, SIZE);

        startTime = MPI_Wtime();
    }

    currentVector = malloc(currentSize * sizeof(int));

    MPI_Scatter(vector1, currentSize, MPI_INT,
                currentVector, currentSize, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(vector2 , SIZE, MPI_INT, 0, MPI_COMM_WORLD);

    currentResult = mult(currentVector, vector2, currentSize);

    MPI_Reduce(&currentResult, &totalResult, 1,
               MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Total result = %lld\n", totalResult);
        endTime = MPI_Wtime();
        printf("Total time is %f\n", endTime - startTime);
    }

    MPI_Finalize();

    return (0);
}
