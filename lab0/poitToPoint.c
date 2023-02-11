//mpicc example.c –o example
//mpiexec –n 4 ./
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "limits.h"

const int SIZE = 99000;

void createVector(int* vector1, int* vector2, int size)
{
    for (unsigned int i = 0; i < size; i++) {
        vector1[i] = (rand() % 100);
        vector2[i] = (rand() % 100);
    }
}

long long int mult(int* vector1, int* vector2, int size)
{
    long long int result = 0;
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
    int *currentVector;

    long long totalResult, currentResult;
    double startTime, endTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &numProc);

    int currentSize = SIZE / (numProc - 1) + 1;

    currentVector = malloc(currentSize * sizeof(int));
    vector1 = malloc(SIZE * sizeof(int));
    vector2 = malloc(SIZE * sizeof(int));

    if (rank == 0)
    {
        totalResult = 0;

        createVector(vector1, vector2, SIZE);

        startTime = MPI_Wtime();
        for (int i = 1; i < numProc; i++)
        {
            for (int j = 0; j < currentSize - 1; j++)
            {
                currentVector[j] = vector1[(i - 1) * currentSize + j];
            }

            if (i <=  SIZE % (numProc - 1))
            {
                currentVector[currentSize] = vector1[(numProc - 1) * (currentSize - 1) + (i - 1)];
            }

            else
            {
                currentVector[currentSize] = 0;
            }

            MPI_Send(currentVector, currentSize, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(vector2, SIZE, MPI_INT, i, 1, MPI_COMM_WORLD);
        }

        for (int i = 0; i < numProc; i++)
        {
            MPI_Recv(&currentResult, 1, MPI_LONG_LONG,i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            totalResult += currentResult;
        }
    }

    else
    {
        currentResult = 0;

        MPI_Recv(currentVector, currentSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(vector2, SIZE, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        currentResult = mult(currentVector, vector2, currentSize);

        MPI_Send(&currentResult, 1, MPI_LONG_LONG, 0, 1, MPI_COMM_WORLD);

    }

    if (rank == 0)
    {
        printf("Total result = %lld\n", totalResult);
        endTime = MPI_Wtime();
        printf("Total time is %f\n", endTime - startTime);

        free(currentVector);
        free(vector1);
        free(vector2);
    }

    MPI_Finalize();

    return (0);
}
