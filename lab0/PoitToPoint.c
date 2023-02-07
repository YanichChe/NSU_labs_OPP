//mpicc example.c –o example
//mpiexec –n 4 ./
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
    int *currentVector1;
    int* currentVector2;

    long totalResult, currentResult;
    double startTime, endTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &numProc);

    int currentSize = SIZE / (numProc - 1) + 1;

    currentVector1 = malloc(currentSize * sizeof(int));
    currentVector2 = malloc(currentSize * sizeof(int));

    if (rank == 0)
    {
        totalResult = 0;

        vector1 = malloc(SIZE * sizeof(int));
        vector2 = malloc(SIZE * sizeof(int));

        createVector(vector1, SIZE);
        createVector(vector2, SIZE);

        startTime = MPI_Wtime();
        for (int i = 1; i < numProc; i++)
        {
            for (int j = 0; j < currentSize - 1; j++)
            {
                currentVector1[j] = vector1[(i - 1) * currentSize + j];
                currentVector2[j] = vector2[(i - 1) * currentSize + j];
            }

            if (i <=  SIZE % (numProc - 1))
            {
                currentVector1[currentSize] = vector1[(numProc - 1) * (currentSize - 1) + (i - 1)];
            }

            else
            {
                currentVector1[currentSize] = 0;
                currentVector2[currentSize] = 0;
            }

            MPI_Send(currentVector1, currentSize, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(currentVector2, currentSize, MPI_INT, i, 1, MPI_COMM_WORLD);

            MPI_Recv(&currentResult, 1, MPI_LONG,i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            totalResult += currentResult;

        }
    }

    else
    {
        currentResult = 0;

        MPI_Recv(currentVector1, currentSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(currentVector2, currentSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        currentResult = mult(currentVector1, currentVector2, currentSize);

        MPI_Send(&currentResult, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD);

    }

    if (rank == 0)
    {
        //printf("Total result = %ld\n", mult(vector1, vector2, SIZE));
        printf("Total result = %ld\n", totalResult);
        endTime = MPI_Wtime();
        printf("Total time is %f\n", endTime - startTime);

        free(currentVector1);
        free(currentVector2);
        free(vector1);
        free(vector2);
    }

    MPI_Finalize();

    return (0);
}
