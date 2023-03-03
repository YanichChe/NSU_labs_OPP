#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#define SIZE 99000

int calcArraySize(int numProc) {
    if (SIZE % numProc == 0) {
        return SIZE;
    }
    else {
        return SIZE + numProc - (SIZE % numProc);
    }
}

void createVector(int* vector1, int* vector2, int size)
{
    for (unsigned int i = 0; i < SIZE; i++) {
        vector1[i] = (rand() % 100);
        vector2[i] = (rand() % 100);
    }

    for(unsigned int i = SIZE; i < size; i++)
    {
        vector1[i] = 0;
        vector1[i] = 0;
    }
}

long long mult(int* vector1, int* vector2, int size) {
    long long result = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < SIZE; ++j) {
            result += vector1[i] * vector2[j];
        }
    }
    return result;
}

int main(int argc, char** argv) {

    int rank, numProc = 0;

    long long totalResult, currentResult;
    double startTime, endTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &numProc);


    const int currentArraySize = calcArraySize(numProc);

    int* vector1 = malloc(currentArraySize* sizeof(int));
    int* vector2 = malloc(currentArraySize * sizeof(int));

    const int currentSize = currentArraySize / numProc;

    if (rank == 0) {
        createVector(vector1, vector2, currentArraySize);

        startTime = MPI_Wtime();

        int shift = currentSize;
        for (int i = 1; i < numProc; ++i) {
            MPI_Send(vector1 + shift, currentSize, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(vector2, SIZE, MPI_INT, i, 1, MPI_COMM_WORLD);
            shift += currentSize;
        }

        totalResult = mult(vector1, vector2, currentSize);
        currentResult = 0;

        for (int i = 1; i < numProc; ++i) {
            MPI_Recv(&currentResult, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            totalResult += currentResult;
        }

    }
    else {
        MPI_Recv(vector1, currentSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(vector2, SIZE, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        currentResult = mult(vector1, vector2, currentSize);
        MPI_Send(&currentResult, 1, MPI_LONG_LONG, 0, 1, MPI_COMM_WORLD);
    }

    if (rank == 0)
    {
        printf("Total result = %lld\n", totalResult);
        endTime = MPI_Wtime();
        printf("Total time is %f\n", endTime - startTime);

        free(vector1);
        free(vector2);
    }

    MPI_Finalize();
    return 0;
}
