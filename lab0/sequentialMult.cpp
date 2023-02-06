#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "limits.h"

const unsigned int SIZE = INT_MAX / 2;

void createVector(int* vector,unsigned int size)
{
    for (unsigned int i = 0; i < size; i++) {
        vector[i] = rand() % 100;
    }
}

long mult(int* vector1, int* vector2, unsigned int size)
{
    long result = 0;
    for (unsigned int i = 0; i < size; i++) {
        result += vector1[i] * vector2[i];
    }
    return result;
}

int main(int argc, char *argv[])
{

    int* vector1 = malloc(SIZE * sizeof(int));
    int* vector2 = malloc(SIZE * sizeof(int));

    time_t begin = time(NULL);

    createVector(vector1, SIZE);
    createVector(vector2, SIZE);

    printf("Result %ld\n", mult(vector1, vector2, SIZE));
    time_t end = time(NULL);
    printf("Total time is %ld seconds\n", (end - begin));

    free(vector1);
    free(vector2);

    return 0;
}
