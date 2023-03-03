#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "limits.h"

int SIZE = 99000;

void createVector(int* vector1, int* vector2, int size)
{
    for (long long int i = 0; i < size; i++) {
        vector1[i] = rand() % 100;
        vector2[i] = rand() % 100;
    }
}

long long int mult(int* vector1, int* vector2, int size)
{
    long long int result = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++){
            result += vector1[i] * vector2[i];
        }
    }
    return result;
}

int main(int argc, char *argv[])
{
    int* vector1 = malloc(SIZE * sizeof(int));
    int* vector2 = malloc(SIZE * sizeof(int));

    createVector(vector1, vector2, SIZE);

    time_t begin = time(NULL);
    printf("Result %lld\n", mult(vector1, vector2, SIZE));
    time_t end = time(NULL);
    printf("Total time is %ld seconds\n", (end - begin));

    free(vector1);
    free(vector2);

    return 0;
}
