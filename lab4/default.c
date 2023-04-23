#include <stdio.h>
#include <time.h>
#include <stdlib.h>

const int n1 = 2400;
const int n2 = 4500;
const int n3 = 2400;

void generate_matrix(double* matrix, const int size1, int size2){
    for(int i = 0; i < size1; i++){
        for(int j = 0; j < size2; j++){
            matrix[i * size2 + j] = (double)rand() / RAND_MAX * 20.0 - 10.0;
        }
    }
}

void mul_matrix(const double* matrix1, const double* matrix2, double* result, int size1, int size2, int size3){
    for (int i = 0; i < size1; i++)
        for (int j = 0; j < size3; j++)
            for (int k = 0; k < size2; k++)
                result[i*size3 + j] += matrix1[i*size2 + k] * matrix2[k*size3 + j];
}

int main(int argc, char** argv){

    double* a_matrix;
    double* b_matrix;
    double* c_matrix;

    a_matrix = malloc(sizeof(double) * n1 * n2);
    b_matrix = malloc(sizeof(double) * n2 * n3);
    c_matrix = calloc(n1 * n3, sizeof(double));
    generate_matrix(a_matrix, n1, n2);
    generate_matrix(b_matrix, n2, n3);

    time_t begin = time(NULL);
    mul_matrix(a_matrix, b_matrix, c_matrix, n1, n2, n3);
    time_t end = time(NULL);
    printf("Total time is %ld seconds\n", (end - begin));
    free(a_matrix);
    free(b_matrix);
    free(c_matrix);

    return EXIT_SUCCESS;
}