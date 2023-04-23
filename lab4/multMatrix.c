#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"

#define X 0
#define Y 1

const int n1 = 2400;
const int n2 = 5500;
const int n3 = 4800;

int p1;
int p2;

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

void print_matrix(const double* matrix, int row_count, int col_count) {
    for (int i = 0; i < row_count; i++) {
        for (int j = 0; j < col_count; j++)
            printf("%2.4f ", matrix[i * col_count + j]);
        printf("\n");
    }
}


void create_cartesian_lattice(MPI_Comm* cartesian_lattice){
    int dimensions_size[] = {p1, p2};
    int periodic_status[] = {1, 1};

    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions_size, periodic_status,
                    1, cartesian_lattice);

}

void create_sub_comm(MPI_Comm* cartesian_lattice, MPI_Comm* row, MPI_Comm* column){
    int varying_coords[2];
    varying_coords[X] = 0;
    varying_coords[Y] = 1;
    MPI_Cart_sub(*cartesian_lattice, varying_coords, row);

    varying_coords[X] = 1;
    varying_coords[Y] = 0;
    MPI_Cart_sub(*cartesian_lattice, varying_coords, column);
}

void gather_matrix_c(double* c_matrix, double* c_part, int a_part_size, int b_part_size, int proc_number){
    MPI_Datatype block_not_resized, block_resized;

    MPI_Type_vector(a_part_size, b_part_size, n3, MPI_DOUBLE, &block_not_resized);
    MPI_Type_commit(&block_not_resized);

    MPI_Type_create_resized(block_not_resized, 0, b_part_size * sizeof(double), &block_resized);
    MPI_Type_commit(&block_resized);

    int *recv_counts = (int*)malloc(sizeof(int) * proc_number);
    int *displs = (int*)malloc(sizeof(int) * proc_number);

    int block_size = a_part_size * b_part_size;

    for (int i = 0; i < p1; ++i)
        for (int j = 0; j < p2; ++j)
        {
            recv_counts[i * p2 + j] = 1;
            displs[i * p2 + j] = j + i * p2 * a_part_size;
        }

    MPI_Gatherv(c_part, block_size, MPI_DOUBLE, c_matrix,
                recv_counts, displs, block_resized, 0, MPI_COMM_WORLD);

    MPI_Type_free(&block_resized);
    MPI_Type_free(&block_not_resized);
}

int main(int argc, char** argv){

    int proc_number, proc_rank;

    double* a_matrix;
    double* b_matrix;
    double* c_matrix;

    double* a_part;
    double* b_part;
    double* c_part;

    MPI_Comm cartesian_lattice;
    MPI_Comm row;
    MPI_Comm column;

    int a_part_size;
    int b_part_size;

    p1 = atoi(argv[1]);
    p2 = atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_number);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    a_part_size = n1 / p1;
    b_part_size = n3 / p2;

    a_part = (double*)malloc(sizeof(double) * n2 * a_part_size);
    b_part = (double*)malloc(sizeof(double) * n2 * b_part_size);
    c_part = (double*)calloc(a_part_size * b_part_size, sizeof(double));

    double start_time, end_time;

    if (proc_rank == 0){
        a_matrix = (double*)malloc(sizeof(double) * n1 * n2);
        b_matrix = (double*)malloc(sizeof(double) * n2 * n3);
        c_matrix = (double*)calloc(n1 * n3, sizeof(double));
        generate_matrix(a_matrix, n1, n2);
        generate_matrix(b_matrix, n2, n3);
    }

    start_time = MPI_Wtime();
    //create communicators
    create_cartesian_lattice(&cartesian_lattice);

    int coordinates[2];
    MPI_Cart_coords(cartesian_lattice, proc_rank, 2,coordinates);

    create_sub_comm(&cartesian_lattice, &row, &column);


    // distribute data

    if (coordinates[Y] == 0){
        MPI_Scatter(a_matrix, a_part_size * n2, MPI_DOUBLE, a_part,
                    a_part_size * n2, MPI_DOUBLE, 0, column);
    }

    MPI_Bcast(a_part, a_part_size * n2, MPI_DOUBLE, 0, row);

    MPI_Datatype column_not_resized, column_resized;
    MPI_Type_vector(n2, b_part_size, n3, MPI_DOUBLE, &column_not_resized);
    MPI_Type_commit(&column_not_resized);

    MPI_Type_create_resized(column_not_resized, 0, b_part_size * sizeof(double), &column_resized);
    MPI_Type_commit(&column_resized);

    if (coordinates[X] == 0) {
        MPI_Scatter(b_matrix, 1, column_resized, b_part, n2 * b_part_size,
                    MPI_DOUBLE, 0, row);
    }

    MPI_Bcast(b_part, b_part_size * n2, MPI_DOUBLE, 0, column);


    mul_matrix(a_part, b_part, c_part, a_part_size, n2, b_part_size);


    //gather matrix
    gather_matrix_c(c_matrix, c_part, a_part_size, b_part_size, proc_number);

    if (proc_rank == 0){
        end_time = MPI_Wtime();
        /*print_matrix(c_matrix, n1, n3);
        printf("_----------------------------------------------------------------\n");
        double* new_matrix = calloc(n1 * n3, sizeof(double));
        mul_matrix(a_matrix, b_matrix, new_matrix, n1, n2,n3);
        print_matrix(new_matrix, n1, n3);*/
        printf("%d %d\n", p1, p2);
        printf("Total time: %f\n", end_time - start_time);
        free(a_matrix);
        free(b_matrix);
        free(c_matrix);
    }

    free(a_part);
    free(b_part);
    free(c_part);

    MPI_Type_free(&column_resized);
    MPI_Type_free(&column_not_resized);
    MPI_Finalize();
    return EXIT_SUCCESS;
}