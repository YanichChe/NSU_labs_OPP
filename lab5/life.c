#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int ROWS;
int COLS;
#define STEPS 1000

void init_grid(int *grid) {
    grid[1] = 1;
    grid[COLS + 2] = 1;
    grid[COLS * 2 + 0] = 1;
    grid[COLS * 2 + 1] = 1;
    grid[COLS * 2 + 2] = 1;
}

void fill(int i, const int *curr_grid, int *next_grid) {
    for (int j = 0; j < COLS; j++) {
        int live_neighbors = 0;
        if (curr_grid[i * COLS + (j + 1) % COLS] == 1) live_neighbors++; // right
        if (curr_grid[i * COLS + (j - 1 + COLS) % COLS] == 1) live_neighbors++; // left

        if (curr_grid[(i + 1) * COLS + j] == 1) live_neighbors++; // down
        if (curr_grid[(i - 1) * COLS + j] == 1) live_neighbors++; // up

        if (curr_grid[(i - 1) * COLS + (j + 1) % COLS] == 1) live_neighbors++;// up right
        if (curr_grid[(i - 1) * COLS + (j - 1 + COLS) % COLS] == 1) live_neighbors++;// up left

        if (curr_grid[(i + 1) * COLS + (j - 1 + COLS) % COLS] == 1) live_neighbors++; // down left
        if (curr_grid[(i + 1) * COLS + (j + 1) % COLS] == 1) live_neighbors++; // down right

        if (curr_grid[i * COLS + j] == 1) { // cell is alive
            if (live_neighbors < 2 || live_neighbors > 3) {
                next_grid[i * COLS + j] = 0; // cell dies
            } else {
                next_grid[i * COLS + j] = 1; // cell stays alive
            }
        } else { // cell is dead
            if (live_neighbors == 3) {
                next_grid[i * COLS + j] = 1; // cell becomes alive
            } else {
                next_grid[i * COLS + j] = 0; // cell stays dead
            }
        }

    }
}

int get_flag(int index, int rows_per_proc, const int *copy_array, const int *curr_grid) {

    for (int i = 0; i < index; i++) {
        int equals = 1;
        for (int k = 1; k <= rows_per_proc; k++) {
            for (int j = 0; j < COLS; j++) {
                if (copy_array[i * ROWS * COLS + k * COLS + j] != curr_grid[k * COLS + j]) {
                    equals = 0;
                    break;
                }
            }
        }

        if (equals) { return 1; }

    }


    return 0;
}

void set_matrix_part(int *send_counts, int *displs, int size, int num_proc) {
    int offset = 0;
    for (int i = 0; i < num_proc; ++i) {
        int line_counts = size / num_proc;

        if (i < size % num_proc) {
            ++line_counts;

        }

        displs[i] = offset * COLS;
        offset += line_counts;
        send_counts[i] = line_counts * COLS;

    }
}

bool check_flags(int *recv_buffer_flags, int size, int index) {
    for (int i = 0; i < size; i++) {
        if (recv_buffer_flags[i * STEPS + (index)] == 0) {
            return false;
        }
    }
    return true;
}

void create_cart_comm(int size, MPI_Comm *cart_comm) {
    int dims[2] = {size, 1};
    int periods[2] = {true, false};

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, cart_comm);
}

int main(int argc, char **argv) {

    int rank, size, step;
    int index = 0;

    int *flags;
    int *copy_array;
    int *recv_buffer_flags;
    int *grid;
    int *send_counts;
    int *displs;

    int rows_per_proc, real_size;

    MPI_Comm cart_comm;

    double start_time, end_time;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    ROWS = atoi(argv[1]);
    COLS = atoi(argv[2]);

    send_counts = malloc(sizeof(int) * size);
    displs = malloc(sizeof(int) * size);

    set_matrix_part(send_counts, displs, ROWS, size);

    rows_per_proc = send_counts[rank] / COLS;
    real_size = rows_per_proc + 2;

    if (rank == 0) grid = (int *) calloc(ROWS * COLS, sizeof(int));

    int *curr_grid = (int *) malloc(sizeof(int) * real_size * COLS);
    int *next_grid = (int *) malloc(sizeof(int) * real_size * COLS);

    recv_buffer_flags = (int *) malloc(sizeof(int) * size * STEPS);
    copy_array = (int *) malloc(sizeof(int) * STEPS * ROWS * COLS);

    flags = calloc(STEPS, sizeof(int));

    srand(rank);
    if (rank == 0) init_grid(grid);

    MPI_Scatterv(grid, send_counts, displs, MPI_INT, (curr_grid + COLS),
                 send_counts[rank], MPI_INT, 0, MPI_COMM_WORLD);

    create_cart_comm(size, &cart_comm);

// Determine neighbor ranks
    int left_rank, right_rank;
    MPI_Cart_shift(cart_comm, 0, 1, &left_rank, &right_rank);

    start_time = MPI_Wtime();

    for (step = 0; step < STEPS; step++) {
        MPI_Request reqs[5];
        MPI_Isend(curr_grid + COLS, COLS,
                  MPI_INT, left_rank,
                  0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(curr_grid + (real_size - 2) * COLS, COLS, MPI_INT, right_rank,
                  0, MPI_COMM_WORLD, &reqs[1]);

        MPI_Irecv(curr_grid + (real_size - 1) * COLS, COLS, MPI_INT, right_rank,
                  0, MPI_COMM_WORLD, &reqs[3]);
        MPI_Irecv(curr_grid, COLS, MPI_INT, left_rank,
                  0, MPI_COMM_WORLD, &reqs[2]);

        flags[index] = get_flag(index, rows_per_proc, copy_array, curr_grid);

        MPI_Iallgather(flags, STEPS, MPI_INT, recv_buffer_flags, STEPS,
                       MPI_INT, MPI_COMM_WORLD, &reqs[4]);

        for (int i = 2; i <= real_size - 3; i++) fill(i, curr_grid, next_grid);

        MPI_Wait(&reqs[0], &status);
        MPI_Wait(&reqs[2], &status);

        fill(1, curr_grid, next_grid);

        MPI_Wait(&reqs[1], &status);
        MPI_Wait(&reqs[3], &status);

        fill(real_size - 2, curr_grid, next_grid);

        MPI_Wait(&reqs[4], &status);
        if (check_flags(recv_buffer_flags, size, index)) break;

        for (int i = 1; i <= rows_per_proc; i++)
            for (int j = 0; j < COLS; j++)
                copy_array[index * ROWS * COLS + i * COLS + j] = curr_grid[i * COLS + j];

        index++;

        int *temp = curr_grid;
        curr_grid = next_grid;
        next_grid = temp;
    }

    if (rank == 0) {
        end_time = MPI_Wtime();
        printf("Total steps: %d\n", index);
        printf("Total time: %f\n", end_time - start_time);
        free(grid);
    }

    free(curr_grid);
    free(next_grid);
    free(flags);
    free(copy_array);
    free(recv_buffer_flags);
    free(send_counts);
    free(displs);

    MPI_Finalize();
    return 0;
}

