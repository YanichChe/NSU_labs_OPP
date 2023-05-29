#include <iostream>
#include <mpi.h>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <pthread.h>
#include <cstdio>

struct Task {
    int repeat_num;
};

const int THREAD_NUM = 2;
const int L = 10000;
const int ITERATIONS_COUNT = 8;

pthread_mutex_t mutex;
pthread_t threads[THREAD_NUM];
std::vector<Task> tasks;
MPI_Datatype TASK_TYPE;

int rank;
int size;
int iter_counter;
int task_index;

int get_new_task(int proc_rank);
double calc_task(Task &task);
void generate_tasks(int count_tasks);

void *calc_thread(__attribute__((unused)) void *args);
void *data_thread(__attribute__((unused)) void *args);
int run();

int main(int argc, char **argv) {
    setlocale(LC_ALL, "Russian");

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (provided != MPI_THREAD_MULTIPLE) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    /* Создание нового типа для отправки структуры Task _____________________________________________________*/
    int count = 1;
    int block_lengths[] = {1};
    MPI_Datatype types[] = {MPI_INT};
    MPI_Aint displs = 0;
    MPI_Type_create_struct(count, block_lengths, &displs, types,
                           &TASK_TYPE);
    MPI_Type_commit(&TASK_TYPE);
    /*------------------------------------------------------------------------------------------------------*/

    pthread_mutex_init(&mutex, nullptr);
    double start = MPI_Wtime();

    if (run() != EXIT_SUCCESS) {
        pthread_mutex_destroy(&mutex);
        MPI_Type_free(&TASK_TYPE);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    double end = MPI_Wtime();
    printf("Процесс = %d. Время: %lf\n", rank, end - start);

    pthread_mutex_destroy(&mutex);
    MPI_Type_free(&TASK_TYPE);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

int run() {
    //атрибуты потока
    pthread_attr_t attrs;

    //инициализация атрибутов потока
    if(pthread_attr_init(&attrs) != 0)
    {
        perror("Не могу инициализировать атрибуты");
        return EXIT_FAILURE;
    }

    //установка атрибута "присоединенный"
    if(pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE) != 0)
    {
        perror("Ошибка в установке атрибута");
        abort();
    }

    pthread_create(&threads[0], &attrs, calc_thread, nullptr);
    pthread_create(&threads[1], &attrs, data_thread, nullptr);

    //освобождение ресурсов, занимаемых описателем атрибутов
    pthread_attr_destroy(&attrs);

    for (auto & thread : threads) {
        if (pthread_join(thread, nullptr) != 0) {
            perror("Не могу подсоединиться к потоку");
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

void *calc_thread(__attribute__((unused)) void *args) {
    int *able_get_task_flags = new int[size];
    for (iter_counter = 0; iter_counter < ITERATIONS_COUNT; iter_counter++) {
        task_index = 0;
        double result = 0;
        int task_count = 0;
        int count = 100 * size;

        generate_tasks(count);
        double start = MPI_Wtime();

        for (int i = 0; i < size; i++) {
            able_get_task_flags[i] = 1;
        }
        able_get_task_flags[rank] = 0;

        while (true) {
            while (task_index < tasks.size()) {
                pthread_mutex_lock(&mutex);
                Task task_to_do = tasks[task_index];
                pthread_mutex_unlock(&mutex);

                result += calc_task(task_to_do);

                task_count++;
                task_index++;
            }

            int i = 0;
            for (; i < size;) {
                if (!able_get_task_flags[i]) {
                    i++;
                    continue;
                } else {
                    task_count++;
                    able_get_task_flags[i] = get_new_task(i);
                    break;
                }
            }
            if (i == size) {
                break;
            }
        }

        double end = MPI_Wtime();
        double time = end - start;
        printf("Процесс: %d. Конец итерации %d, время: %lf, "
               "результат: %lf, сделано заданий: %d\n", rank, iter_counter, time, result, task_count);

        MPI_Barrier(MPI_COMM_WORLD);

        double max_time = 0;
        double min_time = 0;
        MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                   MPI_COMM_WORLD);
        MPI_Reduce(&time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0,
                   MPI_COMM_WORLD);

        if(rank == 0) {
            printf("Итерация = %d Время дисбаланса: %lf\n"
                   "Процент дисбаланса: %lf\n",
                   iter_counter, max_time - min_time, (max_time - min_time) / max_time * 100.0);
        }
    }

    int request = 0;

    MPI_Send(&request, 1, MPI_INT, rank, 1, MPI_COMM_WORLD);

    delete[] able_get_task_flags;
    return nullptr;
}

void *data_thread(__attribute__((unused)) void *args) {
    while (iter_counter < ITERATIONS_COUNT) {
        MPI_Status status;
        int res;

        MPI_Recv(&res, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,
                 &status);

        if (res == 0) break;

        if (task_index >= tasks.size()) {
            int answer = 0;
            MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, 2,
                     MPI_COMM_WORLD);
        } else {
            pthread_mutex_lock(&mutex);
            Task task_to_send = tasks.back();
            tasks.pop_back();
            pthread_mutex_unlock(&mutex);

            int answer = 1;
            MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, 2,
                     MPI_COMM_WORLD);
            MPI_Send(&task_to_send, 1, TASK_TYPE, status.MPI_SOURCE, 3,
                     MPI_COMM_WORLD);

        }
    }
    return nullptr;
}

int get_new_task(int proc_rank) {

    int request = 1;

    MPI_Send(&request, 1, MPI_INT, proc_rank, 1, MPI_COMM_WORLD);
    MPI_Recv(&request, 1, MPI_INT, proc_rank, 2, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    if (request == 0) return 0;

    Task task{};
    MPI_Recv(&task, 1, TASK_TYPE, proc_rank, 3, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    pthread_mutex_lock(&mutex);
    tasks.push_back(task);
    pthread_mutex_unlock(&mutex);

    return 1;
}

double calc_task(Task &task) {
    double res = 0;
    for (int i = 0; i < task.repeat_num; i++) {
        res += sin(i);
    }

    return res;
}

void generate_tasks(int count_tasks)
{
    pthread_mutex_lock(&mutex);
    tasks.clear();
    for (int i = 0; i < count_tasks; i++) {
        Task task{};
        task.repeat_num = std::abs(50 - i % 100) * std::abs(rank - (iter_counter % size)) * L;
        tasks.push_back(task);
    }
    pthread_mutex_unlock(&mutex);
}