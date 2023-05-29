#include <iostream>
#include <mpi.h>
#include <cmath>
#include <ctime>
#include <cstdarg>
#include <unistd.h>
#include <vector>
#include <pthread.h>
#include <stdio.h>

struct Task {
    int repeatNum;
};

pthread_mutex_t mutex;
pthread_t threads[2];
std::vector<Task> tasks;
MPI_Datatype TASK_INFO_TYPE;
int rank;
int size;
int currentIter;
int currentTask;
int doneWeight = 0;
int getNewTask(int procRank);

double doTask(Task &task) {
    double res = 0;
    for (int i = 0; i < task.repeatNum; i++) {
        res += sin(i);
    }
    doneWeight += task.repeatNum / 1000;
    return res;
}

void generateTasks(std::vector<Task> &tasks, int countTasks, int currentIter)
{
    int temp = std::abs(rank - (currentIter % size));
    pthread_mutex_lock(&mutex);
    tasks.clear();
    for (int i = 0; i < countTasks; i++) {
        Task task;
        task.repeatNum = std::abs(50 - i%100) * temp * 10000;
        tasks.push_back(task);
    }
    pthread_mutex_unlock(&mutex);
}

void reset(int *arr) {
    for (int i = 0; i < size; i++) {
        arr[i] = 1;
    }
    arr[rank] = 0;
}

void *calcThread(void *args) {
    int *otherProcesses = new int[size];
    for (currentIter = 0; currentIter < 4; currentIter++) {
        double start1 = MPI_Wtime();
        int count = 100*size;
        generateTasks(tasks, count, currentIter);

        reset(otherProcesses);

        currentTask = 0;

        int i = 0;
        int askRank;
        double result = 0;
        int tCount = 0;
        while (true) {
            pthread_mutex_lock(&mutex);
            while (currentTask < tasks.size()) {
                Task taskToDo = tasks[currentTask];
                pthread_mutex_unlock(&mutex);

                result += doTask(taskToDo);
                tCount++;

                currentTask++;
                pthread_mutex_lock(&mutex);
            }
            pthread_mutex_unlock(&mutex);

            for (; i < size;) {
                askRank = (rank + i) % size;
                if (!otherProcesses[askRank]) {
                    i++;
                    continue;
                } else {
                    tCount++;
                    otherProcesses[askRank] = getNewTask(askRank);
                    break;
                }
            }
            if (i == size) {
                break;
            }
        }
        double end1 = MPI_Wtime();
        double time = end1 - start1;
        printf("Process: %d (thread %d). End iteration %d, time spend: %lf, Result: %lf, tasks done: %d\n", rank, 0, currentIter, time, result, tCount);
        MPI_Barrier(MPI_COMM_WORLD);

        double maxTime = 0;
        double minTime = 0;

        MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0,
                   MPI_COMM_WORLD);
        MPI_Reduce(&time, &minTime, 1, MPI_DOUBLE, MPI_MIN, 0,
                   MPI_COMM_WORLD);

        if(rank == 0) {
            printf("(iteration %d)Time of imbalance: %lf\nPercentage of imbalance: %lf\n", currentIter, maxTime - minTime, (maxTime - minTime) /
            maxTime * 100.0);
        }
    }

    int req = 0;

    MPI_Send(&req, 1, MPI_INT, rank, 1, MPI_COMM_WORLD);

    delete[] otherProcesses;
    return NULL;
}

void *dataThread(void *args) {
    while (currentIter < 4) {
        MPI_Status status;
        int res;

        MPI_Recv(&res, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,
                 &status);
        if (res == 0) {
            break;
        }

        pthread_mutex_lock(&mutex);
        if (currentTask >= tasks.size()) {
            pthread_mutex_unlock(&mutex);
            int answer = 0;

            MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, 2,
                     MPI_COMM_WORLD);
        } else {
            Task taskToSend = tasks.back();
            tasks.pop_back();
            pthread_mutex_unlock(&mutex);
            int answer = 1;

            MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, 2,
                     MPI_COMM_WORLD);
            MPI_Send(&taskToSend, 1, TASK_INFO_TYPE, status.MPI_SOURCE, 3,
                     MPI_COMM_WORLD);

        }
    }
    return NULL;
}

int work() {
    pthread_attr_t attrs;
    pthread_attr_init(&attrs);
    pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);

    pthread_create(&threads[0], &attrs, calcThread, NULL);
    pthread_create(&threads[1], &attrs, dataThread, NULL);

    pthread_attr_destroy(&attrs);

    for (int i = 0; i < 2; i++) {
        if (pthread_join(threads[i], NULL) != 0) {
            perror("Cannot join a thread");
            return -1;
        }
    }
    return 0;
}


int getNewTask(int procRank) {

    int req = 1;
    MPI_Send(&req, 1, MPI_INT, procRank, 1, MPI_COMM_WORLD);
    MPI_Recv(&req, 1, MPI_INT, procRank, 2, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    if (req == 0) {
        return 0;
    }

    Task task;
    MPI_Recv(&task, 1, TASK_INFO_TYPE, procRank, 3, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    pthread_mutex_lock(&mutex);
    tasks.push_back(task);
    pthread_mutex_unlock(&mutex);

    return 1;
}

int main(int argc, char **argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (provided != MPI_THREAD_MULTIPLE) {
        printf("Process: %d. Too low MPI thread level provided\n", rank);
        MPI_Finalize();
        return -1;
    }

    int count = 1;
    int blocklengths[] = {1};
    MPI_Datatype types[] = {MPI_INT};
    MPI_Aint displs = 0;
    MPI_Type_create_struct(count, blocklengths, &displs, types,
                           &TASK_INFO_TYPE);
    MPI_Type_commit(&TASK_INFO_TYPE);

    pthread_mutex_init(&mutex, NULL);

    double start = MPI_Wtime();

    if (work() != 0) {
        pthread_mutex_destroy(&mutex);
        MPI_Type_free(&TASK_INFO_TYPE);
        MPI_Finalize();
        perror("work()");
        return -1;
    }

    double end = MPI_Wtime();
    printf("Process: %d. Time: %lf\n", rank, end - start);

    pthread_mutex_destroy(&mutex);
    MPI_Type_free(&TASK_INFO_TYPE);
    MPI_Finalize();
    return 0;
}