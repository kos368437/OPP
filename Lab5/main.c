#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <mpi.h>
#include <math.h>

typedef struct s_TaskList {
    int * taskList;
    unsigned int remainedTasksCount;
    unsigned int taskListLength;
    pthread_mutex_t mutex_taskList;
}TaskList;

void createTaskList(TaskList * taskList, unsigned int newListLength);
void initTaskList(TaskList * taskList, int commRank, int commSize, unsigned int iterNum);
void createAndStartTaskPullThread(pthread_t *thread_taskPull, TaskList * taskList);
double doTasks(TaskList *taskList, double *workTime_return);
void * pullTasks(void * taskListPointer);
void deleteTaskList(TaskList * taskList);

int main(int argc, char * argv[]) {
    int required = MPI_THREAD_MULTIPLE;
    int provided;
    double result;

    struct timespec tm_start, tm_finish;
    double best_time = 0;
    double totalTime;
    double jobTime = 0;
    int commRank, commSize;

    TaskList taskList;

    int numOfWaves = atoi(argv[1]);
    int taskListLengthPerProcess = atoi(argv[2]);

    MPI_Init_thread(&argc, &argv, required, &provided);

    pthread_t thread_taskPull;

    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if (required != provided) {
        perror("MPI thread initiation failure");
        MPI_Finalize();
        abort();
    }

    createTaskList(&taskList, taskListLengthPerProcess);

    createAndStartTaskPullThread(&thread_taskPull, &taskList);

    clock_gettime(CLOCK_REALTIME, &tm_start);

    for (int i = 0; i < numOfWaves; i++) {
        initTaskList(&taskList, commRank, commSize, i);
        result = doTasks(&taskList, &jobTime);
        if (commRank == 0) {
            printf("Wave %d\n %lf\n", i, jobTime);
        }
    }

    clock_gettime(CLOCK_REALTIME, &tm_finish);

    if (0 != pthread_join(thread_taskPull, NULL)) {
        perror("Cannot join a thread");
        MPI_Finalize();
        abort();
    }

    totalTime = tm_finish.tv_sec - tm_start.tv_sec + 1e-9 * (tm_finish.tv_nsec - tm_start.tv_nsec);

    if (commRank == 0) {
        printf("Total time:\n%lf\n", totalTime);
    }
    deleteTaskList(&taskList);

    MPI_Finalize();

    return 0;
}

double doTasks(TaskList *taskList, double *workTime_return) {
    double time;
    struct timespec tm_start, tm_finish;
    int i = 0;
    double result = 0;

    clock_gettime(CLOCK_REALTIME, &tm_start);

    while (taskList->remainedTasksCount > 0) {
        for (int j = 0; j < taskList->taskList[i]; j++) {
            result += sin(j);
        }
        i++;
        pthread_mutex_lock(&(taskList->mutex_taskList));
        taskList->remainedTasksCount -= 1;
        pthread_mutex_unlock(&(taskList->mutex_taskList));
    }

    clock_gettime(CLOCK_REALTIME, &tm_finish);

    time = tm_finish.tv_sec - tm_start.tv_sec + 1e-9 * (tm_finish.tv_nsec - tm_start.tv_nsec);

    *workTime_return = time;
    return result;
}

void createTaskList(TaskList * taskList, unsigned int newListLength) {
    pthread_mutex_init(&(taskList->mutex_taskList), NULL);
    pthread_mutex_lock(&(taskList->mutex_taskList));
    taskList->taskList = (int *)malloc(taskList->taskListLength * sizeof(int));
    taskList->remainedTasksCount = newListLength;
    taskList->taskListLength = newListLength;
    pthread_mutex_unlock(&(taskList->mutex_taskList));
}

void initTaskList(TaskList * taskList, int commRank, int commSize, unsigned int iterNum) {
    pthread_mutex_lock(&(taskList->mutex_taskList));
    taskList->remainedTasksCount = taskList->taskListLength;
    for (int i = 0; i < taskList->taskListLength; i++) {
        taskList->taskList[i] = abs(50 - i % 100) * abs(commRank - (iterNum % commSize)) * 10000;
//        taskList->taskList[i] = (1 + i % 100) * abs(commRank - (iterNum % commSize));
    }
    pthread_mutex_unlock(&(taskList->mutex_taskList));
}

void deleteTaskList(TaskList * taskList) {
    pthread_mutex_lock(&(taskList->mutex_taskList));
    free(taskList->taskList);
    pthread_mutex_unlock(&(taskList->mutex_taskList));
    pthread_mutex_destroy(&(taskList->mutex_taskList));
}

void createAndStartTaskPullThread(pthread_t *thread_taskPull, TaskList * taskList) {
    pthread_attr_t attrs;

    if (0 != pthread_attr_init(&attrs)) {
        perror("Cannot initialize attributes");
        MPI_Finalize();
        abort();
    };

    if (0 != pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE)) {
        perror("Error in setting attributes");
        MPI_Finalize();
        abort();
    }

    if (0 != pthread_create(thread_taskPull, &attrs, pullTasks, NULL)) {
        perror("Cannot create a thread");
        MPI_Finalize();
        abort();
    }

    pthread_attr_destroy(&attrs);
}
void * pullTasks(void * taskListPointer) {
    TaskList * taskList = (TaskList *)taskListPointer;
}