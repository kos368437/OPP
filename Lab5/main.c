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

#define TAG_REQUEST 42
#define TAG_ANSWER 420
#define QUIT_THREAD (-1)
#define MIN_REMAINED_TASKS_TO_START_SHARING 5

void createTaskList(TaskList * taskList, unsigned int newListLength);
void initTaskList(TaskList * taskList, int commRank, int commSize, unsigned int iterNum);
void createAndStartTaskPullThread(pthread_t *thread_taskPull, TaskList * taskList);
double doTasks(TaskList *taskList, double *workTime_return, int commRank, int commSize);
void * sendTasks(void * taskListPointer);
void deleteTaskList(TaskList * taskList);
int askForExtraTasks(TaskList *taskList, int commRank, int commSize);

int main(int argc, char * argv[]) {
    int required = MPI_THREAD_MULTIPLE;
    int provided;
    double result;

    struct timespec tm_start, tm_finish;
    double totalTime;
    double jobTime = 0;
    double minJobTime = 0;
    double maxJobTime = 0;

    int commRank, commSize;
    int quitMessage = QUIT_THREAD;

    TaskList taskList;

    int numOfWaves = atoi(argv[1]);
    int taskListLengthPerProcess = atoi(argv[2]);

    MPI_Init_thread(&argc, &argv, required, &provided);

    pthread_t thread_taskSend;

    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if (required != provided) {
        perror("MPI thread initiation failure");
        MPI_Finalize();
        abort();
    }

    createTaskList(&taskList, taskListLengthPerProcess);

    createAndStartTaskPullThread(&thread_taskSend, &taskList);

    clock_gettime(CLOCK_REALTIME, &tm_start);

    for (int i = 0; i < numOfWaves; i++) {
        if (commRank == 0) printf("\nWave %d\n", i);

        initTaskList(&taskList, commRank, commSize, i);
        MPI_Barrier(MPI_COMM_WORLD);
        result = doTasks(&taskList, &jobTime, commRank, commSize);

        MPI_Allreduce(&jobTime,&minJobTime,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(&jobTime,&maxJobTime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

        if (commRank == 0) {
            printf("MIN Time\n%lf\n", minJobTime);
            printf("MAX Time\n%lf\n", maxJobTime);
            printf("Imbalance\n%lf\n", maxJobTime - minJobTime);
            printf("Imbalance %%\n%lf\n", ((maxJobTime - minJobTime) / maxJobTime) * 100);
        }
    }

    clock_gettime(CLOCK_REALTIME, &tm_finish);

    MPI_Send(&quitMessage, 1, MPI_INT, commRank, TAG_REQUEST, MPI_COMM_WORLD);
    if (0 != pthread_join(thread_taskSend, NULL)) {
        perror("Cannot join a thread");
        MPI_Finalize();
        abort();
    }

    totalTime = tm_finish.tv_sec - tm_start.tv_sec + 1e-9 * (tm_finish.tv_nsec - tm_start.tv_nsec);

    if (commRank == 0) {
        printf("\nTOTAL TIME:\n%lf\n", totalTime);
    }
    deleteTaskList(&taskList);

    MPI_Finalize();

    return 0;
}

double doTasks(TaskList *taskList, double *workTime_return, int commRank, int commSize) {
    double time;
    struct timespec tm_start, tm_finish;
    int i = 0;
    double result = 0;
    int remainedTasksCount;

    pthread_mutex_lock(&(taskList->mutex_taskList));
    remainedTasksCount = taskList->remainedTasksCount;
    pthread_mutex_unlock(&(taskList->mutex_taskList));

    clock_gettime(CLOCK_REALTIME, &tm_start);

    while (1) {
        i = 0;
        while (remainedTasksCount > 0) {
            for (int j = 0; j < taskList->taskList[i]; j++) {
                result += sin(j);
            }
            i++;
            remainedTasksCount -= 1;
        }
        pthread_mutex_lock(&(taskList->mutex_taskList));
        taskList->remainedTasksCount = remainedTasksCount;
        pthread_mutex_unlock(&(taskList->mutex_taskList));
        if (askForExtraTasks(taskList, commRank, commSize) < 1) {
            break;
        }
        else {
            pthread_mutex_lock(&(taskList->mutex_taskList));
            remainedTasksCount = taskList->remainedTasksCount;
            pthread_mutex_unlock(&(taskList->mutex_taskList));
        }
    }

    clock_gettime(CLOCK_REALTIME, &tm_finish);

    time = tm_finish.tv_sec - tm_start.tv_sec + 1e-9 * (tm_finish.tv_nsec - tm_start.tv_nsec);

//    if (time < 1) {
//        printf("IM HERE");
//        if (askForExtraTasks(taskList, commRank, commSize) < 1) {
//
//            /////////////////////////////////////////////////////
//            printf("%d TASK ARRAY: ", commRank);
//            for (int i = 0; i < 10; i++) {
//                printf(" %d ", taskList->taskList[i]);
//            }
//            printf("\n");
//            printf("TIME: %lf\n", time);
//            printf("TaskList remained: %d\n", taskList->remainedTasksCount);
//            /////////////////////////////////////////////////////
//        }
//    }


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

    if (0 != pthread_create(thread_taskPull, &attrs, sendTasks, taskList)) {
        perror("Cannot create a thread");
        MPI_Finalize();
        abort();
    }

    pthread_attr_destroy(&attrs);
}

void * sendTasks(void * taskListPointer) {
    TaskList * taskList = (TaskList *)taskListPointer;

    int requestRank = QUIT_THREAD;
    int sendTasksCount;
    int * sendTaskListPointer;
    int commRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    while (1) {
        MPI_Recv(&requestRank, 1, MPI_INT, MPI_ANY_SOURCE, TAG_REQUEST, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        if (requestRank == QUIT_THREAD) {
            break;
        }
        else {
            pthread_mutex_lock(&(taskList->mutex_taskList));
            if (taskList->remainedTasksCount >= MIN_REMAINED_TASKS_TO_START_SHARING) {

                sendTasksCount = taskList->remainedTasksCount / 2;
                sendTaskListPointer = taskList->taskList + sendTasksCount;
                if (taskList->remainedTasksCount % 2) sendTaskListPointer += 1;
                taskList->remainedTasksCount -= sendTasksCount;
                MPI_Send(&sendTasksCount, 1, MPI_INT, requestRank, TAG_ANSWER, MPI_COMM_WORLD);
                MPI_Send(sendTaskListPointer, sendTasksCount, MPI_INT, requestRank, TAG_ANSWER, MPI_COMM_WORLD);
            }
            else {
                sendTasksCount = 0;
                MPI_Send(&sendTasksCount, 1, MPI_INT, requestRank, TAG_ANSWER, MPI_COMM_WORLD);

            }
            pthread_mutex_unlock(&(taskList->mutex_taskList));

//            /////////////////////////////////////////////////////
//            printf("I, %d, MY TASK ARRAY: ", commRank);
//            for (int i = 0; i < taskList->taskListLength; i++) {
//                printf(" %d ", taskList->taskList[i]);
//            }
//            printf("\n");
//            /////////////////////////////////////////////////////
//            printf("I, %d,  SEND %d tasks TO PROC %d, my REMAIED TASKS %d  \n", commRank, sendTasksCount, requestRank, taskList->remainedTasksCount);
        }

    }

    return NULL;
}

int askForExtraTasks(TaskList *taskList, int commRank, int commSize) {
    int answer = -1;
    int sendRank = commRank;

    for (int i = (commRank + 1) % commSize; i != commRank; i = (i + 1) % commSize) {
        if (i == commRank) break;
        MPI_Send(&sendRank, 1, MPI_INT, i, TAG_REQUEST, MPI_COMM_WORLD);
        MPI_Recv(&answer, 1, MPI_INT, i, TAG_ANSWER, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if (answer > 0) {
            pthread_mutex_lock(&(taskList->mutex_taskList));

            MPI_Recv(taskList->taskList, answer, MPI_INT, i, TAG_ANSWER, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

//            /////////////////////////////////
//            printf("I, %d, I RECEIVE: ", commRank);
//            for (int i = 0; i < taskList->taskListLength; i++) {
//                printf(" %d ", taskList->taskList[i]);
//            }
//            printf("\n");
//            ////////////////////////////////////

            taskList->remainedTasksCount = answer;
            pthread_mutex_unlock(&(taskList->mutex_taskList));

            return answer;
        }
//        else {
//            printf("I, %d, DIDNT GET TASKS FROM PROC %d. Answer: %d\n", commRank, i, answer);
//        }
    }

    return 0;
}