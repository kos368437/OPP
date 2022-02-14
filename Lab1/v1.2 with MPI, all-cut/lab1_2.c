#include <stdio.h>
#include "../Matrix.h"
#include <time.h>
#include "NCGMP.h"
#include <mpi.h>
#include <malloc.h>


int main(int argc, char * argv[]) {
    
    MPI_Init(&argc, &argv);
    const double eps = 1e-5;
    int rank, size;
    Matrix A;
    Matrix x;
    Matrix b;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    initiate(&A, &b, "../inputA.txt", "../inputb.txt", rank, size);
    if (rank == 0) {
        x = createMatrix(A.height, 1);
    }

    int iterationCounter;

    struct timespec tm_start, tm_finish;
    double best_time = 0;
    double time;

    for (int i = 0; i < 5; ++i) {
        clock_gettime(CLOCK_REALTIME, &tm_start);

        solveSystemUsingNCGM(A, x, b, eps, &iterationCounter);

        clock_gettime(CLOCK_REALTIME, &tm_finish);

        time = tm_finish.tv_sec - tm_start.tv_sec + 1e-9 * (tm_finish.tv_nsec - tm_start.tv_nsec);

        if (best_time == 0 || time < best_time) {
            best_time = time;
        }
    }

    
    if (rank == 0) {
        printf("Iterations count: %d\n", iterationCounter);
        printf("Time elapsed: %lf sec.\n", best_time);

        writeMatrixToFile(x, "output.txt");

        deleteMatrix(x);
    }
    deleteMatrix(A);
    deleteMatrix(b);
    MPI_Finalize();
    return 0;
}

