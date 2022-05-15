#include <stdio.h>
#include "../Matrix.h"
#include <time.h>
#include "ParallelMultiplication.h"
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char * argv[]) {
    
    MPI_Init(&argc, &argv);

    const int N1 = atoi(argv[0]), N2 = atoi(argv[1]), N3 = atoi(argv[2]);

    Matrix A;
    Matrix B_trans;
    Matrix C;

    int dims[2]={0,0},periods[2]={0,0},coords[2],reorder = 1;
    int size, rank, sizey, sizex;

    struct timespec tm_start, tm_finish;
    double best_time = 0;
    double time;

    MPI_Comm GridComm;
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    MPI_Dims_create(size,2,dims);
    sizey = dims[0]; sizex = dims[1];

    MPI_Cart_create(MPI_COMM_WORLD,2, dims, periods, reorder,&GridComm);

    MPI_Comm_rank(GridComm, &rank);

    MPI_Comm RowComm, ColComm;
    int subdims[2];

    subdims[0] = 0;
    subdims[1] = 1;
    MPI_Cart_sub(GridComm, subdims, &RowComm);

    subdims[0] = 1;
    subdims[1] = 0;
    MPI_Cart_sub(GridComm, subdims, &ColComm);

    if (rank == 0) {
        randomInitiate(&A, &B_trans, &C, N1, N2, N3);
    }

    clock_gettime(CLOCK_REALTIME, &tm_start);

    parallelMatrixMultiplication(C, A, B_trans, sizey, sizex, dims, GridComm, RowComm, ColComm, N1, N2, N3);

    clock_gettime(CLOCK_REALTIME, &tm_finish);

    time = tm_finish.tv_sec - tm_start.tv_sec + 1e-9 * (tm_finish.tv_nsec - tm_start.tv_nsec);

    if (best_time == 0 || time < best_time) {
        best_time = time;
    }

    if (rank == 0) {
        printf("Time elapsed: %lf sec.\n", best_time);

        // writeMatrixToConsole(A);
        // writeMatrixToConsole(B_trans);
        // writeMatrixToConsole(C);

        deleteMatrix(A);
        deleteMatrix(B_trans);
        deleteMatrix(C);
    }

    MPI_Finalize();
    return 0;
}

