
#include <mpi.h>
#include <stdlib.h>
#include "BoolMatrix.h"
#include "LifeGameParallel.h"
#include <time.h>
#include <stdio.h>

void initGliderGrid(BoolMatrix *grid, unsigned int height, unsigned int width, int commRank, int commSize);
void initOrdinairyGrid(BoolMatrix *grid, unsigned int height, unsigned int width, int commRank, int commSize);

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);
    BoolMatrix grid;
    int height = atoi(argv[1]), width = atoi(argv[2]);
    int rank, size;
    int iterationsToStop = atoi(argv[3]);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


     struct timespec tm_start, tm_finish;
    double best_time = 0;
    double time;
    int iterCount = 0;

    //initGliderGrid(&grid, height, width, rank, size);
    initOrdinairyGrid(&grid, height, width, rank, size);

        clock_gettime(CLOCK_REALTIME, &tm_start);

    iterCount = gameOfLifeParallel(grid, rank, size, iterationsToStop);

        clock_gettime(CLOCK_REALTIME, &tm_finish);

    time = tm_finish.tv_sec - tm_start.tv_sec + 1e-9 * (tm_finish.tv_nsec - tm_start.tv_nsec);

    if (best_time == 0 || time < best_time) {
        best_time = time;
    }

    if (rank == 0) {
        printf("%d\n", iterCount);
        printf("%lf\n", best_time);
    }

    deleteBoolMatrix(grid);
MPI_Finalize();
    return 0;
}

void initGliderGrid(BoolMatrix *grid, unsigned int height, unsigned int width, int commRank, int commSize) {

    *grid = createParallelGrid(height, width, commRank, commSize);
    fillBoolMatrix(*grid, false);

    globalSet(*grid, height, 0, 1, true, commRank, commSize);
    globalSet(*grid, height, 1, 2, true, commRank, commSize);
    globalSet(*grid, height, 2, 0, true, commRank, commSize);
    globalSet(*grid, height, 2, 1, true, commRank, commSize);
    globalSet(*grid, height, 2, 2, true, commRank, commSize);
}

void initOrdinairyGrid(BoolMatrix *grid, unsigned int height, unsigned int width, int commRank, int commSize) {
    *grid = createParallelGrid(height, width, commRank, commSize);
    bool val = false;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if ((i % (j+1)) % 2) {
                val = true;
            }
            else {
                val = false;
            }
            globalSet(*grid, height, i, j, val, commRank, commSize);
        }
    }
}
