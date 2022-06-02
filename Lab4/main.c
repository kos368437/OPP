#include <mpi.h>
#include <stdlib.h>
#include "BoolMatrix.h"
#include "LifeGameParallel.h"

void initGliderGrid(BoolMatrix *grid, unsigned int height, unsigned int width, int commRank, int commSize);

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);
    BoolMatrix grid;
    int height = atoi(argv[1]), width = atoi(argv[2]);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    initGliderGrid(&grid, height, width, rank, size);
    gameOfLifeParallel(grid, rank, size);

    deleteBoolMatrix(grid);

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