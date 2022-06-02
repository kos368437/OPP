#include "LifeGameParallel.h"
#include <mpi.h>
#include <malloc.h>
#include <unistd.h>
#include "BoolMatrix.h"

const int UPPER_TAG = 42;
const int LOWER_TAG = 420;

void updateGridState(BoolMatrix * activeGrid, BoolMatrix * bufferGrid, BoolMatrix upperNeighbourRowBuffer, BoolMatrix lowerNeighbourRowBuffer, int commRank, int commSize);
bool newCellState(BoolMatrix * activeGrid, int row, int col, BoolMatrix upperEdgeRow, BoolMatrix lowerEdgeRow);
void swap(BoolMatrix ** first, BoolMatrix ** second);
void updateEdgeRow(BoolMatrix * activeGrid, BoolMatrix  * bufferGrid, unsigned int row, BoolMatrix upperNeighbourRowBuffer, BoolMatrix lowerNeighbourRowBuffer, MPI_Request * request);
bool isSame(BoolMatrix fMatrix, BoolMatrix sMatrix);
bool isTimeToStop(bool * stopFlags, int commSize);

void gameOfLifeParallel(BoolMatrix grid, int commRank, int commSize) {

    bool stopFlags[commSize];
    bool stopFlag;
    BoolMatrix * activeGrid;
    BoolMatrix * bufferGrid;
    BoolMatrix gridBuffer = createBoolMatrix(grid.height, grid.width);
    BoolMatrix upperNeighbourRow = createBoolMatrix(1, grid.width);
    BoolMatrix lowerNeighbourRow = createBoolMatrix(1, grid.width);

    activeGrid = &grid;
    bufferGrid = &gridBuffer;

    do {
        updateGridState(activeGrid, bufferGrid, upperNeighbourRow, lowerNeighbourRow, commRank, commSize);
        stopFlag = isSame(*activeGrid, *bufferGrid);
        MPI_Allgather(&stopFlag, 1, MPI_C_BOOL, stopFlags, 1, MPI_C_BOOL, MPI_COMM_WORLD);

        if (commRank == 0) printf("\n");
        for (int i = 0; i < commSize; i++) {
            if (i == commRank) {
                writeBoolMatrixToConsole(*activeGrid);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }


        swap(&activeGrid, &bufferGrid);
        sleep(1);
    } while (!isTimeToStop(stopFlags, commSize));

    deleteBoolMatrix(gridBuffer);
    deleteBoolMatrix(lowerNeighbourRow);
    deleteBoolMatrix(upperNeighbourRow);
}

void updateGridState(BoolMatrix * activeGrid, BoolMatrix * bufferGrid, BoolMatrix upperNeighbourRowBuffer, BoolMatrix lowerNeighbourRowBuffer, int commRank, int commSize) {
    MPI_Request upperSendRequest;
    MPI_Request lowerSendRequest;
    MPI_Request upperRecvRequest;
    MPI_Request lowerRecvRequest;
    MPI_Status status;

    MPI_Isend(activeGrid->arr, activeGrid->width, MPI_C_BOOL, ((commRank - 1 + commSize) % commSize), LOWER_TAG, MPI_COMM_WORLD, &lowerSendRequest);
    MPI_Isend(activeGrid->arr + (activeGrid->height - 1) * activeGrid->width, activeGrid->width, MPI_C_BOOL, ((commRank + 1 + commSize) % commSize), UPPER_TAG, MPI_COMM_WORLD, &upperSendRequest);

    MPI_Irecv(upperNeighbourRowBuffer.arr, activeGrid->width, MPI_C_BOOL, ((commRank - 1 + commSize) % commSize), UPPER_TAG, MPI_COMM_WORLD, &upperRecvRequest);
    MPI_Irecv(lowerNeighbourRowBuffer.arr, activeGrid->width, MPI_C_BOOL, ((commRank + 1 + commSize) % commSize), LOWER_TAG, MPI_COMM_WORLD, &lowerRecvRequest);

    for (int i = 1; i < activeGrid->height - 1; i++) {
        for (int j = 0; j < activeGrid->width; j++) {
            set(*bufferGrid, i, j, newCellState(activeGrid, i, j, upperNeighbourRowBuffer, lowerNeighbourRowBuffer));
        }
    }

    updateEdgeRow(activeGrid, bufferGrid, 0, upperNeighbourRowBuffer, lowerNeighbourRowBuffer, &upperRecvRequest);
    updateEdgeRow(activeGrid, bufferGrid, activeGrid->height - 1, upperNeighbourRowBuffer, lowerNeighbourRowBuffer, &lowerRecvRequest);

    MPI_Wait(&upperSendRequest, &status);
    MPI_Wait(&lowerSendRequest, &status);
}

void swap(BoolMatrix ** first, BoolMatrix ** second) {
    BoolMatrix * temp;
    temp = *first;
    *first = *second;
    *second = temp;
}

bool newCellState(BoolMatrix * activeGrid, int row, int col, BoolMatrix upperEdgeRow, BoolMatrix lowerEdgeRow) {

    int alive = 0;
    bool newState;
    int activeCol, activeRow;
    BoolMatrix * currentBuffer;

    for (int i = row - 1; i <= row + 1; i++) {
        activeRow = i;
        currentBuffer = activeGrid;
        if (i < 0) {
            currentBuffer = &upperEdgeRow;
            activeRow = 0;
        }
        if (i == activeGrid->height) {
            currentBuffer = &lowerEdgeRow;
            activeRow = 0;
        }
        for (int j = col - 1; j <= col + 1; j++) {
            activeCol = (j + activeGrid->width) % activeGrid->width;

            if ((i == row) && (j == col)) continue;

            if (get(*currentBuffer, activeRow, activeCol)) {
                alive++;
            }
        }
    }

    if (newState = get(*activeGrid, row, col)) {
        if ((alive < 2) || (alive > 3)) {
            newState = false;
        }
    }
    else {
        if (alive == 3) {
            newState = true;
        }
    }

    return newState;
}

BoolMatrix createParallelGrid(unsigned int height, unsigned int width, int commRank, int commSize) {

    unsigned int localHeight = height / commSize;

    if (commRank < height % commSize) {
        localHeight += 1;
    }

    return createBoolMatrix(localHeight, width);
}

void globalSet(BoolMatrix grid, unsigned int globalHeight, unsigned int row, unsigned int col, bool value, int commRank, int commSize) {

    if (commRank == getRowOwnerRank(globalHeight, row, commSize)) {
        set(grid, getLocalRowNumber(row, globalHeight, commRank, commSize), col, value);
    }
}
int getRowOwnerRank(unsigned int globalHeight, unsigned int row, int commSize) {
    int owner;

    if (row < ((globalHeight / commSize) + 1) * (globalHeight % commSize)) {
        owner = row / ((globalHeight / commSize) + 1);
    }
    else {
        row = row - ((globalHeight / commSize) + 1) * (globalHeight % commSize);
        owner = (globalHeight % commSize) + (row / (globalHeight / commSize));
    }

    return owner;
}

unsigned int getRowNumberInGlobalGrid(unsigned int rowNumber, unsigned int globalGridHeight, int commRank, int commSize) {
    unsigned int height = globalGridHeight / commSize;
    unsigned int globalRowNumber;

    if (commRank < globalGridHeight % commSize) {
        height += 1;

        globalRowNumber = commRank * height + rowNumber;
    }
    else {
        globalRowNumber = commRank * height + (globalGridHeight % commSize) + rowNumber;
    }

    return globalRowNumber;
}

unsigned int getLocalRowNumber(unsigned int globalRowNumber, unsigned int globalGridHeight, int commRank, int commSize) {
    int localRowNumber;

    if (globalRowNumber < ((globalGridHeight / commSize) + 1) * (globalGridHeight % commSize)) {
        localRowNumber = globalRowNumber % ((globalGridHeight / commSize) + 1);
    }
    else {
        localRowNumber = globalRowNumber - ((globalGridHeight / commSize) + 1) * (globalGridHeight % commSize);
        localRowNumber = localRowNumber % (globalGridHeight / commSize);
    }

    return localRowNumber;
}

void updateEdgeRow(BoolMatrix * activeGrid, BoolMatrix  * bufferGrid, unsigned int row, BoolMatrix upperNeighbourRowBuffer, BoolMatrix lowerNeighbourRowBuffer, MPI_Request * request) {
    MPI_Status status;
    MPI_Wait(request, &status);
    for (int j = 0; j < activeGrid->width; j++) {
        set(*bufferGrid, row, j, newCellState(activeGrid, row, j, upperNeighbourRowBuffer, lowerNeighbourRowBuffer));
    }
}

bool isSame(BoolMatrix fMatrix, BoolMatrix sMatrix) {

    if ((fMatrix.height != sMatrix.height) || (fMatrix.width != sMatrix.width)) return false;

    for (int i = 0; i < fMatrix.height; i++) {
        for (int j = 0; j < fMatrix.width; j++) {
            if (get(fMatrix, i, j) != get(sMatrix, i, j)) return false;
        }
    }

    return true;
}

bool isTimeToStop(bool * stopFlags, int commSize) {
    for (int i = 0; i < commSize; i++) {
        if (!stopFlags[i]) return false;
    }
    return true;
}