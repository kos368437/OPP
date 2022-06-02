#include "BoolMatrix.h"
#include <stdio.h>
#include <stdlib.h>

int realIndex(BoolMatrix boolMatrix, unsigned int row, unsigned int column) {
    return row * boolMatrix.width + column;
}

BoolMatrix createBoolMatrix(unsigned int height, unsigned int width) {
    BoolMatrix boolMatrix;

    boolMatrix.height = height;
    boolMatrix.width = width;
    boolMatrix.arr = (bool *)malloc(height * width * sizeof(bool));

    return boolMatrix;
}

void deleteBoolMatrix(BoolMatrix del) {
    free(del.arr);
}

void fillBoolMatrix(BoolMatrix boolMatrix, bool fillEl) {
    for(int i = 0; i < boolMatrix.height * boolMatrix.width; ++i) {
        boolMatrix.arr[i] = fillEl;
    }
}

void writeBoolMatrixToConsole(BoolMatrix boolMatrix) {

    for (int i = 0; i < boolMatrix.height; i++) {
        for (int j = 0; j < boolMatrix.width; j++) {
            if (get(boolMatrix, i, j)) {
                printf("# ");
            } else {
                printf("O ");
            }
        }
        printf("\n");
    }
}

bool get(BoolMatrix boolMatrix, unsigned int row, unsigned int col) {
    return boolMatrix.arr[realIndex(boolMatrix, row, col)];
}

void set(BoolMatrix boolMatrix, unsigned int row, unsigned int col, double value) {
    boolMatrix.arr[realIndex(boolMatrix, row, col)] = value;
}