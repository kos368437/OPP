#ifndef LAB4_BOOLMATRIX_H
#define LAB4_BOOLMATRIX_H

#include <stdbool.h>

typedef struct {
    bool * arr;
    unsigned int height;
    unsigned int width;
} BoolMatrix;

BoolMatrix createBoolMatrix(unsigned int height, unsigned int width);
void deleteBoolMatrix(BoolMatrix del);
void fillBoolMatrix(BoolMatrix boolMatrix, bool fillEl);
bool get(BoolMatrix boolMatrix, unsigned int row, unsigned int col);
void set(BoolMatrix boolMatrix, unsigned int row, unsigned int col, bool value);
void writeBoolMatrixToConsole(BoolMatrix boolMatrix);

#endif