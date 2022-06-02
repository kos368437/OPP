#ifndef LAB_4_LIFEGAMEPARALLEL_H
#define LAB_4_LIFEGAMEPARALLEL_H

#include "BoolMatrix.h"

void gameOfLifeParallel(BoolMatrix grid, int commRank, int commSize);
void globalSet(BoolMatrix grid, unsigned int globalHeight, unsigned int row, unsigned int col, bool value, int commRank, int commSize);
int getRowOwnerRank(unsigned int globalHeight, unsigned int row, int commSize);
BoolMatrix createParallelGrid(unsigned int height, unsigned int width, int commRank, int commSize);
unsigned int getLocalRowNumber(unsigned int globalRowNumber, unsigned int globalGridHeight, int commRank, int commSize);

#endif //LAB_4_LIFEGAMEPARALLEL_H
