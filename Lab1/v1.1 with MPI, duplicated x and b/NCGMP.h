#ifndef LAB1_NCGM_H
#define LAB1_NCGM_H
#include "../Matrix.h"
//#include <mpi.h>

void solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps, int * counter);
void parallelMultiplyMatrix(Matrix dest, Matrix first, Matrix second, Matrix tempResult);

#endif //LAB1_NCGM_H
