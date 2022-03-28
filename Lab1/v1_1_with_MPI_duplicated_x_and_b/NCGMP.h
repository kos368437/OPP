#ifndef LAB1_NCGM_H
#define LAB1_NCGM_H
#include "../Matrix.h"
//#include <mpi.h>

int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps);

unsigned int * initCounts(int commSize, unsigned int height, unsigned int width);
unsigned int * initDispls(int commSize, const unsigned int * counts);

#endif //LAB1_NCGM_H
