#ifndef LAB1_NCGM_H
#define LAB1_NCGM_H
#include "../Matrix.h"
//#include <mpi.h>

int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps);

unsigned int * initCounts(int commSize, unsigned int height, unsigned int width);
unsigned int * initDispls(int commSize, const unsigned int * counts);

Matrix getStandardSymmetricResolvableMatrix(unsigned int N, int commRank, int commSize);
Matrix getStandardResolvableVector(unsigned int N);
Matrix getStandardRandomResolvableVector(unsigned int N, int commRank, int commSize, Matrix A);
Matrix getRandomSymmetricMatrix(unsigned int N, int commRank, int commSize);

#endif //LAB1_NCGM_H
