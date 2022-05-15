#ifndef LAB1_NCGM_H
#define LAB1_NCGM_H
#include "Matrix.h"

int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps);

Matrix getStandardSymmetricResolvableMatrix(unsigned int N);
Matrix getStandardResolvableVector(unsigned int N);
Matrix getStandardRandomResolvableVector(unsigned int N, Matrix A);
Matrix getRandomSymmetricMatrix(unsigned int N);
#endif //LAB1_NCGM_H
