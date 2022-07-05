#ifndef LAB1_NCGMP_H
#define LAB1_NCGMP_H
#include "../Matrix.h"


int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps);

void initiate(Matrix *A, Matrix *b, unsigned int N, int commRank, int commSize);
void randomInitiate(Matrix *A, Matrix *b, unsigned int N, int commRank, int commSize);

int * initCounts(int N, int M, int commSize);
int * initDispls(int commSize, int *counts);

Matrix getStandardSymmetricResolvableMatrix(unsigned int N, int commRank, int commSize);
Matrix getStandardResolvableVector(unsigned int N, int commRank, int commSize);
Matrix getStandardRandomResolvableVector(unsigned int N, int commRank, int commSize, Matrix A);
Matrix getRandomSymmetricMatrix(unsigned int N, int commRank, int commSize);
#endif //LAB1_NCGMP_H
