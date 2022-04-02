#ifndef LAB1_NCGMP_H
#define LAB1_NCGMP_H
#include "../Matrix.h"


int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps);

void initiate(Matrix *A, Matrix *b, unsigned int N, int commRank, int commSize);
void randomInitiate(Matrix *A, Matrix *b, unsigned int N, int commRank, int commSize);

unsigned int * initCounts(int N, int M, int commSize);
unsigned int * initDispls(int commSize, unsigned int * counts);

Matrix getStandardSymmetricResolvableMatrix(unsigned int N, int commRank, int commSize);
Matrix getStandardResolvableVector(unsigned int N, int commRank, int commSize);
Matrix getStandardRandomResolvableVector(unsigned int N, int commRank, int commSize, Matrix A);
#endif //LAB1_NCGMP_H
