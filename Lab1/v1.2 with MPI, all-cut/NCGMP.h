#ifndef LAB1_NCGMP_H
#define LAB1_NCGMP_H
#include "../Matrix.h"


void solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps, int * counter);
void parallelMultiplyMatrix(Matrix dest, Matrix first, Matrix second);
double parallelMultiplyVector(Matrix first, Matrix second);

void initiate(Matrix * A, Matrix * b, char * inputA, char * inputb, int commRank, int commSize);
int * initCounts(int N, int M, int commSize);
int * initDispls(int commSize, int * counts);

#endif //LAB1_NCGMP_H
