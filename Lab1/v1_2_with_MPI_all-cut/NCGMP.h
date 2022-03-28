#ifndef LAB1_NCGMP_H
#define LAB1_NCGMP_H
#include "../Matrix.h"


int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps);

void initiate(Matrix *A, Matrix *b, char *inputA, char *inputb, int commRank, int commSize, unsigned int *counts,
              unsigned int *displs);

unsigned int * initCounts(int N, int M, int commSize);
unsigned int * initDispls(int commSize, unsigned int * counts);

#endif //LAB1_NCGMP_H
