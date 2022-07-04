#include <mpi.h>
#include "Matrix.h"


int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps);

void initiate(Matrix *A, Matrix *B_trans, Matrix *C, unsigned int N1, unsigned int N2, unsigned int N3);
void randomInitiate(Matrix *A, Matrix *B_trans, Matrix *C, unsigned int N1, unsigned int N2, unsigned int N3);

int * initCounts(unsigned int height, unsigned int width, int commSize);
unsigned int * initDispls(int commSize, unsigned int * counts);

int
parallelMatrixMultiplication(Matrix result, Matrix f_op, Matrix s_op_trans, int sizey, int sizex, int *dims,
                             MPI_Comm GridComm, MPI_Comm RowComm, MPI_Comm ColComm, unsigned int N1,
                             unsigned int N2, unsigned int N3);

//Matrix getStandardSymmetricResolvableMatrix(unsigned int N, int commRank, int commSize);
//Matrix getStandardResolvableVector(unsigned int N, int commRank, int commSize);
//Matrix getStandardRandomResolvableVector(unsigned int N, int commRank, int commSize, Matrix A);
//Matrix getRandomSymmetricMatrix(unsigned int N, int commRank, int commSize);

Matrix getRandomSymmetricMatrix(unsigned int N, unsigned int M);
Matrix getMatrix(unsigned int N, unsigned int M);


