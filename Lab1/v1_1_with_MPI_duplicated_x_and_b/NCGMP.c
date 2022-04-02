#include "NCGMP.h"
#include "../Matrix.h"
#include <malloc.h>
#include <mpi.h>
#include <math.h>

void parallelMultiplyMatrix(Matrix dest, Matrix first, Matrix second, Matrix tempResult, unsigned int *recvcounts,
                            unsigned int *displs);

int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps) {
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Matrix tempResult = createMatrix(A.height, 1);
    Matrix z = createMatrix(A.width, 1);

    Matrix tmpVector;
    Matrix rn;

    double tmpScalar;
    double alpha;
    double beta;
    double rn1Xrn1;
    double rnXrn;
    double bXb;

    int notDone = 1;
    int counter = 0;

    unsigned int * vectorCounts = initCounts(size, A.width, 1);
    unsigned int * vectorDispls = initDispls(size, vectorCounts);

    if (rank == 0) {
        rn = createMatrix(b.height, 1);

        tmpVector = createMatrix(b.height, 1);

        bXb = scalarMultiplicationOfVectors(b, b);

        fillMatrix(x, 0); // x.0 = 0

        /*  What algorithm requires:

            multiplyMatrix(tmpVector, A, x);
            sumMatrix(r, 1, b, -1, tmpVector); // r.0 = b - A * x.0
        */

        copyMatrix(rn, b); // That's the same with the requirement when x.0 == 0

        copyMatrix(z, rn); // z0 = r0

        rnXrn = scalarMultiplicationOfVectors(rn, rn);
   }
    counter = 0;
    do {
	    MPI_Bcast(z.arr, z.height, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        parallelMultiplyMatrix(tmpVector, A, z, tempResult, vectorCounts, vectorDispls); // tmpVector = A * z.n

        if (rank == 0) {
            tmpScalar = scalarMultiplicationOfVectors(tmpVector, z); // (A * z.n, z.n)
            alpha = rnXrn / tmpScalar; // alpha.n+1 = (r.n, r.n) / (A * z.n, z.n)

            sumMatrix(x, 1, x, alpha, z); // x.n+1 = x.n + alpha.n+1 * z.n

            sumMatrix(rn, 1, rn, -alpha, tmpVector); // r.n+1 = rn - alpha.n+1 * A * z.n
            rn1Xrn1 = scalarMultiplicationOfVectors(rn, rn);

            beta = rn1Xrn1 / rnXrn; // beta.n+1 = (r.n+1, r.n+1) / (r.n, r.n)
            rnXrn = rn1Xrn1;

            sumMatrix(z, 1, rn, beta, z); // z.n+1 = r.n+1 + beta.n+1 * z.n

            notDone = (rnXrn > bXb * (eps * eps));
        }

        counter++;

        MPI_Bcast(&notDone, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
    } while (notDone);
    
    if (rank == 0) {
        deleteMatrix(tmpVector);
        deleteMatrix(rn);
        deleteMatrix(z);
    }

    free(vectorCounts);
    free(vectorDispls);

    return counter;
}

void parallelMultiplyMatrix(Matrix dest, Matrix first, Matrix second, Matrix tempResult, unsigned int *recvcounts,
                            unsigned int *displs) {
    multiplyMatrix(tempResult, first, second);

    MPI_Gatherv(tempResult.arr, tempResult.height, MPI_DOUBLE, dest.arr, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

unsigned int * initCounts(int commSize, unsigned int height, unsigned int width) {
    unsigned int * counts = (unsigned int *)malloc(commSize * sizeof(unsigned int));
    for (int i = 0; i < commSize; ++i) {
        counts[i] = height / commSize;
        if (i < height % commSize) {
            counts[i] += 1;
        }
        counts[i] *= width;
    }
    return counts;
}

unsigned int * initDispls(int commSize, const unsigned int * counts) {
    unsigned int * displs = (unsigned int *)malloc(commSize * sizeof(unsigned int));

    displs[0] = 0;
    for (int i = 1; i < commSize; ++i) {
        displs[i] = displs[i-1] + counts[i-1];
    }

    return displs;
}

unsigned int getActualStartRowNumber(unsigned int N, int commRank, int commSize) {
    unsigned int height = N / commSize;
    unsigned int actualStartRowNumber;

    if (commRank < N % commSize) {
        height += 1;

        actualStartRowNumber = commRank * height;
    }
    else {
        actualStartRowNumber = commRank * height + (N % commSize);
    }

    return actualStartRowNumber;
}

Matrix getStandardSymmetricResolvableMatrix(unsigned int N, int commRank, int commSize) {
    unsigned int height = N / commSize;
    unsigned int actualStartRowNumber = getActualStartRowNumber(N, commRank, commSize);

    if (commRank < N % commSize) {
        height += 1;
    }

    Matrix matrix = createMatrix(height, N);

    for (int i = 0; i < matrix.height; i++) {
        for (int j = 0; j < matrix.width; j++) {
            if (actualStartRowNumber + i == j) {
                set(matrix, i, j, 2.0);
            }
            else {
                set(matrix, i, j, 1.0);
            }
        }
    }

    return matrix;
}

Matrix getStandardResolvableVector(unsigned int N) {
    Matrix vector = createMatrix(N, 1);

    for (int i = 0; i < vector.height; i++) {
        set(vector, i, 0, N + 1);
    }

    return vector;
}

Matrix getStandardRandomResolvableVector(unsigned int N, int commRank, int commSize, Matrix A) {
    Matrix tempResult = createMatrix(A.height, 1);
    Matrix  vector = createMatrix(N, 1);
    for (int i = 0; i < vector.height; i++) {
        set(vector, i, 0, sin((2 * M_PI * i) / N));
    }
    int * vectorCounts = initCounts(commSize, N, 1);
    int * vectorDispls = initDispls(commSize, vectorCounts);

    parallelMultiplyMatrix(vector, A, vector, tempResult, vectorCounts, vectorDispls);

    free(vectorCounts);
    free(vectorDispls);
    deleteMatrix(tempResult);

    return vector;
}
