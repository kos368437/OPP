#include "NCGMP.h"
#include "../Matrix.h"
#include <malloc.h>
#include <mpi.h>
#include <stdio.h>

void parallelMultiplyMatrixOnVector(Matrix dest, Matrix matrix, Matrix vector, unsigned int *sendcounts,
                                    unsigned int *displs);
double parallelScalarMultiplication(Matrix first, Matrix second);

int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Matrix z = createMatrix(b.height, 1);

    Matrix tmpVector;
    Matrix rn;
    Matrix localx;

    double tmpScalar;
    double alpha;
    double beta;
    double rn1Xrn1;
    double rnXrn;
    double bXb;

    int counter = 0;

    int notDone = 1;

    unsigned int * vectorSendcounts = initCounts(A.width, 1, size);
    unsigned int * vectorDispls = initDispls(size, vectorSendcounts);

    rn = createMatrix(b.height, 1);
    tmpVector = createMatrix(b.height, 1);

    bXb = parallelScalarMultiplication(b, b);

    localx = createMatrix(b.height, 1);
    fillMatrix(localx, 0); // x.0 = 0

    /*  What algorithm requires:

        multiplyMatrix(tmpVector, A, x);
        sumMatrix(r, 1, b, -1, tmpVector); // r.0 = b - A * x.0
    */

    copyMatrix(rn, b); // That's the same with the requirement when x.0 == 0

    copyMatrix(z, rn); // z0 = r0

    rnXrn = parallelScalarMultiplication(rn, rn);

    counter = 0;

    do {
        parallelMultiplyMatrixOnVector(tmpVector, A, z, vectorSendcounts, vectorDispls); // tmpVector = A * z.n

        tmpScalar = parallelScalarMultiplication(tmpVector, z); // (A * z.n, z.n)
        alpha = rnXrn / tmpScalar; // alpha.n+1 = (r.n, r.n) / (A * z.n, z.n)

        sumMatrix(localx, 1, localx, alpha, z); // x.n+1 = x.n + alpha.n+1 * z.n

        sumMatrix(rn, 1, rn, -alpha, tmpVector); // r.n+1 = rn - alpha.n+1 * A * z.n
        rn1Xrn1 = parallelScalarMultiplication(rn, rn);

        beta = rn1Xrn1 / rnXrn; // beta.n+1 = (r.n+1, r.n+1) / (r.n, r.n)
        rnXrn = rn1Xrn1;

        sumMatrix(z, 1, rn, beta, z); // z.n+1 = r.n+1 + beta.n+1 * z.n

        counter++;

        notDone = (rnXrn > bXb * (eps * eps));
        
    } while (notDone);

    free(vectorSendcounts);
    free(vectorDispls);

    int * recvcounts = initCounts(A.width, 1, size);
    int * displs = initDispls(size, recvcounts);

    MPI_Gatherv(localx.arr, recvcounts[rank], MPI_DOUBLE, x.arr, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(recvcounts);
    free(displs);

    deleteMatrix(localx);
    deleteMatrix(tmpVector);
    deleteMatrix(rn);
    deleteMatrix(z);

    return counter;
}

void parallelMultiplyMatrixOnVector(Matrix dest, Matrix matrix, Matrix vector, unsigned int *sendcounts,
                                    unsigned int *displs) {
    Matrix fullResultVector;

    vector.width = vector.height;
    vector.height = 1;

    Matrix localResultVector = createMatrix(vector.height, matrix.width);

    multiplyMatrix(localResultVector, vector, matrix);

    localResultVector.height = localResultVector.width;
    localResultVector.width = 1;

    vector.height = vector.width;
    vector.width = 1;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        fullResultVector = createMatrix(matrix.width, 1);
    }

    MPI_Reduce(localResultVector.arr, fullResultVector.arr, localResultVector.height, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Scatterv(fullResultVector.arr, sendcounts, displs, MPI_DOUBLE, dest.arr, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

double parallelScalarMultiplication(Matrix first, Matrix second) {
    double fullResult = 0;
    double localResult = scalarMultiplicationOfVectors(first, second);

    MPI_Allreduce(&localResult, &fullResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return fullResult;
}

void initiate(Matrix *A, Matrix *b, char *inputA, char *inputb, int commRank, int commSize, unsigned int *counts,
              unsigned int *displs) {
    Matrix AFull;
    Matrix ATemp;
    Matrix bFull;

    int N = 0;

    if (commRank == 0) {
        AFull = readMatrixFromFile(inputA);
        bFull = readMatrixFromFile(inputb);
        N = AFull.height;
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    counts = initCounts(N, N, commSize);
    displs = initDispls(commSize, counts);

    (*A) = createMatrix(counts[commRank]/N, N);

    (*b) = createMatrix(counts[commRank]/N, 1);

    MPI_Scatterv(AFull.arr, counts, displs, MPI_DOUBLE, (*A).arr, counts[commRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    unsigned int * vectorSendCounts;
    unsigned int * vectorDispls;

    vectorSendCounts = initCounts(N, 1, commSize);
    vectorDispls = initDispls(commSize, vectorSendCounts);

    MPI_Scatterv(bFull.arr, vectorSendCounts, vectorDispls, MPI_DOUBLE, (*b).arr, vectorSendCounts[commRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(vectorSendCounts);
    free(vectorDispls);
}

unsigned int * initCounts(int N, int M, int commSize) {
    unsigned int * counts;
    counts = (unsigned int*)malloc(commSize * sizeof(unsigned int));
    for (int i = 0; i < commSize; ++i) {
        counts[i] = N / commSize;
        if (i < N % commSize) {
            counts[i] += 1;
        }
        counts[i] *= M;
    }
    return counts;
}
unsigned int * initDispls(int commSize, unsigned int * counts) {
    unsigned int * displs;
    displs = (unsigned int*)malloc(commSize * sizeof(unsigned int));
    displs[0] = 0;
    for (int i = 1; i < commSize; i++) {
        displs[i] = displs[i-1] + counts[i-1];
    }
    return displs;
}