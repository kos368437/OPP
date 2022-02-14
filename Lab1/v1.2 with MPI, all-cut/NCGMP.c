#include "NCGMP.h"
#include "../Matrix.h"
#include <malloc.h>
#include <mpi.h>
#include <stdio.h>

void solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps, int * counter) {
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

    int notDone = 1;

    rn = createMatrix(b.height, 1);

    tmpVector = createMatrix(b.height, 1);

    bXb = parallelMultiplyVector(b, b);

    localx = createMatrix(b.height, 1);
    fillMatrix(localx, 0); // x.0 = 0

    /*  What algorithm requires:

        multiplyMatrix(tmpVector, A, x);
        sumMatrix(r, 1, b, -1, tmpVector); // r.0 = b - A * x.0
    */

    copyMatrix(rn, b); // That's the same with the requirement when x.0 == 0

    copyMatrix(z, rn); // z0 = r0

    rnXrn = parallelMultiplyVector(rn, rn);

    *counter = 0;

    do {
        parallelMultiplyMatrix(tmpVector, A, z); // tmpVector = A * z.n

        tmpScalar = parallelMultiplyVector(tmpVector, z); // (A * z.n, z.n)
        alpha = rnXrn / tmpScalar; // alpha.n+1 = (r.n, r.n) / (A * z.n, z.n)

        sumMatrix(localx, 1, localx, alpha, z); // x.n+1 = x.n + alpha.n+1 * z.n

        sumMatrix(rn, 1, rn, -alpha, tmpVector); // r.n+1 = rn - alpha.n+1 * A * z.n
        rn1Xrn1 = parallelMultiplyVector(rn, rn);

        beta = rn1Xrn1 / rnXrn; // beta.n+1 = (r.n+1, r.n+1) / (r.n, r.n)
        rnXrn = rn1Xrn1;

        sumMatrix(z, 1, rn, beta, z); // z.n+1 = r.n+1 + beta.n+1 * z.n

        ++(*counter);

        notDone = (rnXrn > bXb * (eps * eps));
        
    } while (notDone);

    int * recvcounts = initCounts(A.height, 1, size);
    int * displs = initDispls(size, recvcounts);

    MPI_Gatherv(localx.arr, localx.height, MPI_DOUBLE, x.arr, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    deleteMatrix(localx);

    if (rank == 0) {

        deleteMatrix(tmpVector);
        deleteMatrix(rn);
        deleteMatrix(z);
    }

}

void parallelMultiplyMatrix(Matrix dest, Matrix first, Matrix second) {

    Matrix fullResult;

    Matrix localResult = createMatrix(first.height, second.width);
    multiplyMatrix(localResult, first, second);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        fullResult = createMatrix(first.height, second.width);
    }

    MPI_Reduce(localResult.arr, fullResult.arr, localResult.height, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    int * sendcounts = initCounts(first.height, second.width, size);
    int * displs = initDispls(size, sendcounts);

    MPI_Scatterv(fullResult.arr, sendcounts, displs, MPI_DOUBLE, dest.arr, second.height, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

double parallelMultiplyVector(Matrix first, Matrix second) {
    double fullResult = 0;
    double localResult = multiplyVector(first, second);

    MPI_Allreduce(&localResult, &fullResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return fullResult;
}

void initiate(Matrix * A, Matrix * b, char * inputA, char * inputb, int commRank, int commSize) {

    Matrix AFull;
    Matrix ATemp;
    Matrix bFull;

    int * sendcounts;
    int * displs;
    int N = 0;

    if (commRank == 0) {
        AFull = readMatrixFromFile(inputA);
        bFull = readMatrixFromFile(inputb);
        N = AFull.height;
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    sendcounts = initCounts(N, N, commSize);

    displs = initDispls(commSize, sendcounts);

    (*A) = createMatrix(N, sendcounts[commRank]/N);
    ATemp = createMatrix(sendcounts[commRank]/N, N);

    (*b) = createMatrix(sendcounts[commRank]/N, 1);

    MPI_Scatterv(AFull.arr, sendcounts, displs, MPI_DOUBLE, ATemp.arr, sendcounts[commRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < sendcounts[commRank] / N; ++i) {
        for (int j = 0; j < N; ++j) {
            (*A).arr[realIndex(*A, j, i)] = ATemp.arr[realIndex(ATemp, i, j)];
        }
    }


    free(sendcounts);
    free(displs);
    sendcounts = initCounts(N, 1, commSize);
    displs = initDispls(commSize, sendcounts);

    MPI_Scatterv(bFull.arr, sendcounts, displs, MPI_DOUBLE, (*b).arr, sendcounts[commRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(sendcounts);
    free(displs);
}

int * initCounts(int N, int M, int commSize) {
    int * counts;
    counts = (int*)calloc(commSize, sizeof(int));
    for (int i = 0; i < commSize; ++i) {
        counts[i] = N / commSize;
        if (i < N % commSize) {
            counts[i] += 1;
        }
        counts[i] *= M;
    }
    return counts;
}
int * initDispls(int commSize, int * counts) {
    int * displs;
    displs = (int*)calloc(commSize, sizeof(int));
    displs[0] = 0;
    for (int i = 1; i < commSize; ++i) {
        displs[i] = displs[i-1] + counts[i-1];
    }
    return displs;
}