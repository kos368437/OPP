#include "NCGMP.h"
#include "../Matrix.h"
#include <malloc.h>
#include <mpi.h>

void solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps, int * counter) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

    if (rank == 0) {
        rn = createMatrix(b.height, 1);

        tmpVector = createMatrix(b.height, 1);

        bXb = multiplyVector(b, b);

        fillMatrix(x, 0); // x.0 = 0

        /*  What algorithm requires:

            multiplyMatrix(tmpVector, A, x);
            sumMatrix(r, 1, b, -1, tmpVector); // r.0 = b - A * x.0
        */

        copyMatrix(rn, b); // That's the same with the requirement when x.0 == 0

        copyMatrix(z, rn); // z0 = r0

        rnXrn = multiplyVector(rn, rn);

        *counter = 0;
    }

    do {
	
        MPI_Bcast(z.arr, z.height, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        parallelMultiplyMatrix(tmpVector, A, z, tempResult); // tmpVector = A * z.n

        if (rank == 0) {
            tmpScalar = multiplyVector(tmpVector, z); // (A * z.n, z.n)
            alpha = rnXrn / tmpScalar; // alpha.n+1 = (r.n, r.n) / (A * z.n, z.n)

            sumMatrix(x, 1, x, alpha, z); // x.n+1 = x.n + alpha.n+1 * z.n

            sumMatrix(rn, 1, rn, -alpha, tmpVector); // r.n+1 = rn - alpha.n+1 * A * z.n
            rn1Xrn1 = multiplyVector(rn, rn);

            beta = rn1Xrn1 / rnXrn; // beta.n+1 = (r.n+1, r.n+1) / (r.n, r.n)
            rnXrn = rn1Xrn1;

            sumMatrix(z, 1, rn, beta, z); // z.n+1 = r.n+1 + beta.n+1 * z.n

            ++(*counter);

            notDone = (rnXrn > bXb * (eps * eps));

        }
        MPI_Bcast(&notDone, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
    } while (notDone);
    
    if (rank == 0) {

        deleteMatrix(tmpVector);
        deleteMatrix(rn);
        deleteMatrix(z);
    }

}

void parallelMultiplyMatrix(Matrix dest, Matrix first, Matrix second, Matrix tempResult) {

    multiplyMatrix(tempResult, first, second);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int * recvcounts = (int*)calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        recvcounts[i] = first.width / size;
        if (i < first.width % size) {
            recvcounts[i] += 1;
        }
    }
    
    int * displs = (int*)calloc(size, sizeof(int));
    displs[0] = 0;
    for (int i = 1; i < size; ++i) {
        displs[i] = displs[i-1] + recvcounts[i-1];
    }
    MPI_Gatherv(tempResult.arr, tempResult.height, MPI_DOUBLE, dest.arr, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

