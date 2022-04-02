#include "NCGM.h"
#include "../Matrix.h"
#include <math.h>
#include <stdio.h>

int solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps) {
    Matrix rn = createMatrix(b.height, 1);
    Matrix z = createMatrix(b.height, 1);

    double alpha;
    double beta;

    Matrix tmpVector = createMatrix(b.height, 1);
    double tmpScalar;

    double bXb;
    bXb = scalarMultiplicationOfVectors(b, b);

    fillMatrix(x, 0); // x.0 = 0

    /*  What algorithm requires:
        multiplyMatrix(tmpVector, A, x);
        sumMatrix(r, 1, b, -1, tmpVector); // r.0 = b - A * x.0
    */

    copyMatrix(rn, b); // That's the same with the requirement when x.0 == 0
    copyMatrix(z, rn); // z0 = r0

    double rnXrn;
    rnXrn = scalarMultiplicationOfVectors(rn, rn);

    double rn1Xrn1; // (r.n+1, r.n+1)

    int counter = 0;

    do {
        multiplyMatrix(tmpVector, A, z); // tmpVector = A * z.n

        tmpScalar = scalarMultiplicationOfVectors(tmpVector, z); // (A * z.n, z.n)
        alpha = rnXrn / tmpScalar; // alpha.n+1 = (r.n, r.n) / (A * z.n, z.n)

        sumMatrix(x, 1, x, alpha, z); // x.n+1 = x.n + alpha.n+1 * z.n

        sumMatrix(rn, 1, rn, -alpha, tmpVector); // r.n+1 = rn - alpha.n+1 * A * z.n
        rn1Xrn1 = scalarMultiplicationOfVectors(rn, rn);

        beta = rn1Xrn1 / rnXrn; // beta.n+1 = (r.n+1, r.n+1) / (r.n, r.n)
        rnXrn = rn1Xrn1;

        sumMatrix(z, 1, rn, beta, z); // z.n+1 = r.n+1 + beta.n+1 * z.n

        counter++;
    } while (rnXrn > bXb * (eps * eps));

    deleteMatrix(tmpVector);
    deleteMatrix(rn);
    deleteMatrix(z);

    return counter;
}

Matrix getStandardSymmetricResolvableMatrix(unsigned int N) {
    Matrix matrix = createMatrix(N, N);

    for (int i = 0; i < matrix.height; i++) {
        for (int j = 0; j < matrix.width; j++) {
            if (i == j) {
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

Matrix getStandardRandomResolvableVector(unsigned int N, Matrix A) {
    Matrix vector = createMatrix(N, 1);
    Matrix tempResult = createMatrix(N, 1);

    for (int i = 0; i < tempResult.height; i++) {
        set(tempResult, i, 0, sin((2 * M_PI * i) / N));
    }

    multiplyMatrix(vector, A, tempResult);

    deleteMatrix(tempResult);

    return vector;
}
