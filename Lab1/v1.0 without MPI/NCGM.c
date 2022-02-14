#include "NCGM.h"
#include "../Matrix.h"

void solveSystemUsingNCGM(Matrix A, Matrix x, Matrix b, double eps, int * counter) {
    Matrix rn = createMatrix(b.height, 1);
    Matrix z = createMatrix(b.height, 1);

    double alpha;
    double beta;

    Matrix tmpVector = createMatrix(b.height, 1);
    double tmpScalar;

    double bXb;
    bXb = multiplyVector(b, b);

    fillMatrix(x, 0); // x.0 = 0

    /*  What algorithm requires:

        multiplyMatrix(tmpVector, A, x);
        sumMatrix(r, 1, b, -1, tmpVector); // r.0 = b - A * x.0
    */

    copyMatrix(rn, b); // That's the same with the requirement when x.0 == 0

    copyMatrix(z, rn); // z0 = r0

    double rnXrn;
    rnXrn = multiplyVector(rn, rn);

    double rn1Xrn1; // (r.n+1, r.n+1)

    *counter = 0;

    do {
        multiplyMatrix(tmpVector, A, z); // tmpVector = A * z.n

        tmpScalar = multiplyVector(tmpVector, z); // (A * z.n, z.n)
        alpha = rnXrn / tmpScalar; // alpha.n+1 = (r.n, r.n) / (A * z.n, z.n)

        sumMatrix(x, 1, x, alpha, z); // x.n+1 = x.n + alpha.n+1 * z.n

        sumMatrix(rn, 1, rn, -alpha, tmpVector); // r.n+1 = rn - alpha.n+1 * A * z.n
        rn1Xrn1 = multiplyVector(rn, rn);

        beta = rn1Xrn1 / rnXrn; // beta.n+1 = (r.n+1, r.n+1) / (r.n, r.n)
        rnXrn = rn1Xrn1;

        sumMatrix(z, 1, rn, beta, z); // z.n+1 = r.n+1 + beta.n+1 * z.n

        ++(*counter);
    } while (rnXrn > bXb * (eps * eps));

    deleteMatrix(tmpVector);
    deleteMatrix(rn);
    deleteMatrix(z);

}

