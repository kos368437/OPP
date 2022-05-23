#include <stdio.h>
#include "Matrix.h"
#include <time.h>
#include "NCGM.h"

#include <omp.h>
#include <stdlib.h>


int main(int argc, char * argv[]) {
    const double eps = 1e-5;
    const unsigned int N = 300;

    omp_set_num_threads(atoi(argv[1]));

    Matrix A = getRandomSymmetricMatrix(N);
    Matrix b = getStandardResolvableVector(N);

    int iterationCounter;
    Matrix x = createMatrix(b.height, 1);

    struct timespec tm_start, tm_finish;
    double best_time = 0;
    double time;

    clock_gettime(CLOCK_REALTIME, &tm_start);

    iterationCounter = solveSystemUsingNCGM(A, x, b, eps);

    clock_gettime(CLOCK_REALTIME, &tm_finish);

    time = tm_finish.tv_sec - tm_start.tv_sec + 1e-9 * (tm_finish.tv_nsec - tm_start.tv_nsec);

    if (best_time == 0 || time < best_time) {
        best_time = time;
    }

    printf("Iterations count: %d\n", iterationCounter);
    printf("Time elapsed: %lf sec.\n", best_time);

    //writeMatrixToConsole(x);

    deleteMatrix(A);
    deleteMatrix(b);
    deleteMatrix(x);
    return 0;
}
