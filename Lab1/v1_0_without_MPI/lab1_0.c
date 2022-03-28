#include <stdio.h>
#include "../Matrix.h"
#include <time.h>
#include "NCGM.h"

int main() {
    const double eps = 1e-5;
    Matrix A = readMatrixFromFile("/home/chaos/Programming/OPP/Lab1/v1_0_without_MPI/InputA.txt");
    Matrix b = readMatrixFromFile("/home/chaos/Programming/OPP/Lab1/v1_0_without_MPI/Inputb.txt");

    int iterationCounter;
    Matrix x = createMatrix(b.height, 1);

    struct timespec tm_start, tm_finish;
    double best_time = 0;
    double time;

    for (int i = 0; i < 5; ++i) {
        clock_gettime(CLOCK_REALTIME, &tm_start);

        solveSystemUsingNCGM(A, x, b, eps, &iterationCounter);

        clock_gettime(CLOCK_REALTIME, &tm_finish);

        time = tm_finish.tv_sec - tm_start.tv_sec + 1e-9 * (tm_finish.tv_nsec - tm_start.tv_nsec);

        if (best_time == 0 || time < best_time) {
            best_time = time;
        }
    }

    printf("Iterations count: %d\n", iterationCounter);
    printf("Time elapsed: %lf sec.\n", best_time);

    writeMatrixToFile(x, "output.txt");

    deleteMatrix(A);
    deleteMatrix(b);
    deleteMatrix(x);
    return 0;
}
