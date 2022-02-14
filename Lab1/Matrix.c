#include "Matrix.h"
#include <malloc.h>
#include <stdio.h>

Matrix createMatrix(int height, int width) {
    Matrix matrix;

    matrix.height = height;
    matrix.width = width;
    matrix.arr = (double *)calloc(height * width, sizeof(double));

    return matrix;
}

Matrix readMatrixFromFile(char * filename) {
    FILE * fin;
    fin = fopen(filename, "r");

    Matrix matrix;
    matrix.height = 0;
    matrix.width = 0;

    fscanf(fin, "%d %d", &matrix.height, &matrix.width);

    matrix = createMatrix(matrix.height, matrix.width);
    for (int i = 0; i < matrix.height * matrix.width; ++i) {
        fscanf(fin, "%lf", &matrix.arr[i] );
    }
    fclose(fin);

    return matrix;
}

void multiplyMatrix( Matrix dest,  Matrix firstOp,  Matrix secondOp) {
    for (int i = 0; i < dest.height; ++i) {
        for (int j = 0; j < dest.width; ++j) {
            dest.arr[i * dest.width + j] = 0;
            for (int k = 0; k < firstOp.width; ++k) {
                dest.arr[realIndex(dest, i, j)] +=
                        firstOp.arr[realIndex(firstOp, i, k)] *
                        secondOp.arr[realIndex(secondOp, k, j)];
            }
        }
    }
}

int realIndex(Matrix matrix, int row, int column) {
    return row * matrix.width + column;
}

void copyMatrix( Matrix dest,  Matrix source) {
    dest.height = source.height;
    dest.width = source.width;

    for (int i = 0; i < dest.height * dest.width; ++i) {
        dest.arr[i] = source.arr[i];
    }
}

void deleteMatrix(Matrix del) {
    free(del.arr);
}

void fillMatrix(Matrix matrix, double fillEl) {
    for(int i = 0; i < matrix.height * matrix.width; ++i) {
        matrix.arr[i] = fillEl;
    }
}

void sumMatrix(Matrix result, double a, Matrix firstOp, double b, Matrix secondOp) { // a*firstOp + b*secondOp
    for (int i = 0; i < firstOp.height * firstOp.width; ++i) {
        result.arr[i] = a * firstOp.arr[i] + b * secondOp.arr[i];
    }
}

void multiplyMatrixOnScalar(Matrix result, Matrix matrix, double scal) {
    for (int i = 0; i < matrix.height * matrix.width; ++i) {
        result.arr[i] = matrix.arr[i] * scal;
    }
}

double multiplyVector(Matrix first, Matrix second) {
    double result = 0;
    for (int i = 0; i < first.height; ++i) {
        result += first.arr[i] * second.arr[i];
    }
    return result;
}

void writeMatrixToFile(Matrix matrix, char * filename) {
    FILE * fout = fopen("output.txt", "w+");

    fprintf(fout, "%d %d\n", matrix.height, matrix.width);
    for (int i = 0; i < matrix.height; ++i) {
        for (int j = 0; j < matrix.width; ++j) {
            fprintf(fout, "%lf ", matrix.arr[realIndex(matrix, i, j)]);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
}
