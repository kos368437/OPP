#include "Matrix.h"
#include <malloc.h>
#include <stdio.h>

int realIndex(Matrix matrix, unsigned int row, unsigned int column) {
    return row * matrix.width + column;
}

Matrix createMatrix(unsigned int height, unsigned int width) {
    Matrix matrix;

    matrix.height = height;
    matrix.width = width;
    matrix.arr = (double *)calloc(height * width, sizeof(double));

    return matrix;
}

Matrix readMatrixFromFile(const char * filename) {
    Matrix matrix;
    FILE * fin;
    fin = fopen(filename, "r");

    matrix.height = 0;
    matrix.width = 0;

    if (fin == NULL) {
        matrix.arr = NULL;
        return matrix;
    }

    fscanf(fin, "%d %d", &matrix.height, &matrix.width);

    matrix = createMatrix(matrix.height, matrix.width);
    for (int i = 0; i < matrix.height * matrix.width; ++i) {
        fscanf(fin, "%lf", &matrix.arr[i] );
    }
    fclose(fin);

    return matrix;
}

void multiplyMatrix(Matrix dest,  Matrix firstOp,  Matrix secondOp) {
    for (int i = 0; i < dest.height; ++i) {
        for (int j = 0; j < dest.width; ++j) {
            dest.arr[realIndex(dest, i, j)] = 0;
            for (int k = 0; k < firstOp.width; ++k) {
                dest.arr[realIndex(dest, i, j)] +=
                        firstOp.arr[realIndex(firstOp, i, k)] *
                        secondOp.arr[realIndex(secondOp, k, j)];
            }
        }
    }
}

void copyMatrix(Matrix dest,  Matrix source) {
    if (dest.height * dest.width != source.height * source.width) {
        free(dest.arr);
        dest.arr = (double *) malloc(source.height * source.width * sizeof(double));
    }
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

double scalarMultiplicationOfVectors(Matrix first, Matrix second) {
    double result = 0;
    for (int i = 0; i < first.height; ++i) {
        result += first.arr[i] * second.arr[i];
    }
    return result;
}

int writeMatrixToFile(Matrix matrix, const char *filename) {
    FILE * fout = fopen(filename, "w");
    if (fout == NULL) return -1;

    fprintf(fout, "%d %d\n", matrix.height, matrix.width);
    for (int i = 0; i < matrix.height; i++) {
        for (int j = 0; j < matrix.width; j++) {
            fprintf(fout, "%lf ", matrix.arr[realIndex(matrix, i, j)]);
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
    return 0;
}

void transposeMatrix(Matrix * out, Matrix in) {
    (*out) = createMatrix(in.width, in.height);

    for (int i = 0; i < in.height; i++) {
        for (int j = 0; j < in.width; j++) {
            out->arr[realIndex(*out, j, i)] = in.arr[realIndex(in, i, j)];
        }
    }
}

double get(Matrix matrix, unsigned int row, unsigned int col) {
    return matrix.arr[realIndex(matrix, row, col)];
}

void set(Matrix matrix, unsigned int row, unsigned int col, double value) {
    matrix.arr[realIndex(matrix, row, col)] = value;
}