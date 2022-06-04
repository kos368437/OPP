#ifndef LAB1_MATRIX_H
#define LAB1_MATRIX_H

typedef struct {
    double * arr;
    unsigned int height;
    unsigned int width;
} Matrix;

Matrix createMatrix(unsigned int height, unsigned int width);
Matrix readMatrixFromFile(const char * filename);

void multiplyMatrix(Matrix dest,  Matrix firstOpMatrix,  Matrix secondOpMatrix);
void copyMatrix( Matrix dest,  Matrix source);
void deleteMatrix(Matrix del);
void sumMatrix(Matrix result, double a, Matrix firstOp, double b, Matrix secondOp); //result = a * first + b * second
void multiplyMatrixOnScalar(Matrix result, Matrix matrix, double scal);
void fillMatrix(Matrix matrix, double fillEl);

double get(Matrix matrix, unsigned int row, unsigned int col);
void set(Matrix matrix, unsigned int row, unsigned int col, double value);

double scalarMultiplicationOfVectors(Matrix first, Matrix second);

int writeMatrixToFile(Matrix matrix, const char *filename);
void writeMatrixToConsole(Matrix matrix);
void matrixTransposedMultiplication(Matrix dest,  Matrix firstOp,  Matrix secondOpTransposed);

void transposeMatrix(Matrix * out, Matrix in);


#endif //LAB1_MATRIX_H
