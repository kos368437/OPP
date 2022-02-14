#ifndef LAB1_MATRIX_H
#define LAB1_MATRIX_H

typedef struct {
    double * arr;
    int height;
    int width;
} Matrix;

Matrix createMatrix(int height, int width);
Matrix readMatrixFromFile(char * filename);

void multiplyMatrix( Matrix dest,  Matrix firstOpMatrix,  Matrix secondOpMatrix);
void copyMatrix( Matrix dest,  Matrix source);
void deleteMatrix(Matrix del);
void sumMatrix(Matrix result, double a, Matrix firstOp, double b, Matrix secondOp); //result = a * first + b * second
void multiplyMatrixOnScalar(Matrix result, Matrix matrix, double scal);
void fillMatrix(Matrix matrix, double fillEl);

int realIndex(Matrix matrix, int row, int column);
double multiplyVector(Matrix first, Matrix second);

void writeMatrixToFile(Matrix matrix, char * filename);


#endif //LAB1_MATRIX_H
