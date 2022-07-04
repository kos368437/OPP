#include "ParallelMultiplication.h"
#include <mpi.h>
#include <stdlib.h>

void broadcastMatrix(Matrix matrix, MPI_Comm comm, int root_rank);
int * initResultAssembleCounts(int * f_sendcounts, unsigned int N2, unsigned int N3, int sizey);
void assembleResult(Matrix result, Matrix local_result, int *row_assemble_counts, int *assemble_result_counts,
                    int *f_sendcounts, int ranky, int rankx, int sizey, int sizex, unsigned int N2, unsigned int N3,
                    MPI_Comm RowComm, MPI_Comm ColComm);
int * initRowAssembleCounts(int * f_sendcounts, int * s_sendcounts, unsigned int N2, int ranky, int sizex);

void scatterMatrix(Matrix source, Matrix * dest, int * sendcounts, int * displs, MPI_Comm comm, int root_rank, unsigned int width);

Matrix createLocalOp(int * counts, unsigned int width, int line_rank);

int parallelMatrixMultiplication(Matrix result, Matrix f_op, Matrix s_op_trans, int sizey, int sizex, int *dims,
                             MPI_Comm GridComm, MPI_Comm RowComm, MPI_Comm ColComm, unsigned int N1,
                             unsigned int N2, unsigned int N3) {

    int periods[2]={0,0}, coords[2];
    int ranky, rankx;

    Matrix f_op_local;
    Matrix s_op_trans_local;

    MPI_Cart_get(GridComm, 2, dims, periods, coords);

    ranky = coords[0];
    rankx = coords[1];

    int rank, col_rank, row_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_rank(ColComm, &col_rank);
    MPI_Comm_rank(RowComm, &row_rank);

    int * f_sendcounts = initCounts(N1, N2, sizey);
    int * f_displs = initDispls(sizey, f_sendcounts);

    f_op_local = createLocalOp(f_sendcounts, N2, col_rank);

    int * s_sendcounts = initCounts(N3, N2, sizex);
    int * s_displs = initDispls(sizex, s_sendcounts);

    s_op_trans_local = createLocalOp(s_sendcounts, N2, row_rank);

    Matrix local_result = createMatrix(f_sendcounts[col_rank] / N2, s_sendcounts[row_rank] / N2);

    if (rankx == 0) {
        scatterMatrix(f_op, &f_op_local, f_sendcounts, f_displs, ColComm, 0, N2);
    }

    broadcastMatrix(f_op_local, RowComm, 0);

    if (ranky == 0) {

        scatterMatrix(s_op_trans, &s_op_trans_local, s_sendcounts, s_displs, RowComm, 0, N2);
    }
    broadcastMatrix(s_op_trans_local, ColComm, 0);

    matrixTransposedMultiplication(local_result, f_op_local, s_op_trans_local);

    int * row_assemble_counts = initRowAssembleCounts(f_sendcounts, s_sendcounts, N2, ranky, sizex);
    int * assemble_result_counts = initResultAssembleCounts(f_sendcounts, N2, N3, sizey);

    assembleResult(result, local_result, row_assemble_counts, assemble_result_counts, f_sendcounts, ranky, rankx, sizey,
                   sizex, N2, N3, RowComm, ColComm);

    deleteMatrix(local_result);
    deleteMatrix(f_op_local);
    deleteMatrix(s_op_trans_local);

    free(row_assemble_counts);
    free(assemble_result_counts);
    free(f_sendcounts);
    free(f_displs);
    free(s_sendcounts);
    free(s_displs);

    return 0;
}

void broadcastMatrix(Matrix matrix, MPI_Comm comm, int root_rank) {

    MPI_Bcast(matrix.arr, matrix.height * matrix.width, MPI_DOUBLE, root_rank, comm);
}

int * initRowAssembleCounts(int * f_sendcounts, int * s_sendcounts, unsigned int N2, int ranky, int sizex) {

    int * counts = (int *)malloc(sizeof(int) * sizex);

    for (int i = 0; i < sizex; i++) {
        counts[i] = (f_sendcounts[ranky] * s_sendcounts[i]) / (N2 * N2);
    }

    return counts;
}

int * initResultAssembleCounts(int * f_sendcounts, unsigned int N2, unsigned int N3, int sizey) {

    int * counts = (int *)malloc(sizeof(int) * sizey);

    for (int i = 0; i < sizey; i++) {
        counts[i] = N3 * f_sendcounts[i] / N2;
    }
    return counts;
}


void assembleResult(Matrix result, Matrix local_result, int *row_assemble_counts, int *assemble_result_counts,
                    int *f_sendcounts, int ranky, int rankx, int sizey, int sizex, unsigned int N2, unsigned int N3,
                    MPI_Comm RowComm, MPI_Comm ColComm) {
    Matrix row;
    Matrix transposed_row;
    Matrix local_result_transposed;
    int * row_assemble_displs = initDispls(sizex, row_assemble_counts);
    int * assemble_result_displs = initDispls(sizey, assemble_result_counts);

    if (rankx == 0) {
        row = createMatrix(N3, f_sendcounts[ranky] / N2);
    }

    transposeMatrix(&local_result_transposed, local_result);

    MPI_Gatherv(local_result_transposed.arr, row_assemble_counts[rankx], MPI_DOUBLE, row.arr, row_assemble_counts, row_assemble_displs, MPI_DOUBLE, 0, RowComm);

    if (rankx == 0) {
        transposeMatrix(&transposed_row, row);

        deleteMatrix(row);

        MPI_Gatherv(transposed_row.arr, assemble_result_counts[ranky], MPI_DOUBLE, result.arr, assemble_result_counts, assemble_result_displs, MPI_DOUBLE, 0, ColComm);

        deleteMatrix(transposed_row);
    }

    deleteMatrix(local_result_transposed);

    free(row_assemble_displs);
    free(assemble_result_displs);
}

Matrix createLocalOp(int * counts, unsigned int width, int line_rank) {

    return createMatrix(counts[line_rank] / width, width);
}

void scatterMatrix(Matrix source, Matrix * dest, int * sendcounts, int * displs, MPI_Comm comm, int root_rank, unsigned int width) {

    int rank;
    MPI_Comm_rank(comm, &rank);

    MPI_Scatterv(source.arr, sendcounts, displs, MPI_DOUBLE, (*dest).arr, sendcounts[rank], MPI_DOUBLE, root_rank, comm);
}

void randomInitiate(Matrix *A, Matrix *B_trans, Matrix *C, unsigned int N1, unsigned int N2, unsigned int N3) {

    (*A) = getRandomSymmetricMatrix(N1, N2);
    (*B_trans) = getRandomSymmetricMatrix(N3, N2);
    (*C) = createMatrix(N1, N3);
}

int * initCounts(unsigned int height, unsigned int width, int commSize) {
    unsigned int * counts;
    counts = (unsigned int*)malloc(commSize * sizeof(unsigned int));
    for (int i = 0; i < commSize; ++i) {
        counts[i] = height / commSize;
        if (i < height % commSize) {
            counts[i] += 1;
        }
        counts[i] *= width;
    }

    return counts;
}

unsigned int * initDispls(int commSize, unsigned int * counts) {

    unsigned int * displs;

    displs = (unsigned int*)malloc(commSize * sizeof(unsigned int));
    displs[0] = 0;
    for (int i = 1; i < commSize; i++) {
        displs[i] = displs[i-1] + counts[i-1];
    }

    return displs;
}

Matrix getRandomSymmetricMatrix(unsigned int N, unsigned int M) {

    Matrix matrix = createMatrix(N, M);

    for (int i = 0; i < matrix.height; i++) {
        for (int j = 0; j < matrix.width; j++) {
            srandom(i + j);
            set(matrix, i, j, random() % N + 1);
        }
    }

    return matrix;
}