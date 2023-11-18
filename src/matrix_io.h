#pragma once

#define MATRIX_OUTPUT 1

typedef int32_t INT_T;

#define SPMATRIX_TYPE_CSC        (0)
#define SPMATRIX_TYPE_CSR        (1)
#define SPMATRIX_TYPE_INDEX0     (0<<4)
#define SPMATRIX_TYPE_INDEX1     (1<<4)
#define SPMATRIX_TYPE_ASYMMETRIC (0<<8)
#define SPMATRIX_TYPE_SYMMETRIC  (1<<8)
#define SPMATRIX_TYPE_DISTRIBUTE (1<<12)

int save_matrix_csr(char *filename, INT_T neq, INT_T ndim, INT_T* pointers, INT_T* indice, double* value, INT_T flags);
int save_vector(char *filename, INT_T size, double* v);
