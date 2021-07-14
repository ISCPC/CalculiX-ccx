/*
 * vesolver.h
 *
 *  Created on: Jun 11, 2021
 *      Author: uno
 */

#ifndef VESOLVER_H_
#define VESOLVER_H_
#include <stdint.h>

typedef int32_t INT_T;

/*
 * From Matrix.h
 */
typedef struct matrix_desc {
	INT_T neq;
	INT_T nnz;
	uint32_t flags;
	INT_T* pointers;
	INT_T* indice;
	double* values;
	double* bptr;
	double* xptr;
	int64_t handle;
} matrix_desc_t;

#define MATRIX_TYPE_INDEX0     (0)
#define MATRIX_TYPE_INDEX1     (1)
#define MATRIX_TYPE_ASYMMETRIC (0<<4)
#define MATRIX_TYPE_SYMMETRIC  (1<<4)
#define MATRIX_TYPE_LOWER      (1<<5)
#define MATRIX_TYPE_UPPER      (0)
#define MATRIX_TYPE_UNIT       (1<<6)
#define MATRIX_TYPE_NON_UNIT   (0)
#define MATRIX_TYPE_CSC        (0<<8)
#define MATRIX_TYPE_CSR        (1<<8)

#if 0
#define MATRIX_INDEX_TYPE(A)    (((A)->flags)&0xf)
#define MATRIX_IS_SYMMETRIC(A)  (((A)->flags)&MATRIX_TYPE_SYMMETRIC)
#define MATRIX_IS_LOWER(A)      (((A)->flags)&MATRIX_TYPE_LOWER)
#define MATRIX_IS_UNIT(A)       (((A)->flags)&MATRIX_TYPE_UNIT)
#endif

/*
 * From LinearSolver.h
 */
#define VESOLVER_OPTION_SOLVER       1

#define VESOLVER_ITER_CG             1
#define VESOLVER_ITER_BICGSTAB2      2
#define VESOLVER_DIRECT_HS           3
#define VESOLVER_DIRECT_PARDISO      4

/*
 * VE solver API
 */
typedef int vesolver_handle_t;

#ifdef __cplusplus
extern "C" {
#endif

int vesolver_init();
void vesolver_finalize();
vesolver_handle_t vesolver_activate();
int vesolver_deactivate(vesolver_handle_t hdl);
int vesolver_set_option(vesolver_handle_t hdl, int id, int value);
#if 0
int vesolver_set_matrix_csr(vesolver_handle_t hdl, INT_T neq,
	INT_T* pointers, INT_T* indice, double* value, uint32_t flags);
int vesolver_set_matrix_csc(vesolver_handle_t hdl, INT_T neq,
	INT_T* pointers, INT_T* indice, double* value, uint32_t flags);
#else
matrix_desc_t* vesolver_alloc_matrix(vesolver_handle_t hdl, INT_T neq,
	INT_T nnz, uint32_t flags);
int vesolver_set_matrix(vesolver_handle_t hdl, matrix_desc_t *desc);
#endif
int vesolver_free_matrix(vesolver_handle_t hdl);
int vesolver_solve_sync(vesolver_handle_t hdl, double* b, double* x, double res);
//int vesolver_solve_async();
//int vesolver_solve_wait();

#ifdef __cplusplus
}
#endif

#endif /* VESOLVER_H_ */
