/*

  pcgsolver.c	--	The preconditioned conjugate gradient solver(s)
  
  				Implemented solver/preconditioner:
  				    preconditioned cg-solver with
  				       * diagonal scaling only
  				       * incomplete Cholesky on full matrix
  				
  				Most functions are based on:
  				H.R. Schwarz: FORTRAN-Programme zur Methode 
                                             der finiten Elemente,
  				              Teubner, 1981 
 
                               The present version is based on the c-code by
                               Martin Ruecker and Ernst Rank of the Technical
                               University of Munich 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <strings.h>
#include <sys/time.h>
#include "CalculiX.h"
#ifdef SBLAS
#include <sblas.h>
#endif

#ifdef DISABLE_INLINE
#define INLINE
#else
#define INLINE inline
#endif

/* Prototyping	*/
uint64_t cgve_solve(ITG neq, ITG len, uint64_t aptr, uint64_t aind, uint64_t val, uint64_t bptr, uint64_t xptr);
static INLINE void CG (double *A, double *x, double *b, ITG neq, ITG len, ITG *ia, ITG *iz,
        double *eps, ITG *niter);
static INLINE void Scaling (double *A, double *b, ITG neq, ITG *ia, ITG *iz, double *d);
static INLINE void MatVecProduct(const double * const A, const double * const p, const ITG neq, const ITG * const ia, const ITG * const iz, double *z);
//static INLINE void MatVecProduct(double *A, double *p, ITG neq, ITG *ia, ITG *iz, double *z);
static INLINE void InnerProduct(double *a, double *b, ITG n, double *Sum);


/* **********************************************************************

   Scaling the equation system A x + b = 0

   The equation system is scaled in consideration of keeping the symmetry in
   such a way that the diagonal elements of matrix A are 1. This is performed
   by the diagonal matrix Ds with the diagonal elements d_i = 1/sqrt(a_ii).
   The given equation system Ax+b=0 is transformed into
   -1                     - -   -
   Ds A Ds  Ds x + Ds b = 0   or    A x + b = 0
   _                                             _ 
   with the scaled Matrix A= Ds A Ds  and  the scaled right hand side b= Ds b.
   The scaling factor Ds is needed later for backscaling of the solution 
   vector
   _    -1                       _ 
   x = Ds x      resp.    x = Ds x

parameter:
A       compact row oriented storage of lower left of matrix A
b       right hand side
neq     order of A, number of equations
ia      column indices
iz      row indices (diagonal elements indices)

The function corresponds to function SCALKO() in H.R. Schwarz: FORTRAN-Pro-
gramme zur Methode der finiten Elemente, p.105, Teubner, 1981

 ********************************************************************** 
 */
void Scaling (double *A, double *b, ITG neq, ITG *ia, ITG *iz, double *d)
{
    ITG			i=0, j=0, jlo=0, jup=0;

    /*  extract diagonal vector from matrix A  */
    for (i=0; i<neq; i++) {
        d[i] =  1.0/sqrt(A[iz[i]]);
    }

    /*  scale right hand side (Ax=b -> Ax+b=0: negative sign)  */
    for (i=0; i<neq; i++) {
        b[i] *= -d[i];
    }

    /*  scale matrix A  */
    for (i=0; i<neq; i++) {
        jlo = iz[i];
        jup = iz[i+1];
        for (j=jlo; j<jup; j++) {
            A[j] *= d[i]*d[ia[j]];
        }
    }
}


/* **********************************************************************

   Computation of matrix vector product z = A p

parameter:
A       compact row oriented storage of lower left of matrix A
p       vector
neq     order of A, number of equations
ia      column indices
iz      row indices (diagonal elements indices)

The function corresponds to function APZ() in H.R. Schwarz: FORTRAN-
Programme zur Methode der finiten Elemente, p.111, Teubner, 1981

 ********************************************************************** 
 */
void MatVecProduct(const double * const A, const double * const p, const ITG neq, const ITG * const ia, const ITG * const iz, double *z)
{
    ITG		i=0, j=0, jlo=0, jup=0, k=0;

    /*  matrix vector product  */
    for (i=neq-1; i>=0; i--) {
        z[i] = A[iz[i]]*p[i];
        jlo = iz[i]+1;
        jup = iz[i+1];
        for (j=jlo; j<jup; j++) {
            k = ia[j];
            z[i] += A[j]*p[k];
            z[k] += A[j]*p[i];
        }
    }
    return;
}

/*--Calculation of the inner product of two (distributed) vectors------------------	*/
/*---------------------------------------------------------------------------------	*/
void InnerProduct(double *a, double *b, ITG n, double *Sum)
{
    ITG 	i=0;

    *Sum=0.;
    for (i=0; i<n; i++)	 *Sum += a[i]*b[i];
    return;
}
/*---------------------------------------------------------------------------------	*/



/*--(parallel) conjugate gradient solver-------------------------------------------	*/
/*---------------------------------------------------------------------------------	*/
/*--parameter:                                                                   --	*/
/*--           A       compact row oriented storage of lower left of the scaled  --	*/
/*--                   matrix A                                                  --	*/
/*--           x       solution vector                                           --	*/
/*--           b       right hand side                                           --	*/
/*--           neq     order of A, number of equations                           --	*/
/*--           len     number of non zero entries in Matrix A                    --	*/
/*--           ia      column indices                                            --	*/
/*--           iz      row indices (diagonal elements indices)                   --	*/
/*--           eps     required accuracy -> residual                             --	*/
/*--           niter   maximum number of iterations -> number of iterations      --	*/
/*---------------------------------------------------------------------------------	*/
void CG (double *A, double *x, double *b, ITG neq, ITG len, ITG *ia, ITG *iz,
        double *eps, ITG *niter)
{
    ITG	i=0, k=0, ncg=0,iam;
    double ekm1=0.0,c1=0.005,qam,ram=0.,err;
    double rr=0.0, pz=0.0, qk=0.0, rro=0.0;

    double* r = (double*)calloc(sizeof(double), neq);
    double* p = (double*)calloc(sizeof(double), neq);
    double* z = (double*)calloc(sizeof(double), neq);

    /* solving the system of equations using the iterative solver */

    //printf("Solving the system of equations using the iterative solver on VE\n\n");

    /*..initialize result, search and residual vectors.................................	*/
    qam=0.;iam=0;
    for (i=0; i<neq; i++) {
        x[i] =  0.0;					/*..start vector x0 = 0....................	*/
        r[i] =  b[i];					/*..residual vector r0 = Ax+b = b..........	*/
        p[i] = -r[i];					/*..initial search vector..................	*/
        err=fabs(r[i]);
        if(err>1.e-20) {
            qam+=err;
            iam++;
        }
    }

    if(iam>0) {
        qam=qam/neq;
    } else {
        *niter=0;
        free(r);free(p);free(z);
        return;
    }

#ifdef SBLAS
    /* Set the number of OpenMP threads */
    //omp_set_num_threads(8);

    /* Creation of a handle from CSR format */
    int ierr;
    sblas_handle_t hdl; 
    ierr = sblas_create_matrix_handle_from_csc_rd(neq, neq, iz, ia, A, SBLAS_INDEXING_0, SBLAS_SYMMETRIC, &hdl);

    if (ierr != SBLAS_OK) {
        if (ierr == SBLAS_ERROR_INPUT) {
            printf("ERROR: sblas_create_matrix_handle with SBLAS_ERROR_INPUT\n");
        } else if (ierr == SBLAS_ERROR_MEMORY_VE) {
            printf("ERROR: sblas_create_matrix_handle with SBLAS_ERROR_MEMORY_VE\n");
        } else {
            printf("ERROR: sblas_create_matrix_handle with unknown error(%d)\n", ierr);
        } 
        exit(1);
    }

    /* Analysis of the sparse matrix A */
    ierr = sblas_analyze_mv_rd(SBLAS_NON_TRANSPOSE, hdl);
    if (ierr != SBLAS_OK) {
        if (ierr == SBLAS_ERROR_INPUT) {
            printf("ERROR: sblas_create_matrix_handle with SBLAS_ERROR_INPUT\n");
        } else if (ierr == SBLAS_ERROR_MEMORY_VE) {
            printf("ERROR: sblas_create_matrix_handle with SBLAS_ERROR_MEMORY_VE\n");
        } else {
            printf("ERROR: sblas_create_matrix_handle with unknown error(%d)\n", ierr);
        } 
        exit(1);
    }
#endif	

    /*else qam=0.01;*/
    /*..main iteration loop............................................................	*/
    for (k=1; k<=(*niter); k++, ncg++) {
        /*......inner product rT r......................................................... */
        InnerProduct(r,r,neq,&rr);
#if 0
        if ((ncg % 100) == 0) {
        //if (ncg < 20) {
            printf("iteration= %" ITGFORMAT ", error= %e, limit=%e\n",ncg,ram,c1*qam);
        }
#endif
        /*......If working with Penalty-terms for boundary conditions you can get nume-....	*/
        /*......rical troubles because RR becomes a very large value. With at least two....	*/
        /*......iterations the result may be better !!!....................................	*/
        /*......convergency check..........................................................	*/
        if (k!=1 && (ram<=c1*qam)) break;
        /*......new search vector..........................................................	*/
        if (k!=1) {
            ekm1 = rr/rro;
            for (i=0; i<neq; i++)	p[i] = ekm1*p[i]-r[i];
        }
        /*......matrix vector product A p = z..............................................	*/
#ifdef SBLAS
        ierr = sblas_execute_mv_rd(SBLAS_NON_TRANSPOSE, hdl, 1.0, p, 0.0, z);
        if (ierr != SBLAS_OK) {
            if (ierr == SBLAS_ERROR_INPUT) {
                printf("ERROR: sblas_create_matrix_handle with SBLAS_ERROR_INPUT\n");
            } else if (ierr == SBLAS_ERROR_MEMORY_VE) {
                printf("ERROR: sblas_create_matrix_handle with SBLAS_ERROR_MEMORY_VE\n");
            } else {
                printf("ERROR: sblas_create_matrix_handle with unknown error(%d)\n", ierr);
            } 
            exit(1);
        }
#else
        MatVecProduct(A,p,neq,ia,iz,z);
#endif
        /*......inner product pT z.........................................................	*/
        InnerProduct(p,z,neq,&pz);
        /*......step size..................................................................	*/
        qk = rr/pz;
        /*......approximated solution vector, residual vector..............................	*/
        ram=0.;
        for (i=0; i<neq; i++) {
            x[i] = x[i] + qk*p[i];
            r[i] = r[i] + qk*z[i];
            err=fabs(r[i]);
            if(err>ram) ram=err;
        }
        /*......store previous residual....................................................	*/
        rro = rr;
    }
    if(k==*niter) {
        printf("*ERROR in CG: no convergence;");
        exit(1);
    } 
    *eps = rr;						/*..return residual............................	*/
    *niter = ncg;					/*..return number of iterations................	*/
    /*..That's it......................................................................	*/

    free(r);free(p);free(z);
#ifdef SBLAS
    ierr = sblas_destroy_matrix_handle(hdl);
    if (ierr != SBLAS_OK) {
        if (ierr == SBLAS_ERROR_HND) {
            printf("ERROR: sblas_create_matrix_handle with SBLAS_ERROR_HND\n");
        } else {
            printf("ERROR: sblas_create_matrix_handle with unknown error(%d)\n", ierr);
        } 
        exit(1);
    }
#endif
    return;
}
/*---------------------------------------------------------------------------------	*/

/* **********************************************************************

   The (preconditioned) conjugate gradient solver

parameter:                                                                
A       compact row oriented storage of lower left of matrix A  
neq     order of A, number of equations                         
len     number of non zero entries in Matrix A                  
ia      column indices of corresponding elements in Matrix A    
iz      row indices (diagonal elements indices)               
x       solution vector                                       
b       right hand side                                       
eps     required accuracy -> residual                         
niter   maximum number of iterations -> number of iterations  
precFlg preconditioning flag                                  

The compact row oriented storage of sparse quadratic matrices is decsribed in
H.R. Schwarz: FORTRAN-Programme zur Methode der finiten Elemente, pp.66-67, 
Teubner, 1981

 ********************************************************************** 
 */
uint64_t cgve_solve(ITG neq, ITG len, uint64_t aptr, uint64_t aind, uint64_t val, uint64_t bptr, uint64_t xptr) {
    ITG* iz = (ITG*)aptr;
    ITG* ia = (ITG*)aind;
    double* A  = (double*)val;
    double* b  = (double*)bptr;
    double* x  = (double*)xptr;
    double *Factor=NULL,*x0; 

#if 0
    printf("INFO: cgsolver neq=%d, len=%d\n", neq, len);
    for (i=0; i<20; i++) {
        printf("DEBUG: (iz,ia,A)[%d]=(%d,%d,%lf)\n", i, iz[i], ia[i], A[i]);
    }
#endif

    /*  Scaling the equation system A x + b = 0  */
    Factor = (double*)calloc(sizeof(double), neq);
    Scaling(A,b,neq,ia,iz,Factor);

    /*  SOLVER/PRECONDITIONING TYPE  */
    /*  Conjugate gradient solver without preconditioning  */
    ITG niter=5000000;
    double eps=1.e-4;
    x0 = (double*)calloc(sizeof(double), neq);
    CG(A,x0,b,neq,len,ia,iz,&eps,&niter);

    /*  Backscaling of the solution vector  */
    for (int i=0; i<neq; i++) {
        x[i] = x0[i]*Factor[i];
    }
    free(x0);

    /*  That's it  */
    free(Factor);

    //printf("error condition (0=good, 1=bad) = %" ITGFORMAT "\n",ier);
    printf("# of iterations = %" ITGFORMAT "\n",niter);

    return 0;
}
