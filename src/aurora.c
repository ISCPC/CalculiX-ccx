/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */
/*              Copyright (C) 2020 Shunji Uno                            */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#ifdef AURORA

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "aurora.h"

static uint64_t veo_handle = 0;
static struct veo_proc_handle *proc;
static struct veo_thr_ctxt *ctx;

static uint64_t aptr, aind, val;
static int mrow, ncol;

inline void set_matrixes(double *ad, double *au, double *adb, double *aub, double *sigma,
        ITG *icol, ITG *irow, ITG *neq, ITG *nzs, ITG *pointers, ITG *indice, double *value, long long ndim){
    ITG i,j,k,l;

    //printf("INFO: aptr: %d, aind: %lld sigma: %lf\n", *neq+1, ndim, *sigma);

    k=ndim;
    l=(*nzs);
    pointers[*neq]=ndim;

    if (*sigma==0.0) {
        for(i=*neq-1;i>=0;--i){
            for(j=0;j<icol[i];++j){
                indice[--k]=irow[--l]-1;
                value[k]=au[l];
            }
            pointers[i]=--k;
            indice[k]=i;
            value[k]=ad[i];
        }
    } else {
        for(i=*neq-1;i>=0;--i){
            for(j=0;j<icol[i];++j){
                indice[--k]=irow[--l]-1;
                value[k]=au[l]-(*sigma)*aub[l];
            }
            pointers[i]=--k;
            indice[k]=i;
            value[k]=ad[i]-(*sigma)*adb[i];
        }
    }
}

#define VEO_ALLOC_MEM(addr, size) veo_alloc_mem_wrapper((addr), (size), __LINE__)
#define VEO_WRITE_MEM(dst, src, size) veo_write_mem_wrapper((dst), (src), (size), __LINE__)
#define VEO_CALL_ASYNC(symname, argp) veo_call_async_wrapper((symname), (argp), __LINE__)
#define VEO_CALL_WAIT(reqid) veo_call_wait_wrapper((reqid), __LINE__)
#define VEO_CALL_SYNC(symname, argp) veo_call_sync((symname), (argp), __LINE__)

static inline void veo_alloc_mem_wrapper(uint64_t* addr, const size_t size, int line) {
    int rc;

    rc = veo_alloc_mem(proc, addr, size);
    if (rc != 0) {
        printf("ERROR: Cannot allocate memory on VE (line:%d)\n", line);
        exit(1);
    }
}

static inline void veo_write_mem_wrapper(uint64_t dst, const void* src, size_t size, int line) {
    int rc;

    rc = veo_write_mem(proc, dst, src, size);
    if (rc != 0) {
        printf("ERROR: Cannot write memory to VE (line:%d)\n", line);
        exit(1);
    }
}

static inline uint64_t veo_call_async_wrapper(const char* symname, struct veo_args* argp, int line) {
    uint64_t id;

    id = veo_call_async_by_name(ctx, veo_handle, symname, argp);
    if (id == VEO_REQUEST_ID_INVALID) {
        printf("ERROR: veo_call_async_by_name(\"%s\") failed (line:%d)\n", symname, line);
        exit(1);
    }
    return id;
}

static inline uint64_t veo_call_wait_wrapper(uint64_t reqid, int line) {
    uint64_t retval;
    
    int rc = veo_call_wait_result(ctx, reqid, &retval);
    if (rc != VEO_COMMAND_OK) {
        printf("ERROR: veo_call_wait_result failed with rc=%d (line:%d)\n", rc, line);
        exit(1);
    }
    return retval;
}

static inline uint64_t veo_call_sync(const char* symname, struct veo_args* argp, int line) {
    uint64_t id, retval;

    id = veo_call_async_by_name(ctx, veo_handle, symname, argp);
    if (id == VEO_REQUEST_ID_INVALID) {
        printf("ERROR: veo_call_async_by_name(\"%s\") failed (line:%d)\n", symname, line);
        exit(1);
    }

    int rc = veo_call_wait_result(ctx, id, &retval);
    if (rc != VEO_COMMAND_OK) {
        printf("ERROR: veo_call_wait_result(\"%s\") failed with rc=%d (line:%d)\n", symname, rc, line);
        exit(1);
    }
    return retval;
}

/*
 * HeteroSolver subroutines for SX-Aurora VE
 */
void aurora_hs_factor(double *ad, double *au, double *adb, double *aub, 
        double *sigma, ITG *icol, ITG *irow, ITG *neq, ITG *nzs) { 
    ITG *pointers=NULL;
    ITG *indice=NULL;
    double *value=NULL;
    long long ndim = (*neq)+(*nzs);

    if (veo_handle <= 0) {
        return;
    }

    printf(" Factoring the system of equations using the Aurora Heterosolver\n\n");
    NNEW(pointers,ITG,*neq+1);
    NNEW(indice,ITG,ndim);
    NNEW(value,double,ndim);
    set_matrixes(ad, au, adb, aub, sigma, icol, irow, neq, nzs, pointers, indice, value, ndim);

    mrow = (*neq);
    ncol = (*neq);

    // Call VE function named "factorize"
    struct veo_args *argp = veo_args_alloc();
    uint64_t retval;

    veo_args_set_i32(argp, 0, mrow);
    veo_args_set_i32(argp, 1, ncol);

    VEO_ALLOC_MEM(&aptr, (*neq+1)*sizeof(int));
    VEO_WRITE_MEM(aptr, pointers, (*neq+1)*sizeof(int));
    veo_args_set_i64(argp, 2, aptr);

    VEO_ALLOC_MEM(&aind, ndim*sizeof(int));
    VEO_WRITE_MEM(aind, indice, ndim*sizeof(int));
    veo_args_set_i64(argp, 3, aind);

    VEO_ALLOC_MEM(&val, ndim*sizeof(double));
    VEO_WRITE_MEM(val, value, ndim*sizeof(double));
    veo_args_set_i64(argp, 4, val);

    retval = VEO_CALL_SYNC("factorize", argp);
    if (retval !=0) {
        printf("ERROR: %s failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    } 
    veo_args_free(argp);

    SFREE(pointers);
    SFREE(indice);
    SFREE(value);

    return;
}

void aurora_hs_solve(double *b){
    /* Solution Phase */
    struct veo_args *argp = veo_args_alloc();
    uint64_t retval;
    uint64_t bptr, xptr;

    // Call VE function named "solve"
    VEO_ALLOC_MEM(&bptr, mrow*sizeof(double));
    VEO_WRITE_MEM(bptr, b, mrow*sizeof(double));
    veo_args_set_i64(argp, 0, bptr);

    VEO_ALLOC_MEM(&xptr, ncol*sizeof(double));
    veo_args_set_i64(argp, 1, xptr);

    retval = VEO_CALL_SYNC("solve", argp);
    if (retval !=0) {
        printf("ERROR: %s failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    } 
    veo_args_free(argp);

    veo_read_mem(proc, b, xptr, ncol*sizeof(double));
    veo_free_mem(proc, bptr);
    veo_free_mem(proc, xptr);

    return;
}

void aurora_hs_cleanup(){
    /* Handle Finalization */
    struct veo_args *argp = veo_args_alloc();
    uint64_t retval;

    // Call VE function named "finalize"
    retval = VEO_CALL_SYNC("finalize", argp);
    if (retval !=0) {
        printf("ERROR: %s failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    } 
    veo_args_free(argp);

    veo_free_mem(proc, aptr);
    veo_free_mem(proc, aind);
    veo_free_mem(proc, val);

    return;
}

void aurora_hs_main(double *ad, double *au, double *adb, double *aub, double *sigma,
        double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs){
    ITG *pointers=NULL;
    ITG *indice=NULL;
    double *value=NULL;
    long long ndim = (*neq)+(*nzs);

    if (veo_handle <= 0) {
        return;
    }

    printf(" Factoring the system of equations using the Aurora Heterosolver\n\n");
    NNEW(pointers,ITG,*neq+1);
    NNEW(indice,ITG,ndim);
    NNEW(value,double,ndim);
    set_matrixes(ad, au, adb, aub, sigma, icol, irow, neq, nzs, pointers, indice, value, ndim);

    mrow = (*neq);
    ncol = (*neq);

    // Call VE function named "factorize_solve"
    struct veo_args *argp1 = veo_args_alloc();
    uint64_t id, retval;
    uint64_t bptr, xptr;

    // Call VE function named "factorize"
    veo_args_set_i32(argp1, 0, mrow);
    veo_args_set_i32(argp1, 1, ncol);

    VEO_ALLOC_MEM(&aptr, (*neq+1)*sizeof(int));
    VEO_WRITE_MEM(aptr, pointers, (*neq+1)*sizeof(int));
    veo_args_set_i64(argp1, 2, aptr);

    VEO_ALLOC_MEM(&aind, ndim*sizeof(int));
    VEO_WRITE_MEM(aind, indice, ndim*sizeof(int));
    veo_args_set_i64(argp1, 3, aind);

    VEO_ALLOC_MEM(&val, ndim*sizeof(double));
    VEO_WRITE_MEM(val, value, ndim*sizeof(double));
    veo_args_set_i64(argp1, 4, val);

    id = VEO_CALL_ASYNC("factorize", argp1);

    // Call VE function named "solve"
    struct veo_args *argp2 = veo_args_alloc();

    VEO_ALLOC_MEM(&bptr, mrow*sizeof(double));
    VEO_WRITE_MEM(bptr, b, mrow*sizeof(double));
    veo_args_set_i64(argp2, 0, bptr);

    VEO_ALLOC_MEM(&xptr, ncol*sizeof(double));
    veo_args_set_i64(argp2, 1, xptr);

    retval = VEO_CALL_WAIT(id);
    if (retval !=0) {
        printf("ERROR: %s:factorize failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    } 
    id = VEO_CALL_ASYNC("solve", argp2);
    SFREE(pointers);
    SFREE(indice);
    SFREE(value);
    veo_args_clear(argp1);

    retval = VEO_CALL_WAIT(id);
    if (retval !=0) {
        printf("ERROR: %s:solve failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    } 
    veo_read_mem(proc, b, xptr, ncol*sizeof(double));

    // Call VE function named "finalize"
    id = VEO_CALL_ASYNC("finalize", argp1);
    veo_args_free(argp2);
    veo_free_mem(proc, bptr);
    veo_free_mem(proc, xptr);
    retval = VEO_CALL_WAIT(id);
    if (retval !=0) {
        printf("ERROR: %s:finalize failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    } 
    veo_args_free(argp1);

    veo_free_mem(proc, aptr);
    veo_free_mem(proc, aind);
    veo_free_mem(proc, val);

    return;
}

/*
 * CG solver subroutines for SX-Aurora VE
 */
void aurora_cg_factor(double *ad, double *au, double *adb, double *aub, 
        double *sigma, ITG *icol, ITG *irow, ITG *neq, ITG *nzs) { 
    printf("ERROR: %s has not been supported yet.\n");
    exit(1);
}

void aurora_cg_solve(double *b){
    printf("ERROR: %s has not been supported yet.\n");
    exit(1);
}

void aurora_cg_cleanup(){
    printf("ERROR: %s has not been supported yet.\n");
    exit(1);
}

void aurora_cg_main(double *ad, double *au, double *adb, double *aub, double *sigma,
        double *b, ITG *icol, ITG *irow, ITG *neqp, ITG *nzsp){
    ITG *pointers=NULL;
    ITG *indice=NULL;
    double *value=NULL;
    ITG ndim = (*neqp)+(*nzsp);

    if (veo_handle <= 0) {
        return;
    }

    printf(" Factoring the system of equations using the CG/VE solver\n\n");
    NNEW(pointers,ITG,*neqp+1);
    NNEW(indice,ITG,ndim);
    NNEW(value,double,ndim);
    set_matrixes(ad, au, adb, aub, sigma, icol, irow, neqp, nzsp, pointers, indice, value, ndim);
    ITG neq = *neqp;

    // Call VE function named "cgve_solve"
    // uint64_t cgve_solve(ITG neq, ITG len, uint64_t aptr, uint64_t aind, uint64_t val, uint64_t bptr, uint64_t xptr);
    struct veo_args *argp = veo_args_alloc();
    uint64_t retval;
    uint64_t bptr, xptr;

    veo_args_set_i32(argp, 0, neq);
    veo_args_set_i32(argp, 1, ndim);

    VEO_ALLOC_MEM(&aptr, (neq+1)*sizeof(int));
    VEO_WRITE_MEM(aptr, pointers, (neq+1)*sizeof(int));
    veo_args_set_i64(argp, 2, aptr);

    VEO_ALLOC_MEM(&aind, ndim*sizeof(int));
    VEO_WRITE_MEM(aind, indice, ndim*sizeof(int));
    veo_args_set_i64(argp, 3, aind);

    VEO_ALLOC_MEM(&val, ndim*sizeof(double));
    VEO_WRITE_MEM(val, value, ndim*sizeof(double));
    veo_args_set_i64(argp, 4, val);

    VEO_ALLOC_MEM(&bptr, neq*sizeof(double));
    VEO_WRITE_MEM(bptr, b, neq*sizeof(double));
    veo_args_set_i64(argp, 5, bptr);

    VEO_ALLOC_MEM(&xptr, neq*sizeof(double));
    veo_args_set_i64(argp, 6, xptr);
    retval = VEO_CALL_SYNC("cgve_solve", argp);
    if (retval !=0) {
        printf("ERROR: %s failed with %ld\n", __FUNCTION__, retval);
        exit(1);
    } 
    SFREE(pointers);
    SFREE(indice);
    SFREE(value);
    veo_args_free(argp);
    veo_read_mem(proc, b, xptr, neq*sizeof(double));
    veo_free_mem(proc, bptr);
    veo_free_mem(proc, xptr);
    veo_free_mem(proc, aptr);
    veo_free_mem(proc, aind);
    veo_free_mem(proc, val);

    return;
}
#define CCX_VEO_DEFAULT_LIBRARY_PATH "/usr/local/lib/CalculiX/libccx.so"

int aurora_init() {
    proc = veo_proc_create(0);
    if (proc == NULL) {
        printf("ERROR: Creating VE offload process failed.\n");
        fflush(stdout);
        return -1;
    }

    char* libpath = getenv("CCX_VEO_LIBRARY_PATH");
    if (libpath == NULL) {
        libpath = CCX_VEO_DEFAULT_LIBRARY_PATH;
    }

    veo_handle = veo_load_library(proc, libpath);
    if (veo_handle == 0) {
        printf("ERROR: Loading %s failed.\n", libpath);
        fflush(stdout);
        return -1;
    }

    ctx = veo_context_open(proc);
    if (ctx == NULL) {
        printf("ERROR: opening VE offload context failed.\n");
        fflush(stdout);
        return -1;
    }

    return 0;
}

int aurora_fini() {
    if (veo_handle != 0) {
        veo_handle = 0;
        return veo_context_close(ctx);
    }
    return -1;
}
#endif /* AURORA */
