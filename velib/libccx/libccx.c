/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 2020 Shunji Uno                            */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; either version 2 of    */
/*     the License, or (at your option) any later version.               */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <omp.h>
#include <sys/time.h>
#include <heterosolver.h>

#define   NRHS   1   /* The number of right-hand side vectors */

HS_handle_t hnd;
HS_int_t *iaptr;
HS_int_t *iaind;
double   *aval;

uint64_t factorize(int row, int col, uint64_t aptr, uint64_t aind, uint64_t val) {
    const HS_int_t mrow = row;
    const HS_int_t ncol = col;
    iaptr = (HS_int_t*)aptr;
    iaind = (HS_int_t*)aind;
    aval  = (double*)val;

    //struct timeval tv1, tv2;

    HS_int_t ierr = HS_init_handle(&hnd, mrow, ncol, HS_SYMMETRIC, HS_CSC);

    omp_set_num_threads(8);

    ierr = HS_set_option(hnd, HS_ORDP, HS_ORDP_METIS);
    if (ierr != HS_RESULT_OK) {
        fprintf(stderr, "ERROR: HS_set_option failed with %d.\n", ierr);
        return ierr;
    }

    //gettimeofday(&tv1, NULL);
    //printf("Call HS_preprocess(): %d, %d, ...\n", iaptr[1], iaind[1]);
    ierr = HS_preprocess_rd(hnd, iaptr, iaind, aval);
    if (ierr) {
        return ierr;
    }

    //printf("Call HS_factorize(): %lf\n", aval[0]);
    ierr = HS_factorize_rd(hnd, iaptr, iaind, aval);

    //gettimeofday(&tv2, NULL);
    //printf("TIME: HS_factorize : %f [sec]\n",
    //    (float)(tv2.tv_sec - tv1.tv_sec)
    //    + ((float)(tv2.tv_usec - tv1.tv_usec))/1000000.0);

    return ierr;
}

uint64_t solve(uint64_t bptr, uint64_t xptr) {
    double* b = (double*)bptr;
    double* x = (double*)xptr;
    //double res = 1.0e-13;  
    double res = 1.0e-10;  

    //struct timeval tv1, tv2;

    omp_set_num_threads(8);

    //gettimeofday(&tv1, NULL);
    //printf("Call HS_solve()\n");
    //for (int i=0; i<7; i++) { x[i] = b[i]; }
    HS_int_t ierr = HS_solve_rd(hnd, iaptr, iaind, aval, NRHS, b, x, &res);
    //gettimeofday(&tv2, NULL);
    //printf("TIME: HS_solve_rd : %f [sec]\n",
    //    (float)(tv2.tv_sec - tv1.tv_sec)
    //    + ((float)(tv2.tv_usec - tv1.tv_usec))/1000000.0);
    if (ierr == HS_ERROR_ACCURACY) {
        if (res < 1.0e-5) {
            printf("WARNING:HS_solve_rd: HS_ERROR_ACCURACY(res=%g)\n", res);
            fflush(stdout);
            return 0;
        }
        printf("WARNING:HS_solve_rd: HS_ERROR_ACCURACY(res=%g)\n", res);
        fflush(stdout);
    }
    return ierr;
}

uint64_t finalize() {
    //printf("Call HS_finalize()\n");
    HS_int_t ierr = HS_finalize_handle(hnd);
    return ierr;
}

uint64_t factorize_solve(int row, int col, uint64_t aptr, uint64_t aind, uint64_t val, uint64_t bptr, uint64_t xptr) {
    HS_int_t ierr = factorize(row, col, aptr, aind, val);

    if (ierr == 0) {
        ierr = solve(bptr, xptr);
    }

    if (ierr == 0) {
        ierr = finalize();
    }

    return ierr;
}
