/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */
/*              Copyright (C) 2021 Shunji Uno                            */
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

#ifndef TESTS_AURORA_H_
#define TESTS_AURORA_H_
#include <stdint.h>
#include "vesolver.h"
#include "CalculiX.h"

//typedef int32_t ITG;

#define SOLVER_TYPE_HS        0
#define SOLVER_TYPE_CG        1
#define SOLVER_TYPE_BICGSTAB2 2

/*
 * Main function
 */
#ifdef __cplusplus
extern "C" {
#endif

int sxat_ve_init();
void sxat_ve_fini();
void sxat_ve_factor(double *ad, double *au, double *adb, double *aub,
        const double sigma, ITG *icol, ITG *irow, const ITG neq, const ITG nzs,
	    const ITG symmetryflag, const ITG inputformat, ITG *jq, const ITG nzs3,
        const int solvertype);
void sxat_ve_solve(double *b);
void sxat_ve_cleanup();

void sxat_ve_main(double *ad, double *au, double *adb, double *aub, const double sigma,
        double *b, ITG *icol, ITG *irow, const ITG neq, const ITG nzs,
	    const ITG symmetryflag, const ITG inputformat, ITG *jq, const ITG nzs3,
        const int solvertype);

#ifdef __cplusplus
}
#endif

#endif /* TESTS_AURORA_H_ */
