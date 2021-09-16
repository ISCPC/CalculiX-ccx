/*
 * sxat.h
 *
 *  Created on: Jun 17, 2021
 *      Author: uno
 */

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
