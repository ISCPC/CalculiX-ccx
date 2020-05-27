/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */
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

#include <ve_offload.h>
#include "CalculiX.h"

int aurora_init();
int aurora_fini();

void aurora_hs_main(double *ad, double *au, double *adb, double *aub, double *sigma,
        double *b, ITG *icol, ITG *irow, ITG *neq, ITG *nzs);

void aurora_hs_factor(double *ad, double *au, double *adb, double *aub, 
        double *sigma,ITG *icol, ITG *irow, ITG *neq, ITG *nzs);

void aurora_hs_solve(double *b);

void aurora_hs_cleanup();

void aurora_cg_factor(double *ad, double *au, double *adb, double *aub,
        double *sigma, ITG *icol, ITG *irow, ITG *neq, ITG *nzs);
void aurora_cg_solve(double *b);
void aurora_cg_cleanup();
void aurora_cg_main(double *ad, double *au, double *adb, double *aub, double *sigma,
        double *b, ITG *icol, ITG *irow, ITG *neqp, ITG *nzsp);
