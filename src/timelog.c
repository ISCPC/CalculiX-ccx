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
#include "timelog.h"

#ifdef _TIMELOG
void _timelog_start(timelog_t* tl) {
    gettimeofday(&(tl->stv), NULL);
}

void _timelog_end(timelog_t* tl, char* str) {
    gettimeofday(&(tl->etv), NULL);
    float diff = (float)(tl->etv.tv_sec - tl->stv.tv_sec)
        + ((float)(tl->etv.tv_usec - tl->stv.tv_usec))/1000000.0;
    printf("TIME: %s : %f [sec]\n", str, diff);
}

float _timelog_gettime(timelog_t* tl) {
    gettimeofday(&(tl->etv), NULL);
    float diff = (float)(tl->etv.tv_sec - tl->stv.tv_sec)
        + ((float)(tl->etv.tv_usec - tl->stv.tv_usec))/1000000.0;
    return diff;
}
#endif
