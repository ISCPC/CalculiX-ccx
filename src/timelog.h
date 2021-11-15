/*     CALCULIX - A 3-dimensional finite element program                 */
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

#ifndef _TIMELOG_H_
#define _TIMELOG_H_

#ifdef _TIMELOG
#ifndef __USE_POSIX199309
#define __USE_POSIX199309
#endif
#include <time.h>

#define TIMELOG(tl) timelog_t tl;

#define TIMELOG_START(tl)   _timelog_start(&tl)
#define TIMELOG_END(tl, str)  _timelog_end(&tl, str)
#define TIMELOG_GETTIME(t,tl)  (t = _timelog_gettime(&tl))

typedef struct timelog {
    struct timespec stv;
    struct timespec etv;
} timelog_t;

static inline void _timelog_start(timelog_t* tl) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &(tl->stv));
}

static inline void _timelog_end(timelog_t* tl, const char* str) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &(tl->etv));
    double diff = (double)(tl->etv.tv_sec - tl->stv.tv_sec)
        + ((double)(tl->etv.tv_nsec - tl->stv.tv_nsec))/1000000000.0;
    printf("TIME: %s : %12.6lf [sec]\n", str, diff);
}

static inline double _timelog_gettime(timelog_t* tl) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &(tl->etv));
    double diff = (double)(tl->etv.tv_sec - tl->stv.tv_sec)
        + ((double)(tl->etv.tv_nsec - tl->stv.tv_nsec))/1000000000.0;
    return diff;
}
#else
#define TIMELOG(tl)

#define TIMELOG_START(tl)
#define TIMELOG_END(tl, str)
#define TIMELOG_GETTIME(t, tl)
#endif

#endif /* _TIMELOG_H_ */
