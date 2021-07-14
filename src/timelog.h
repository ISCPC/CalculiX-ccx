/*
MIT License

Copyright (c) 2021 Shunji Uno <shunji_uno@iscpc.jp>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef _TIMELOG_H_
#define _TIMELOG_H_

#ifdef _TIMELOG
#define __USE_POSIX199309
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
