// Compile with: gcc -shared -fPIC timetravel.c -o timetravel.so
// Use with: export LD_PRELOAD="timetravel.so:${LD_PRELOAD}"

/* Adapted from datefudge
 * Copyright (C) 2001-2003, Matthias Urlichs <smurf@noris.de>
 * Copyright (C) 2008-2011, Robert Luberda <robert@debian.org>
 */

#define _GNU_SOURCE

#include <sys/file.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <assert.h>
#include <features.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>

static time_t travel = 1661140800; // date -d "2022-08-22" '+%s'

time_t real_time(time_t *x);
time_t real_gettimeofday(struct timeval *x, void *y);

static void set_travel(time_t *seconds)
{
    if (!seconds)
        return;
    *seconds = travel;
}

// Proxy to libc time function, setting/returning the "real" time.
time_t real_time(time_t *x)
{
    static time_t (*libc_time)(time_t *) = NULL;

    if (!libc_time)
        libc_time = (typeof(libc_time))dlsym(RTLD_NEXT, "time");
    return libc_time(x);
}

time_t time(time_t *x)
{
    time_t res;

    res = real_time(x);
    set_travel(x);
    set_travel(&res);
    return res;
}

// Proxy to libc __gettimeofday function, setting/returning the "real" time.
time_t real_gettimeofday(struct timeval *x, void *y)
{
    static int (*libc_gettimeofday)(struct timeval *, void *) = NULL;

    if (!libc_gettimeofday)
        libc_gettimeofday = (typeof(libc_gettimeofday))dlsym(RTLD_NEXT, "gettimeofday");
    return libc_gettimeofday(x, y);
}

int __gettimeofday(struct timeval *x, void *y)
{
    int res;

    res = real_gettimeofday(x, y);
    if (res)
        return res;
    set_travel(&x->tv_sec);
    return 0;
}

int gettimeofday(struct timeval *x, void *y)
{
    return __gettimeofday(x, y);
}

int clock_gettime(clockid_t x, struct timespec *y)
{
    static int (*libc_clock_gettime)(clockid_t, struct timespec *);
    int res;

    if (!libc_clock_gettime)
        libc_clock_gettime = (typeof(libc_clock_gettime))dlsym(RTLD_NEXT, __func__);
    res = (*libc_clock_gettime)(x, y);
    if (res || CLOCK_REALTIME != x)
        return res;
    set_travel(&y->tv_sec);
    return 0;
}
