//
//  available_threads.c
//  parallel computation
//
//  Created by Maksat Nurtay on 27.09.2021.
//

#include "available_threads.h"
#include <omp.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#define THREAD_NUM 4

void runAvailableThreadsExample(void) {
    omp_set_num_threads(THREAD_NUM); // set number of threads in "parallel" blocks
#pragma omp parallel
        sleep(omp_get_thread_num()); // do this to avoid race condition while printing
        printf("Number of available threads: %d\n", omp_get_num_threads());
        // each thread can also get its own number
        printf("Current thread number: %d\n", omp_get_thread_num());
}
