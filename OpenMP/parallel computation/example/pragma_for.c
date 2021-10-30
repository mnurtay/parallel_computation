//
//  pragma_for.c
//  parallel computation
//
//  Created by Maksat Nurtay on 27.09.2021.
//

#include "pragma_for.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define ARRAY_SIZE 100000000
#define ARRAY_VALUE 1231

void runPragmaForExample(void) {
    int arr[ARRAY_SIZE];
    
    #pragma omp parallel
    {
        #pragma omp for
        {
            for(int i = 0; i < ARRAY_SIZE; i++) {
                arr[i] = arr[i] / arr[i] + arr[i] / 5 - 14;
            }
        }
    }
    
}
