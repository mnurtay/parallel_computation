//
//  sorting_ algorithm.c
//  parallel computation
//
//  Created by Maksat Nurtay on 04.10.2021.
//

#include "sorting_algorithm.h"

#define N 5000

void sorting_algorithm(void) {
    int a[N], sorted[N];
    
//    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        a[i] = rand() % 1000 + 1;
        sorted[i] = 0;
    }
    
    int i, j, count;
    omp_set_num_threads(50);
    double start_time = omp_get_wtime();
    #pragma omp parallel private(i, j, count)
    {
        #pragma omp for
        for (i = 0; i < N; i++) {
            count = 0;
            for (j = 0; j < N; j++) {
                if (a[i] > a[j]) count++;
            }
            while (sorted[count] != 0) count++;
            sorted[count] = a[i];
        }
    }
    
    double end_time = omp_get_wtime();
    double time_used = end_time - start_time;
    printf("Parallel time: %f s\n", time_used);
    
    for (i = 0; i < N; i++) printf("%d ", sorted[i]);
}
