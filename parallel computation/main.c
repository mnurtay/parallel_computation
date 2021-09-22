//
//  main.c
//  parallel computation
//
//  Created by Maksat Nurtay on 21.09.2021.
//
#include <omp.h>

#include <stdlib.h>
#include <stdio.h>

int main(int argc, const char * argv[]) {
    #pragma omp parallel
    {
        printf("Hello, World!\n");
    }
    return 0;
}
