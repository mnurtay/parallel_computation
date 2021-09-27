//
//  main.c
//  parallel computation
//
//  Created by Maksat Nurtay on 27.09.2021.
//

//#include "example/pragma_for.h"
//#include "example/available_threads.h"

#include "lab_2/matrix_multiplication.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, const char * argv[]) {
    matrix_multiplication();
    return 0;
}
