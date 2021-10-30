// omp_critical.cpp
// compile with: /openmp
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define SIZE 10

int main() {
   int n = omp_get_thread_num();
   printf("%d", n);
   int N = omp_get_num_threads();
   printf("%d", N);
   double t1 = omp_get_wtime();
   printf("%f", t1);
}