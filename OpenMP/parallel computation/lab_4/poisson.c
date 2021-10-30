////
////  poisson.c
////  parallel computation
////
////  Created by Maksat Nurtay on 15.10.2021.
////
//
//#include "poisson.h"
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <time.h>
//#include <omp.h>
//
//#define NX 161
//#define NY 161
//
//int main ( int argc, char *argv[] );
//double r8mat_rms ( int nx, int ny, double a[NX][NY] );
//void rhs ( int nx, int ny, double f[NX][NY] );
//void sweep ( int nx, int ny, double dx, double dy, double f[NX][NY],
//int itold, int itnew, double u[NX][NY], double unew[NX][NY] );
//double u_exact ( double x, double y );
//double uxxyy_exact ( double x, double y );
//
//
//int main (int argc, char *argv[]) {
//    double diff, dx, dy;
//    double f[NX][NY], error;
//    int i, j;
//    int itnew, itold;
//    int nx = NX, ny = NY;
//    double tolerance = 0.000001;
//    double u[NX][NY];
//    double u_norm;
//    double udiff[NX][NY];
//    double uexact[NX][NY];
//    double unew[NX][NY];
//    double unew_norm;
//    double wtime;
//    double x, y;
//
//    dx = 1.0 / (double) (nx - 1);
//    dy = 1.0 / (double) (ny - 1);
//
//    printf("The number of processors is %d\n", omp_get_num_procs());
//
//    rhs(nx, ny, f);
//    for (j = 0; j < ny; j++) {
//        for (i = 0; i < nx; i++) {
//            if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) {
//                unew[i][j] = f[i][j];
//            } else {
//                unew[i][j] = 0.0;
//            }
//        }
//    }
//    unew_norm = r8mat_rms ( nx, ny, unew );
//
//    for (j = 0; j < ny; j++) {
//        y = (double) (j) / (double) (ny - 1);
//        for (i = 0; i < nx; i++) {
//            x = (double) (i) / (double) (nx - 1);
//            uexact[i][j] = u_exact(x, y);
//        }
//    }
//    u_norm = r8mat_rms ( nx, ny, uexact );
//
//    for (j = 0; j < ny; j++) {
//        for (i = 0; i < nx; i++) {
//            udiff[i][j] = unew[i][j] - uexact[i][j];
//        }
//    }
//    error = r8mat_rms(nx, ny, udiff);
//    wtime = omp_get_wtime();
//
//    itnew = 0;
//
//    for ( ; ; ) {
//        itold = itnew;
//        itnew = itold + 500;
//
//        sweep (nx, ny, dx, dy, f, itold, itnew, u, unew);
//
//        u_norm = unew_norm;
//        unew_norm = r8mat_rms (nx, ny, unew);
//
//        for (j = 0; j < ny; j++) {
//            for (i = 0; i < nx; i++) {
//                udiff[i][j] = unew[i][j] - u[i][j];
//            }
//        }
//        diff = r8mat_rms (nx, ny, udiff);
//
//        for (j = 0; j < ny; j++) {
//            for (i = 0; i < nx; i++) {
//                udiff[i][j] = unew[i][j] - uexact[i][j];
//            }
//        }
//        error = r8mat_rms (nx, ny, udiff);
//        if (diff <= tolerance) break;
//    }
//
//    wtime = omp_get_wtime ( ) - wtime;
//    printf("Elapsed seconds = %g\n", wtime);
//    return 0;
//}
//
//double r8mat_rms (int nx, int ny, double a[NX][NY]) {
//    int i, j;
//    double v = 0.0;
//    for (j = 0; j < ny; j++) {
//        for (i = 0; i < nx; i++) {
//            v = v + a[i][j] * a[i][j];
//        }
//    }
//    return sqrt(v / (double) (nx * ny));
//}
//
//void rhs(int nx, int ny, double f[NX][NY]) {
//    double x, y, fnorm;
//    int i, j;
//    for (j = 0; j < ny; j++) {
//        y = (double) (j) / (double) (ny - 1);
//        for (i = 0; i < nx; i++) {
//            x = (double) (i) / (double) (nx - 1);
//            if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) {
//                f[i][j] = u_exact(x, y);
//            } else {
//                f[i][j] = -uxxyy_exact(x, y);
//            }
//        }
//    }
//    fnorm = r8mat_rms ( nx, ny, f );
//    return;
//}
//
//
//void sweep(int nx, int ny, double dx, double dy, double f[NX][NY], int itold, int itnew, double u[NX][NY], double unew[NX][NY]) {
//    int i;
//    int it;
//    int j;
//
//# pragma omp parallel shared (dx, dy, f, itnew, itold, nx, ny, u, unew) private (i, it, j)
//    for (it = itold + 1; it <= itnew; it++) {
//# pragma omp for
//        for (j = 0; j < ny; j++) {
//            for (i = 0; i < nx; i++) {
//                u[i][j] = unew[i][j];
//            }
//        }
//# pragma omp for
//        for (j = 0; j < ny; j++) {
//            for (i = 0; i < nx; i++) {
//                if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1) {
//                    unew[i][j] = f[i][j];
//                } else {
//                    unew[i][j] = 0.25 * (u[i-1][j] + u[i][j+1] + u[i][j-1] + u[i+1][j] + f[i][j] * dx * dy);
//                }
//            }
//        }
//  }
//  return;
//}
//
//
//double u_exact(double x, double y) {
//    double r8_pi = 3.141592653589793;
//    return sin (r8_pi * x * y);
//}
//
//
//// (d/dx * d/dx + d/dy * d/dy)
//double uxxyy_exact(double x, double y) {
//    double r8_pi = 3.141592653589793;
//    return -r8_pi * r8_pi * (x * x + y * y) * sin (r8_pi * x * y);
//}
