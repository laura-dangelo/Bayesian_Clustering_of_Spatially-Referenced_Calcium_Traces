#ifndef MVTNORM_H
#define MVTNORM_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

void F77_NAME(mvtdst)(int *n, int *nu, double *lower, double *upper,
                      int *infin, double *corr, double *delta,
                      int *maxpts, double *abseps, double *releps,
                      double *error, double *value, int *inform);

#endif
