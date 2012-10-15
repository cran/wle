#include <R.h>
#include <Rmath.h>

double F77_SUB(dgammac)(double *x, double *shape, double *scale, int *give_log) { 
  return(dgamma(*x, *shape, *scale, *give_log)); }
