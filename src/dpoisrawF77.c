#include <R.h>
#include <Rmath.h>
#include "nmath.h"
#include "dpq.h"

double F77_SUB(dpoisraw)(double *x, double *lambda, int *give_log)  { 
  return(dpois_raw(*x, *lambda, *give_log)); }
