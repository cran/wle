#include <R.h>
#include <Rmath.h>
     
void F77_SUB(ddbeta)(double *x, double *a, double *b, int *give_log, double *d) {
  *d = dbeta(*x, *a, *b, *give_log);
}
