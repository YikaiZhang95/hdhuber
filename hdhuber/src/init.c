#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(huber_cd)(double *alpha, double *lam2, double *hval, 
  int *nobs,int *nvars, double *x, double *y, int *jd, int *pfncol, 
  double *pf, double *pf2, int *dfmax, int *pmax, int *nlam, double *flmin, 
  double *ulam, double *eps, int *isd, int *maxit, int *istrong, int *nalam, 
  double *b0, double *beta, int *ibeta, int *nbeta, double *alam, 
  int *npass, int *jerr);
// extern void F77_NAME(ssvm_pgd)(double *Xmat, int *nobs, int *np, 
//   double *y, int *nlam, double *ulam, double *eps, int *maxit, double *gam,
//   int *anlam, int *npass, int *jerr, double *btmat);

static const R_FortranMethodDef FortranEntries[] = {
    {"huber_cd",   (DL_FUNC) &F77_SUB(huber_cd),   28},
    // {"ssvm_pgd",   (DL_FUNC) &F77_SUB(ssvm_pgd),   13},
    {NULL, NULL, 0}
};

void R_init_hdhuber(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
