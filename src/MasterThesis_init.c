#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _MasterThesis_dexp_rcpp(SEXP, SEXP, SEXP);
extern SEXP _MasterThesis_dgamma_rcpp(SEXP, SEXP, SEXP);
extern SEXP _MasterThesis_dlnorm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _MasterThesis_dvm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _MasterThesis_dweibull_rcpp(SEXP, SEXP, SEXP);
extern SEXP _MasterThesis_dwrpcauchy_rcpp(SEXP, SEXP, SEXP);
extern SEXP _MasterThesis_mllk_cpp();

static const R_CallMethodDef CallEntries[] = {
  {"_MasterThesis_dexp_rcpp",       (DL_FUNC) &_MasterThesis_dexp_rcpp,       3},
  {"_MasterThesis_dgamma_rcpp",     (DL_FUNC) &_MasterThesis_dgamma_rcpp,     3},
  {"_MasterThesis_dlnorm_rcpp",     (DL_FUNC) &_MasterThesis_dlnorm_rcpp,     3},
  {"_MasterThesis_dvm_rcpp",        (DL_FUNC) &_MasterThesis_dvm_rcpp,        3},
  {"_MasterThesis_dweibull_rcpp",   (DL_FUNC) &_MasterThesis_dweibull_rcpp,   3},
  {"_MasterThesis_dwrpcauchy_rcpp", (DL_FUNC) &_MasterThesis_dwrpcauchy_rcpp, 3},
  {"_MasterThesis_mllk_cpp",        (DL_FUNC) &_MasterThesis_mllk_cpp,        0},
  {NULL, NULL, 0}
};

void R_init_MasterThesis(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}