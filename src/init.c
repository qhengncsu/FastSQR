#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP PlinkMultvC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP solve(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP testwz(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP testquantile(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP max_lambda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"PlinkMultvC", (DL_FUNC) &PlinkMultvC,  9},
    {"solve",       (DL_FUNC) &solve,       37},
    {"testwz",      (DL_FUNC) &testwz,       5},
    {"testquantile",(DL_FUNC) &testquantile, 5},
    {"max_lambda",(DL_FUNC) &max_lambda, 14},
    {NULL, NULL, 0}
};

void R_init_myglmnet(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
