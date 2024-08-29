#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "clogic.h"

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
/*
extern void clogreg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
*/

static const R_CMethodDef CEntries[] = {
    {"clogreg", (DL_FUNC) &clogreg, 15},
    {NULL, NULL, 0}
};

void R_init_LogicReg(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
