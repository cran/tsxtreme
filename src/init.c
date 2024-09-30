/*
 *  tsxtreme : Bayesian Modelling of Extremal Dependence in Time Series
 *  Copyright (C) 2017-2018   Thomas Lugrin
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available at
 *  https://www.R-project.org/Licenses/
 */

#include "tsxtreme.h"
#include <stdlib.h> // for NULL
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n, argTypes) {#name, (DL_FUNC) &name, n, argTypes}

static R_NativePrimitiveArgType et_interface_t[] = {
   REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,
   INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
   REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
   INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
   INTSXP, REALSXP
};

const static R_CMethodDef R_CDef[] = {
   CALLDEF(et_interface, 27, et_interface_t),
   {NULL, NULL, 0}
};

void R_init_tsxtreme(DllInfo *dll)
{
   R_registerRoutines(dll, R_CDef, NULL, NULL, NULL);
   R_useDynamicSymbols(dll, FALSE);
   R_forceSymbols(dll, TRUE);
}
 
