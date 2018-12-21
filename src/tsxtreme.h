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

#ifndef TSEXTREME_H
#define TSEXTREME_H 

void et_interface(double* data,
                   int* n, int* nlag, int* k, int* kred,
                   int* maxit, int* burn, int* thin,
                   int* adapt, int* batchsize,
                   double* sd_propa, double* sd_propb,
                   double* mu_prior, double* nu_prior, double* eta_prior,
                   int* mode, int* spec,
                   double* t_a, double* t_b, double* t_sig, double* t_mu,
                   double* t_w, double* t_gam,
                   int* t_ci, int* t_noo, int* t_noc, double* t_sd);

#endif
