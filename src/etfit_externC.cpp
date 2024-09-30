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


#include "etfit.h"
#include <R.h>
#include<ctime>

extern "C" {
    void et_interface(double* data,
                   int* n, int* nlag, int* k, int* kred,
                   int* maxit, int* burn, int* thin,
                   int* adapt, int* batchsize,
                   double* sd_propa, double* sd_propb,
                   double* mu_prior, double* nu_prior, double* eta_prior,
                   int* mode, int* spec,
                   double* t_a, double* t_b, double* t_sig, double* t_mu,
                   double* t_w, double* t_gam, int* t_ci, int* t_noo, int* t_noc, double* t_sd){

        //////////////////////////////////////////////////
        // FEED C++ CLASS & LAUNCH COMPUTATIONS

        std::clock_t start=std::clock();//start chrono

        tsxtreme::debmode modeCpp;
        switch(*mode){
        case 0: modeCpp = tsxtreme::debug;
            break;
        case 1: modeCpp = tsxtreme::normal;
            break;
        case 2: modeCpp = tsxtreme::silent;
            break;
        default: Rf_error("bad integer initialisation value for _mode_ in [et_interface()]");
        }
        tsxtreme::submodel specCpp;
        switch(*spec){
        case 0: specCpp = tsxtreme::univariateMixture;
            break;
        case 1: specCpp = tsxtreme::firstOrderMarkov;
            break;
        case 2: specCpp = tsxtreme::none;
            break;
        default: Rf_error("bad integer initialisation value for _spec_ in [et_interface()]");
        }

        ETfit fit(data, n, nlag, k, kred, maxit, burn, thin, adapt, batchsize,
                  sd_propa, sd_propb, mu_prior, nu_prior, eta_prior, modeCpp, specCpp);
        fit.run();
        const std::vector<ETpar> tr = fit.getTraces();



        //////////////////////////////////////////////////
        // PRINT TIME INFO
        double duration = (std::clock()-start) / (double) CLOCKS_PER_SEC;
        if(*mode!=2) Rprintf("Running time: %u min %u sec. Per iteration: %.2f msec.",
                        ((int)duration)/60,((int)duration)%60,1000*duration/(*maxit));



        //////////////////////////////////////////////////
        // FEED ARGUMENTS FOR R

        unsigned int const len = ((*maxit) - (*burn))/(*thin);//actual length after burn-in and thinning
        unsigned int const d_max = (specCpp == tsxtreme::none)?(*nlag):1;

        for(unsigned int it = 0; it < len; it++){
            for(int d = 0; d < (*nlag); d++){
                t_a[d*len + it] = tr[it].a[d];
                t_b[d*len + it] = tr[it].b[d];
            }
            for(int c = 0; c < (*kred); c++){
                t_w[c*len + it]   = tr[it].w[c];
                t_noo[c*len + it] = tr[it].noo[c];
                for(int d = 0; d < (*nlag); d++){
                    t_mu[d*(*kred)*len + c*len + it]  = tr[it].mu[c][d];
                    t_sig[d*(*kred)*len + c*len + it] = tr[it].sig[c][d];
                }
            }
            for(int i = 0; i < (*n); i++){
                t_ci[i*len + it] = tr[it].ci[i];
            }
            t_gam[it] = tr[it].gam;
            t_noc[it] = tr[it].noc;
            for(unsigned int reg = 0; reg < 8; reg++){
                for(unsigned int d = 0; d < d_max; d++){
                    t_sd[d*8*len+reg*len+it] = tr[it].sd[d][reg];//flip dimensions here
                }
            }
        }
    }
}

