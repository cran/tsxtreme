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

#ifndef ETFIT_H
#define ETFIT_H

#define R_USE_C99_IN_CXX

#include<vector>
#include<R.h>
#include<Rmath.h>
#include<cstdarg>

//////////////////////////////////////////////////
// STRUCTURE TO BE USED WITHIN THE ETFIT CLASS
// ...either for current parameters or to store traces

struct ETpar{
    std::vector<double> a;                    // alpha
    std::vector<double> b;                    // beta
    std::vector<std::vector<double> > sd;    // standard deviation for alpha and beta proposals (RAMA)
    std::vector<std::vector<double> > mu;     // components' means
    std::vector<std::vector<double> > sig;    // components' standard deviation
    std::vector<double> w;                    // components' weights
    double gam;                               // precision parameter for DP
    std::vector<unsigned int> ci;     // component index attached to each observation
    std::vector<unsigned int> noo;    // number of observations per component
    unsigned int noc;                 // number of (non-empty) components
};


//////////////////////////////////////////////////
// TYPEDEFS
namespace tsxtreme {
enum debmode{
    debug,       // print values of parameters at each step
    normal,      // print status regularly & end-of-process summary
    silent       // print nothing
};

enum submodel{// structure of the parameters of H+T model
    univariateMixture,   // alpha=beta=0, fitting a kind of kernel smoothing (e.g. galaxies)
    firstOrderMarkov,    // alpha geometrically decreasing, beta constant across lags
    none                 // no constraints
};

enum algotype{
    conditional, // blocked Gibbs type
    marginal     // R. M. Neal type (NOT IMPLEMENTED YET!)
};
}


//////////////////////////////////////////////////
// CLASS ETFIT
// no hyper-prior on the components' means (for now)

class ETfit
{
public:
    ETfit(double *data, int *n, int *nlag,
          int *k, int *kred, int *maxit,
          int *burn, int *thin, int *adapt, int* batchsize,
          double *sd_propa, double *sd_propb,
          double *mu_prior, double *nu_prior, double *eta_prior,
          tsxtreme::debmode const& mode, tsxtreme::submodel const& spec);
    ~ETfit();


    //--------------------//
    // CLASS INTERFACE

    void run(tsxtreme::algotype const& type=tsxtreme::conditional);
    std::vector<ETpar> const& getTraces() const;// return by reference for memory efficiency



    //--------------------//
    // STATS METHODS
    static int rmult(std::vector<double> const& probs, const double &sum);
    double mean(std::vector<double> const& x) const;
    double var(std::vector<double> const& x) const;
    double cov(std::vector<double> const& x, std::vector<double> const& y) const;



private:
    //--------------------//
    // ATTRIBUTES

    // posterior samples
    ETpar curr;            //accepted params at latest state (t or t-1)
    std::vector<ETpar> traces;
    std::vector<std::vector<std::vector<double> > > acc_a; // acceptance rates for alpha, beta
    std::vector<std::vector<std::vector<double> > > acc_b;
    // input
    std::vector<std::vector<double> > data;
    const unsigned int n;           // number of rows in _data_
    const unsigned int nlag;        // number of columns -1 in _data_ (number of time lags)
    const unsigned int k;           // number of components
    const unsigned int kred;        // only first _kred_ components' parameters are stored in _traces_
    const unsigned int maxit;
    const unsigned int burn;        // nbr of sweeps burnt
    const unsigned int thin;        // 1 in every _thin_ sweeps saved
    const unsigned int maxadapt;    // beyond this sweep, keep sd fixed
    const unsigned int batchsize;  // length of batches on  which mean of acceptance ratios is computed (RAMA)
    // prior (hyper-)parameters
    double mu[2];    // for the components' means
    double nu[2];    // for the components' variances
    double eta[2];   // for the precision parameter
    // run mode, constants, & others
    std::vector<double> V;         // for stick-breaking procedure
    double sumV;                   // sum_{c=1}^{k-1} log(1-V_c)
    unsigned int nbswaps1;         // counts number of swaps of type 1 (output just through rout)
    unsigned int nbswaps2;         // idem for swaps of type 2
    const tsxtreme::debmode mode;
    const tsxtreme::submodel spec;
    const double tol;              // when computing bounds on (alpha,beta)
    const double v;                // high quantile, same context



    //--------------------//
    // INITIALISATION METHODS
    void initialise(const double* dataArr);                   //initialise _data_ matrix
    void initialise();  //initialise _curr_



    //--------------------//
    // SWEEP UPDATES & RELATED FUNCTIONS
    void update_a(unsigned int const& iter);
    void update_b(unsigned int const& iter);
    double loglik_diff(double const& star, unsigned int const& dim, bool const& alpha);
    void update_mu();
    void update_sig();
    void update_ci();
    void update_comp();
    void update_w();
    void update_gam();
    void swap_1();       //adapted from Papaspiliopoulos & Roberts (2008, Bka)
    void swap_2();       //idem
    void swapcomp(unsigned int const& c1, unsigned int const& c2);//swap if swap_1 or swap_2 accepted



    //--------------------//
    // METHODS FOR BOUNDS ON (A,B)
    void bounds(bool fix_a, const double &val, double *bds, //finds the outer point closest to the true bound,
                unsigned int const& dim=0) const;           //according to _tol_
    bool cond(double const& a, double const& b,
              double const& p, unsigned int const& dim=0) const; //verifies if conditions are satisfied
    double qresid(double const& a, double const& b,
                  double const& p, unsigned int const& dim=0) const;//finds quantile p of residuals



    //--------------------//
    // END-OF-LOOP OPERATIONS & OTHERS
    void savetrace(unsigned int const& it); //store traces (checks for burn-in & thinning)
    void eol_msg(unsigned int const& it) const; //end-of-loop message, e.g., computing time up to _it_
    void rout(const char * msg, ...) const;
};

#endif // ETFIT_H
