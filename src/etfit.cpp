/*
 *  tsxtreme : Bayesian Modelling of Extremal Dependence in Time Series
 *  Copyright (C) 2017   Thomas Lugrin
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


//////////////////////////////////////////////////
// CONSTRUCTOR & DESTRUCTOR

ETfit::ETfit(double* data, int* n, int* nlag,
             int* k, int* kred, int* maxit,
             int *burn, int *thin, int *adapt, int *batchsize,
             double* sd_propa, double* sd_propb, double* mu_prior,
             double *nu_prior, double *eta_prior,
             const debmode &mode, const submodel &spec)
    : n(*n), nlag(*nlag), k(*k), kred(*kred),
      maxit(*maxit), burn(*burn), thin(*thin),
      maxadapt(*adapt),batchsize(*batchsize),
      V(*k,0), sumV(0),
      mode(mode),spec(spec),
      tol(0.001), v(-log(2*(1-0.99999999)))
{
    GetRNGstate();
    rout("DEBUG: Entering constructor...\n");
    nbswaps1=0; nbswaps2=0;

    mu[0] = mu_prior[0]; mu[1] = mu_prior[1];
    nu[0] = nu_prior[0]; nu[1] = nu_prior[1];
    eta[0] = eta_prior[0]; eta[1] = eta_prior[1];

    initialise(data);
    initialise();
    // regions for RAMA: 5 for alpha + 3 for beta = 8
    if(spec==none){
        curr.sd = std::vector<std::vector<double> >(this->nlag,std::vector<double>(8));
        for(unsigned int d=0;d<this->nlag;d++){
            for(unsigned int reg=0; reg<8; reg++){
                curr.sd[d][reg] = (reg<5)?sd_propa[reg]:sd_propb[reg-5];
            }
        }
    }else{//avoid bugs in structures when spec==univariateMixture
        curr.sd = std::vector<std::vector<double> >(1,std::vector<double>(8));
        for(unsigned int reg=0; reg<8; reg++){ curr.sd[0][reg] = (reg<5)?sd_propa[reg]:sd_propb[reg-5]; }
    }
    rout("DEBUG: initialisation finished...\n");

    savetrace(0);
    rout("DEBUG: v=%f\n",v);
}

ETfit::~ETfit(){
    PutRNGstate();
}

//////////////////////////////////////////////////
// INITIALISATION


void ETfit::initialise(const double *dataArr){
    std::vector<double> zeros(nlag+1,0);
    data.clear();

    for(unsigned int i=0; i<n; i++){
        for(unsigned int d=0; d<nlag+1; d++){
            if(d==0){ data.push_back(zeros); }
            data[i][d] = dataArr[n*d + i];
        }
    }
}

void ETfit::initialise(){
    curr.a.clear(); curr.b.clear();
    double bds[2];

    if(spec==univariateMixture){
        for(unsigned int d=0; d<nlag; d++){
            curr.a.push_back(0);
            curr.b.push_back(0);
        }
        acc_a.clear();acc_b.clear();
    }else if(spec==firstOrderMarkov){
        bounds(false, 0, bds); // get bounds ok for beta=0
        curr.a.push_back(runif(bds[0],bds[1]));
        bounds(true, curr.a[0], bds); // satisfied only for first lag
        curr.b.push_back(runif(bds[0], bds[1]));
        for(unsigned int d=1; d<nlag; d++){
            curr.a.push_back(pow(curr.a[0],d+1));
            curr.b.push_back(curr.b[0]);
        }
        // structure more complex than needed due to spec==none cases
        acc_a=std::vector<std::vector<std::vector<double> > >
                (1,std::vector<std::vector<double> >(5));//alpha: RAMA, regions [-1;-0.9;-0.1;0.1;0.9;1]
        acc_b=std::vector<std::vector<std::vector<double> > >
                (1,std::vector<std::vector<double> >(3));//beta: RAMA, regions [0;0.1;0.9;1]
    }else if(spec==none){
        for(unsigned int d=0; d<nlag; d++){
            bounds(false, 0, bds, d); // get bounds ok for beta=0
            curr.a.push_back(runif(bds[0],bds[1]));
            bounds(true, curr.a[d], bds, d);
            curr.b.push_back(runif(bds[0], bds[1]));
        }
        acc_a=std::vector<std::vector<std::vector<double> > >
                (nlag,std::vector<std::vector<double> >(5));//alpha: RAMA, regions [-1;-0.9;-0.1;0.1;0.9;1]
        acc_b=std::vector<std::vector<std::vector<double> > >
                (nlag,std::vector<std::vector<double> >(3));//beta: RAMA, regions [0;0.1;0.9;1]
    }else{
        error("case %u is not defined as a submodel type in ETfit::initialise(const submodel&)",spec);
    }

    curr.mu.clear(); curr.sig.clear();
    curr.w.resize(k);
    std::vector<double> zeros(nlag,0);

    for(unsigned int c=0; c<k; c++){
        for(unsigned int d=0; d<nlag; d++){
            if(d==0){
                curr.mu.push_back(zeros);//create a row
                curr.sig.push_back(zeros);
            }
            curr.mu[c][d]  = rnorm(mu[0], mu[1]);
            curr.sig[c][d] = sqrt(1/rgamma(nu[0], 1/nu[1]));
        }
        curr.w[c] = 1.0/k;
    }

    curr.gam = 1;

    curr.ci.resize(n);
    curr.noo.resize(k, 0);
    unsigned int ind;

    for(unsigned int i=0; i<n; i++){
        ind = rmult(curr.w, 1);
        curr.ci[i]   = ind;
        curr.noo[ind]++;
    }

    curr.noc = 0;
    for(unsigned int c=0; c<k; c++){
        if(curr.noo[c] > 0){ curr.noc++; }
    }
}



//////////////////////////////////////////////////
// SWEEP UPDATES

void ETfit::update_a(const unsigned int &iter){
    if(spec==univariateMixture){
        for(unsigned int d=0; d<nlag; d++){ curr.a[d]=0; }
    }else if(spec==firstOrderMarkov or spec==none){
        unsigned int const d_max = (spec==none)?nlag:1;
        double prop=0, prop_star=0;
        double a_star=0;
        double delta=0;

        // adapt (Regional Adaptive Metropolis Algorithm) [-1;-0.9] [-0.9;-0.1] [-0.1;0.1] [0.1;0.9] [0.9;1]
        if(iter<maxadapt and iter%batchsize==0){
            delta=std::min(0.05,5/(sqrt(iter)));
            for(unsigned int d=0; d<d_max; d++){
                for(unsigned int reg=0; reg<5; reg++){
                if(acc_a[d][reg].size()>0){
                    double weight=(double)acc_a[d][reg].size()/batchsize;
                    if(mean(acc_a[d][reg])>0.44) curr.sd[d][reg]*=exp(delta*weight);
                    else                         curr.sd[d][reg]/=exp(delta*weight);
                }else{
                    curr.sd[d][reg]*=exp(delta*0.01);//as if 1% of occurences with low acceptance
                }
                if(curr.sd[d][reg]<1e-3)    { curr.sd[d][reg]=1e-3;}
                else if(curr.sd[d][reg]>0.5){ curr.sd[d][reg]=0.5; }
                acc_a[d][reg].clear();
                }
            }
        }

        // normal proposal
        for(unsigned int d=0; d<d_max; d++){
            if(curr.a[d]<-0.9)                          a_star = rnorm(curr.a[d],curr.sd[d][0]);
            else if(curr.a[d]>=-0.9 and curr.a[d]<-0.1) a_star = rnorm(curr.a[d],curr.sd[d][1]);
            else if(curr.a[d]>=-0.1 and curr.a[d]<0.1)  a_star = rnorm(curr.a[d],curr.sd[d][2]);
            else if(curr.a[d]>=0.1 and curr.a[d]<0.9)   a_star = rnorm(curr.a[d],curr.sd[d][3]);
            else                                        a_star = rnorm(curr.a[d],curr.sd[d][4]);

            if(cond(a_star,curr.b[d],0,d) and cond(a_star,curr.b[d],1,d)){
                // difference of likelihoods
                double like_diff = loglik_diff(a_star,d,true);
                //proposal densities
                if(a_star<-0.9)                       prop = dnorm(curr.a[d],a_star,curr.sd[d][0],1);
                else if(a_star>=-0.9 and a_star<-0.1) prop = dnorm(curr.a[d],a_star,curr.sd[d][1],1);
                else if(a_star>=-0.1 and a_star<0.1)  prop = dnorm(curr.a[d],a_star,curr.sd[d][2],1);
                else if(a_star>=0.1 and a_star<0.9)   prop = dnorm(curr.a[d],a_star,curr.sd[d][3],1);
                else                                  prop = dnorm(curr.a[d],a_star,curr.sd[d][4],1);

                if(curr.a[0]<-0.9)                          prop_star = dnorm(a_star,curr.a[d],curr.sd[d][0],1);
                else if(curr.a[0]>=-0.9 and curr.a[0]<-0.1) prop_star = dnorm(a_star,curr.a[d],curr.sd[d][1],1);
                else if(curr.a[0]>=-0.1 and curr.a[0]<0.1)  prop_star = dnorm(a_star,curr.a[d],curr.sd[d][2],1);
                else if(curr.a[0]>=0.1 and curr.a[0]<0.9)   prop_star = dnorm(a_star,curr.a[d],curr.sd[d][3],1);
                else                                        prop_star = dnorm(a_star,curr.a[d],curr.sd[d][4],1);

                // acceptance (cst prior)
                double accept = std::min(1.0,exp(like_diff - prop_star + prop)) ;
                if(runif(0,1)<=accept){
                    if(spec==firstOrderMarkov){
                        for(unsigned int j=0; j<nlag; j++){ curr.a[j]=pow(a_star,j+1); }
                    }else{//if spec==none
                        curr.a[d]=a_star;
                    }
                }
                if(iter<maxadapt+batchsize*10){//avoid memory fill
                    if(a_star<-0.9)                       acc_a[d][0].push_back(accept);
                    else if(a_star>=-0.9 and a_star<-0.1) acc_a[d][1].push_back(accept);
                    else if(a_star>=-0.1 and a_star<0.1)  acc_a[d][2].push_back(accept);
                    else if(a_star>=0.1 and a_star<0.9)   acc_a[d][3].push_back(accept);
                    else                                  acc_a[d][4].push_back(accept);
                }
            }
        }
    }else{
        error("case %u is not defined as a submodel type in ETfit::update_a(...)",spec);
    }
}


void ETfit::update_b(const unsigned int &iter){
    if(spec==univariateMixture){
        for(unsigned int d=0; d<nlag; d++){ curr.b[d]=0; }
    }else if(spec==firstOrderMarkov or spec==none){
        unsigned int const d_max = (spec==none)?nlag:1;
        double prop=0, prop_star=0;
        double b_star=0;
        double delta=0;

        // adapt (Regional Adaptive Metropolis Algorithm) [0;0.1] [0.1;0.9] [0.9;1]
        if(iter<maxadapt and iter%batchsize==0){
            delta=std::min(0.05,5/(sqrt(iter)));
            for(unsigned int d=0; d<d_max; d++){
                for(unsigned int reg=0; reg<3; reg++){
                    if(acc_b[d][reg].size()>0){
                        double weight = (double)acc_b[d][reg].size()/batchsize;
                        if(mean(acc_b[d][reg])>0.44) curr.sd[d][reg+5]*=exp(delta*weight);
                        else                         curr.sd[d][reg+5]/=exp(delta*weight);
                    }else{
                        curr.sd[d][reg+5]*=exp(delta*0.01);//as if one occurence with low acceptance
                    }
                    if(curr.sd[d][reg+5]<1e-3)    { curr.sd[d][reg+5]=1e-3;}
                    else if(curr.sd[d][reg+5]>0.5){ curr.sd[d][reg+5]=0.5; }
                    acc_b[d][reg].clear();
                }
            }
        }

        // normal proposal
        for(unsigned int d=0; d<d_max; d++){
            if(curr.b[d]<0.1)      b_star = rnorm(curr.b[d],curr.sd[d][5]);
            else if(curr.b[d]>0.9) b_star = rnorm(curr.b[d],curr.sd[d][7]);
            else                   b_star = rnorm(curr.b[d],curr.sd[d][6]);

            if(cond(curr.a[d],b_star,0,d) and cond(curr.a[d],b_star,1,d)){
                // difference of likelihoods
                double like_diff = loglik_diff(b_star,d,false);
                // proposal densities
                if(b_star<0.1)      prop = dnorm(curr.b[d],b_star,curr.sd[d][5],1);
                else if(b_star>0.9) prop = dnorm(curr.b[d],b_star,curr.sd[d][7],1);
                else                prop = dnorm(curr.b[d],b_star,curr.sd[d][6],1);

                if(curr.b[d]<0.1)      prop_star = dnorm(b_star,curr.b[d],curr.sd[d][5],1);
                else if(curr.b[d]>0.9) prop_star = dnorm(b_star,curr.b[d],curr.sd[d][7],1);
                else                   prop_star = dnorm(b_star,curr.b[d],curr.sd[d][6],1);

                // acceptance (cst prior)
                double accept = std::min(1.0,exp(like_diff - prop_star + prop));
                if(runif(0,1)<=accept){
                    if(spec==firstOrderMarkov){
                        for(unsigned int j=0; j<nlag; j++){ curr.b[j]=b_star; }
                    }else{//if spec==none
                        curr.b[d]=b_star;
                    }
                }
                if(iter<maxadapt+batchsize*10){//avoid memory fill
                    if(b_star<0.1)      acc_b[d][0].push_back(accept);
                    else if(b_star>0.9) acc_b[d][2].push_back(accept);
                    else                acc_b[d][1].push_back(accept);
                }
            }
        }
    }else{
        error("case %u is not defined as a submodel type in ETfit::update_b(...)",spec);
    }
}


double ETfit::loglik_diff(const double &star, const unsigned int &dim, const bool &alpha){
    double like=0,like_star=0;
    unsigned int d_min=(spec==none)?dim:0;//if Markov chain, then likelihood computed in one go
    unsigned int const d_max=(spec==none)?(dim+1):nlag;// otherwise: each dimension separately
    double a_star;

    for(unsigned int d=d_min; d<d_max; d++){
        for(unsigned int i=0; i<n; i++){
            like += dnorm(data[i][d+1],
                     curr.a[d]*data[i][0] + pow(data[i][0], curr.b[d])*curr.mu[curr.ci[i]][d],
                     pow(data[i][0], curr.b[d])*curr.sig[curr.ci[i]][d], 1);
            if(alpha){
                if(spec==none){ a_star = star; }
                else          { a_star = pow(star, d+1); }
                like_star += dnorm(data[i][d+1],
                        a_star*data[i][0] + pow(data[i][0], curr.b[d])*curr.mu[curr.ci[i]][d],
                        pow(data[i][0], curr.b[d])*curr.sig[curr.ci[i]][d], 1);
            }else{
                like_star += dnorm(data[i][d+1],
                        curr.a[d]*data[i][0] + pow(data[i][0], star)*curr.mu[curr.ci[i]][d],
                        pow(data[i][0], star)*curr.sig[curr.ci[i]][d], 1);
            }
        }
    }
    return(like_star-like);
}

void ETfit::update_mu(){
    rout("DEBUG: entering update_mu()...\n");
    double s_star=0, m_star=0; //_s_star_ is a variance, not a standard deviation

    for(unsigned int c=0; c<k; c++){
        if(curr.noo[c]==0){
            for(unsigned int d=0; d<nlag; d++){
                curr.mu[c][d] = rnorm(mu[0], mu[1]);
            }
        }else{
            for(unsigned int d=0; d<nlag; d++){
                s_star = curr.noo[c]/pow(curr.sig[c][d],2) + 1/pow(mu[1],2);
                s_star = 1/s_star;
                m_star = 0;
                for(unsigned int i=0; i<n; i++){
                    if(curr.ci[i]==c){
                        m_star += (data[i][d+1]-curr.a[d]*data[i][0])/pow(data[i][0], curr.b[d]);
                    }
                }
                m_star = s_star * (m_star/pow(curr.sig[c][d],2) + mu[0]/pow(mu[1],2));
                curr.mu[c][d] = rnorm(m_star, sqrt(s_star));
            }
        }
    }
}


void ETfit::update_sig(){
    rout("DEBUG: entering update_sig()...\n");
    double nu2_star=0, m_star=0;

    for(unsigned int c=0; c<k; c++){
        if(curr.noo[c]==0){
            for(unsigned int d=0; d<nlag; d++){
                curr.sig[c][d] = sqrt(1/rgamma(nu[0], 1/nu[1]));
            }
        }else{
            for(unsigned int d=0; d<nlag; d++){
                nu2_star = 0;
                for(unsigned int i=0; i<n; i++){
                    if(curr.ci[i]==c){
                        m_star    = curr.a[d]*data[i][0] + pow(data[i][0],curr.b[d])*curr.mu[c][d];
                        nu2_star += pow(data[i][d+1]-m_star, 2) / pow(data[i][0], 2*curr.b[d]);
                    }
                }
                nu2_star /= 2;
                nu2_star += nu[1];
                curr.sig[c][d] = sqrt(1/rgamma(nu[0] + curr.noo[c]/2.0, 1/nu2_star));
            }
        }
    }
}


void ETfit::update_ci(){
    rout("DEBUG: entering update_ci()...\n");
    double sum=0;
    double m_star=0, s_star=0;
    std::vector<double> w(k,0);

    for(unsigned int i=0; i<n; i++){
        sum  = 0;
        for(unsigned int c=0; c<k; c++){
            w[c] = log(curr.w[c]);
            for(unsigned int d=0; d<nlag; d++){
                m_star = curr.a[d]*data[i][0] + pow(data[i][0], curr.b[d]) * curr.mu[c][d];
                s_star = pow(data[i][0],2*curr.b[d]) * pow(curr.sig[c][d],2);
                w[c]  -= pow(data[i][d+1]-m_star,2) / (2*s_star) + log(s_star)/2;
            }
            w[c] = exp(w[c]);
            sum += w[c];
        }
        curr.ci[i] = rmult(w, sum);
    }
}


void ETfit::update_comp(){
    rout("DEBUG: entering update_comp()...\n");
    curr.noc = 0;
    for(unsigned int c=0; c<k; c++)
        curr.noo[c] = 0;

    for(unsigned int i=0; i<n; i++)
        curr.noo[curr.ci[i]]++;

    for(unsigned int c=0; c<k; c++)
        if(curr.noo[c]>0) curr.noc++;
}


void ETfit::update_w(){
    rout("DEBUG: entering update_w()...\n");
    double par=curr.gam + n;
    double prod=1;  // prod_{k=1}^{c-1} (1-V_k)
    sumV = 0;


    for(unsigned int c=0; c<k-1; c++){
        par  -= curr.noo[c];
        V[c]  = rbeta(1+curr.noo[c], par);
        sumV += log(1-V[c]);
        curr.w[c] = V[c]*prod;
        prod     *= 1-V[c];
    }
    rout("DEBUG: sumV = %f...\n",sumV);
    curr.w[k-1] = prod;
    V[k-1] = 1;
}


void ETfit::update_gam(){
    rout("DEBUG: entering update_gam()...\n");
    //curr.gam = rgamma(k+eta[0]-1, eta[1]/(1-eta[1]*sumV));
    // truncated gamma prior
    unsigned int index=0;
    do{
        curr.gam = rgamma(k+eta[0]-1, eta[1]/(1-eta[1]*sumV));
        index++;
    }while(curr.gam<0.5 and index<10000);
    if(index>=10000) curr.gam = 0.5;
}



//////////////////////////////////////////////////
// SWAPS: LABEL-SWITCHING


void ETfit::swap_1(){
    rout("DEBUG: entering swap_1()...\n");

    std::vector<double> w_ne(curr.noc,1.0/curr.noc);            //sub-vectors of non-empty components
    std::vector<unsigned int> ind_ne(curr.noc,0);
    unsigned int c1=0,c2=0;                          //indices of swap candidates
    double sum=1;                                    //sum of weights of non-empty components

    for(unsigned int i=0, c=0; c<k and i<curr.noc; c++){
        if(curr.noo[c]>0){
            ind_ne[i] = c;
            i++;
        }
    }
    //choose candidates
    c1  = rmult(w_ne, sum);
    sum-= w_ne[c1]; w_ne[c1] = 0;
    c2  = rmult(w_ne, sum);
    //compute log-acceptance ratio
    double lograt=0;
    double Wdown=0,Wup=0;

    c1 = ind_ne[c1];
    c2 = ind_ne[c2];
    if(c1 > c2) std::swap(c1,c2);

    if(c2==k-1){
        lograt = (curr.noo[c2]+curr.gam-1) * log(curr.w[c1]);
        lograt-= (curr.noo[c2]+curr.gam-1) * log(curr.w[c2]);
    }else{
        lograt = curr.noo[c2] * log(curr.w[c1]);
        lograt-= curr.noo[c2] * log(curr.w[c2]);
    }
    lograt += curr.noo[c1] * log(curr.w[c2]);
    lograt -= curr.noo[c1] * log(curr.w[c1]);

    for(unsigned int c=c1+1; c<k; c++) Wup += curr.w[c];
    Wdown = Wup-curr.w[c2]+curr.w[c1];

    for(unsigned int c=c1+1; c<std::min(c2+1,k-1); c++){
        lograt += log(Wup);
        lograt -= log(Wdown);
        Wup   -= curr.w[c];
        Wdown -= curr.w[c];
    }
    //accept-reject
    if(runif(0,1)<exp(lograt)){
        swapcomp(c1,c2);
        nbswaps1++;
    }
}

void ETfit::swap_2(){
    rout("DEBUG: entering swap_2()...\n");
    double lograt=0;        //log acceptance ratio
    double Wup=0,Wdown=0;   //prod_(j=2)^k W(tilde)_j, with W(tilde)_j=sum_(i=j)^K w_i

    std::vector<double> w_ne(curr.noc,1.0/curr.noc);  //sub-vectors of non-empty components
    std::vector<unsigned int> ind_ne(curr.noc,0);
    unsigned int c1=0;                                //index of swap candidate (the other being c1+1)
    double sum=1;                                     //sum of weights of non-empty components
    //select non-empty components
    for(unsigned int i=0, c=0; c<k and i<curr.noc; c++){
        if(curr.noo[c]>0){
            ind_ne[i] = c;
            i++;
        }
    }
    //choose candidates
    if(ind_ne[0] == 0){//ensure c+1>1
        sum    -= w_ne[0];
        w_ne[0] = 0;
    }
    c1 = rmult(w_ne, sum);
    c1 = ind_ne[c1]-1;//could be empty
    //swapped sticks
    double w1 = V[c1+1]/V[c1]                         * curr.w[c1];
    double w2 = V[c1]*(1-V[c1+1])/(V[c1+1]*(1-V[c1])) * curr.w[c1+1];
    //compute log acceptance ratio
    lograt = curr.noo[c1]*log(1-V[c1+1]);
    if(c1+1<k-1){
        lograt -= curr.noo[c1+1]*log(1-V[c1]);
    }else{
        lograt -= (curr.noo[c1+1]+curr.gam-1)*log(1-V[c1]);
    }
    for(unsigned int c=1; c<k; c++){ Wup += curr.w[c]; }   // initialise: W_2
    Wdown = Wup - curr.w[c1] - curr.w[c1+1] + w1 + w2;     // W(tilde)_2
    for(unsigned int c=1; c<std::min(c1+1,k-1); c++){
        lograt -= log(Wdown);
        lograt += log(Wup);
        Wdown -= curr.w[c];
        Wup   -= curr.w[c];
    }
    if(c1+1 < k-1){
        Wdown = Wup - curr.w[c1+1] + w2;
        lograt -= log(Wdown);
        lograt += log(Wup);
    }
    //accept-reject
    if(runif(0,1)<exp(lograt)){
        rout("DEBUG: swap 2 accepted for c=%u...\n",c1);
        curr.w[c1]   = w1;
        curr.w[c1+1] = w2;
        swapcomp(c1,c1+1);
        nbswaps2++;
    }
}


void ETfit::swapcomp(const unsigned int &c1, const unsigned int &c2){
    // swap components' sizes
    std::swap(curr.noo[c1],curr.noo[c2]);
    // swap components' indices
    for(unsigned int i=0; i<n; i++){
        if(curr.ci[i] == c1){
            curr.ci[i] = c2;
        }else if(curr.ci[i] == c2){
            curr.ci[i] = c1;
        }
    }
    // swap components' means and standard deviances
    std::swap(curr.mu[c1],curr.mu[c2]);
    std::swap(curr.sig[c1],curr.sig[c2]);
    // swap components' weights -> NO!
//    std::swap(curr.w[c1],curr.w[c2]);
}



//////////////////////////////////////////////////
// CORE: RUN ALGO


void ETfit::run(const algotype &type){
    rout("DEBUG: entering run()...\n");

    if(type==conditional){
        for(unsigned int it=1; it<maxit; it++){// iteration 0: initialisation in constructor
            rout("DEBUG: beginning sweep %u...\n###################################\n",it);
            update_a(it);
            update_b(it);
            update_mu();
            update_sig();
            update_ci();
            update_comp();
            update_w();
            update_gam();
            if(curr.noc>1){
                swap_1();
                swap_2();
            }

            savetrace(it);
            eol_msg(it);
        }
    }else if(type==marginal){
        error("in ETfit::run(): Marginal method of R. M. Neal not yet implemented...");
    }else{
        error("in ETfit::run(): algotype can be conditional or marginal only");
    }
}

//////////////////////////////////////////////////
// END-OF-LOOP OPERATIONS & OTHERS


void ETfit::savetrace(const unsigned int &it){
    // save only what matters!
    if(it+1 > burn and (it+1-burn)%thin == 0){
        traces.push_back(curr);
        traces.back().mu.resize(kred);
        traces.back().sig.resize(kred);
        traces.back().w.resize(kred);
        traces.back().noo.resize(kred);
    }
    R_CheckUserInterrupt();//receives interruptions from R console
}

void ETfit::eol_msg(const unsigned int &it) const{
    if((it+1)%2000==0){
        rout("Sweep %u reached...\n",it+1);
        if(spec!=univariateMixture){
            rout("         (a) [-1;-0.9] | [-0.9;-0.1] | [-0.1;0.1] | [0.1;0.9] | [0.9;1]\n");
            rout("Acceptance:  %.2f%%  |  %.2f%%   |  %.2f%%   |  %.2f%% |  %.2f%%\n",
                 (acc_a[0][0].size()>0)?mean(acc_a[0][0])*100:0,
                 (acc_a[0][1].size()>0)?mean(acc_a[0][1])*100:0,
                 (acc_a[0][2].size()>0)?mean(acc_a[0][2])*100:0,
                 (acc_a[0][3].size()>0)?mean(acc_a[0][3])*100:0,
                 (acc_a[0][4].size()>0)?mean(acc_a[0][4])*100:0);
            rout("         (b) [0;0.1]  [0.1;0.9]  [0.9;1]\n");
            rout("Acceptance:  %.2f%%  |  %.2f%%   |  %.2f%%\n",
                 (acc_b[0][0].size()>0)?mean(acc_b[0][0])*100:0,
                 (acc_b[0][1].size()>0)?mean(acc_b[0][1])*100:0,
                 (acc_b[0][2].size()>0)?mean(acc_b[0][2])*100:0);
            rout("Adapted sd: %.2e | %.2e | %.2e | %.2e | %.2e || %.2e | %.2e | %.2e\n",
                 curr.sd[0][0],curr.sd[0][1],curr.sd[0][2],curr.sd[0][3],curr.sd[0][4],curr.sd[0][5],curr.sd[0][6],curr.sd[0][7]);
        }
    }
    if(it+1==maxit)
        rout("\nnbr of swaps (1): %u, and (2): %u\n",nbswaps1,nbswaps2);
}

void ETfit::rout(const char *msg,...) const{
    if(mode!=silent){
        unsigned int i=0;
        const char* deb="DEBUG";

        while(*msg!='\0' and i<5){
            if(*msg!=deb[i]) break;//not a debug message => i<5
            i++;msg++;
        }
        if(i<4 or mode==debug){
            msg-=i;
            va_list args;
            va_start(args,msg);
            Rvprintf(msg,args);
            va_end(args);
        }
    }
}



//////////////////////////////////////////////////
// BOUNDS ON (ALPHA,BETA)

void ETfit::bounds(bool fix_a, double const& val, double* bds, const unsigned int &dim) const{
    double in=2,out=2;
    double mid;

    //find bounds for beta (verified for quantiles with p=0 and p=1)...
    if(fix_a){
        //initialise on bounds for beta : {0,1}
        if(cond(val, 0, 1, dim) and cond(val, 0, 0, dim)){ in=0; }
        else{ out=0; }
        if(cond(val, 1, 1, dim) and cond(val, 1, 0, dim)){ in=1; }
        else{ out=1; }
        if(out==2){ // no restrictions due to conditions
            bds[0]=0; bds[1]=1; return;
        }
        else if(in==2){ // restrictions on both sides
            rout("DEBUG: restrictions on both sides for fixed a...\n");
            const double eps=0.01;
            double testb=0;
            do{
                testb+=eps;
            }while((!cond(val,testb,1,dim) or !cond(val,testb,0,dim)) and testb<1);
            if(testb<1) in=testb;
            else error("Refine grid to get proper results in bounds() for fixed a...\n");
            out=0;
            //dichotomy
            do{
                mid = (in+out)/2;
                if(cond(val,mid,1,dim) and cond(val,mid,0,dim)){ in=mid; }
                else{ out=mid; }
                rout("DEBUG: fixed a: in=%f, out=%f\n",in,out);
            }while(fabs(in-out)>tol);
            bds[0] = in;//ensure to sample within the true bounds
            out=1;in=testb;
        }
        else{
            bds[0] = in;
        }
        //dichotomy
        do{
            mid = (in+out)/2;
            if(cond(val, mid, 1, dim) and cond(val, mid, 0, dim)){ in=mid; }
            else{ out=mid; }
            rout("DEBUG: in=%f, out=%f\n",in,out);
        }while(fabs(in-out)>tol);
        bds[1] = in;//ensure to sample within the true bounds
    }

    //...or on alpha
    else{
        //initialise on bounds for alpha : {-1,1}
        if(cond(-1, val, 1, dim) and cond(-1, val, 0, dim)){ in=-1; }
        else{ out=-1; }
        if(cond(1, val, 1, dim) and cond(1, val, 0, dim)){ in=1; }
        else{ out=1; }
        if(out==2){// no restrictions due to conditions
            bds[0]=-1; bds[1]=1; return;
        }
        else if(in==2){ // restrictions on both sides
            rout("DEBUG: restrictions on both sides for fixed b...\n");
            const double eps=0.01;
            double testa=-1;
            do{
                testa+=eps;
            }while((!cond(testa,val,1,dim) or !cond(testa,val,0,dim)) and testa<1);
            if(testa<1) in=testa;
            else error("Refine grid to get proper results in bounds() for fixed b...\n");
            out=-1;
            //dichotomy
            do{
                mid = (in+out)/2;
                if(cond(mid, val, 1, dim) and cond(mid, val, 0, dim)){ in=mid; }
                else{ out=mid; }
                rout("DEBUG: fixed b: in=%f, out=%f\n",in,out);
            }while(fabs(in-out)>tol);
            bds[0] = in;//ensure to sample within the true bounds
            out=1;in=testa;
        }
        else{
            bds[0] = in;
        }
        //dichotomy
        do{
            mid = (in+out)/2;
            if(cond(mid, val, 1, dim) and cond(mid, val, 0, dim)){ in=mid; }
            else{ out=mid; }
            rout("DEBUG: in=%f, out=%f\n",in,out);
        }while(fabs(in-out)>tol);
        bds[1] = in;//ensure to sample within the true bounds
    }
    //ensure bounds are ordered
    if(bds[1] < bds[0]){ std::swap(bds[0],bds[1]); }
}


bool ETfit::cond(double const& a, double const& b, double const& p, const unsigned int &dim) const{
    double z_pos = qresid(1, 0, p, dim);
    double z     = qresid(a, b, p, dim);
    double z_neg = qresid(-1, 0, p, dim);
    rout("DEBUG: a=%.1f, b=%.1f, z_pos=%f, z=%f, z_neg=%f,v=%f\n",a,b,z_pos,z,z_neg,v);
    if(a<-1 or a>1 or b>=1 or b<0){ return(false); }

    // Case I
    if(a > fmin2(1 - b*z*pow(v,b-1),  1 - pow(v,b-1)*z + z_pos/v)){//Case I first part not satisfied
        bool cond1 = (1 - b*z*pow(v,b-1) < a);
        bool cond2 = ( (1-1/b) * pow(b*z,1/(1-b)) * pow(1-a,-b/(1-b)) + z_pos > 0 );
        if(!cond1 or !cond2){ return(false); }//Case I second part not satisfied
    }
    // Case II
    if(-a > fmin2(1 + b*z*pow(v,b-1), 1 + pow(v,b-1)*z - z_neg/v)){//Case II first part not satisfied
        bool cond1 = ( 1 + b*pow(v,b-1)*z < -a );
        bool cond2 = ( (1-1/b) * pow(-b*z,1/(1-b)) * pow(1+a,-b/(1-b)) - z_neg > 0 );
        if(!cond1 or !cond2){ return(false); }//Case II second part not satisfied
    }
    return(true);
}


double ETfit::qresid(const double &a, const double &b, const double &p, const unsigned int &dim) const{
    double Z;//residuals
    double ret=0;

    for(unsigned int i=0; i<n; i++){
        Z = (data[i][dim+1] - a*data[i][0])/pow(data[i][0],b);
        if(i==0){ ret = Z; }
        else{
            if(p==0){ ret = fmin2(ret, Z); }
            else if(p==1){ ret = fmax2(ret, Z); }
            else{ error("only p=0 or 1 implemented in ETfit::qresid"); }
        }
    }
    return(ret);
}



//////////////////////////////////////////////////
// SIMULATE RANDOM VARIABLES

int ETfit::rmult(const std::vector<double> &probs, const double &sum){
    double u = runif(0,sum);
    double p = 0;

    for(unsigned int c=0; c<probs.size(); c++){
        p+=probs[c];
        if(u <= p){
            return(c);
        }
    }
    error("in rmult() (cluster assignment): u = %f, sum = %f, p = %f", u, sum, p);
}



//////////////////////////////////////////////////
// COMPUTE STATISTICS

double ETfit::mean(const std::vector<double> &x) const{
    if(x.size()<1) error("Empty vectors not supported in ETfit::mean.");
    double ret=0;

    for(unsigned int i=0; i<x.size(); i++){ ret+= x[i]; }
    return(ret/x.size());
}

double ETfit::var(const std::vector<double> &x) const{
    if(x.size()<2) error("Empty vectors or singletons not supported in ETfit::var.");
    double ret=0;
    double m=mean(x);

    for(unsigned int i=0; i<x.size(); i++){
        ret+=pow(x[i]-m,2);
    }
    return(ret/(x.size()-1));
}

double ETfit::cov(const std::vector<double> &x, const std::vector<double> &y) const{
    if(x.size()!=y.size()) error("Sizes of vectors in ETfit::cov do not match.");
    if(x.size()<2) error("Empty vectors or singletons not supported in ETfit::cov.");
    double ret=0;
    double mx=mean(x), my=mean(y);
    rout("DEBUG: [ETfit::cov] mean(x)=%.3f, mean(y)=%.3f\n", mx, my);

    for(unsigned int i=0; i<x.size(); i++){
        ret+=(x[i]-mx)*(y[i]-my);
    }
    return(ret/(x.size()-1));
}



//////////////////////////////////////////////////
// ACCESSORS

const std::vector<ETpar>& ETfit::getTraces() const{
    return(traces);
}
