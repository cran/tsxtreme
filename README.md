tsxtreme R package
==================

A New Approach to Fitting Time Series Extremal Dependence
---------------------------------------------------------

Classical approaches to fitting extremes of time series with short-term dependence use a pre-processing stage where independent extremes are filtered out of the original series. A well-known approach is the peaks-over-threshold method which typically involves

1. selecting a high threshold $u$,
2. choosing a minimal distance (*run-length*) beyond which excesses of $u$ are considered independent,
3. defining clusters of excesses based on this run-length,
4. making inference with the cluster maxima only.

This approach provides a simple procedure to estimating the marginal distribution of a series, but it lacks interpretation of the extremal dependence of the series; it is a multi-step procedure which does not recognise uncertainty in the definition of the clusters; the clustering procedure itself is arbitrary and can yield badly-biased estimates of risk measures; it discards all observations below $u$ and all exceedances of $u$ smaller than their cluster's maximum.

We model the sub-asymptotic dependence in stationary time series using the conditional approach to extreme values of Heffernan and Tawn (2004, JRSSB) formalised by Heffernan and Resnick (2007, Ann. Appl. Probab.). This approach has the advantage of covering a much broader class of extremal and asymptotic dependence structures than standard extreme value models, in particular situations where dependence at observed levels decays to asymptotic independence, defined as
$$\lim_{v\rightarrow 1}{\rm Pr}\{F_Y(Y) > v\mid F_X(X) > v\} = 0,\quad X\sim F_X,\ Y\sim F_Y,$$
as is often experienced in applications.

Example on Simulated Data
-------------------------

Generate some data from a stationary AR(1) and a Gaussian dependence structure.

```r
dep <- 0.7
n   <- 2e4
data<- numeric(n)
data[1] <- rnorm(1)
for(i in 2:n)
  data[i] <- rnorm(1, mean=dep*data[i-1], sd=1-dep^2)
```

The series in `data` has Gaussian marginal distribution; we can transform `data` to exponential scale so that the GPD marginal model is exact for any threshold.

```r
data <- qexp(pnorm(data))
u.gpd <- 0.9
```

We can now compare different estimators of the sub-asymptotic extremal index

$$\theta(x,m) = {\rm Pr}(X_1 < x,\ldots,X_m < x\mid X_0 > x), \quad x \geq u,$$

converging to the extremal index $\theta$ as $x,m\rightarrow\infty$ appropriately.

The *tsxtreme* package provides 3 different estimators of $\theta(x,m)$, namely

1. an empirical estimator based on the run-length $m$ with the routine `thetaruns()`,
2. an estimator that needs multi-step inference with `theta2fit()`,
3. an estimator based on our Bayesian semi-parametric approach with `thetafit()`.

A simple comparison of those 3 estimators can be done through

```r
x.empirical <- seq(u.gpd, 0.99, 0.005)
x.model     <- c(x.empirical, seq(0.995,0.9999,length.out=10))
theta.empirical <- thetaruns(ts = data, u.mar = u.gpd, probs = x.empirical)
theta.2step     <- theta2fit(ts = data, u.mar = u.gpd, u.dep = 0.98, probs = x.model)
theta.Bayesian  <- thetafit(ts = data, u.mar = u.gpd, u.dep = 0.98, probs = x.model)
par(mfrow=c(1,3))
plot(theta.empirical, main="Empirical")
plot(theta.2step, main="Stepwise Inference")
plot(theta.Bayesian, main="Semi-parametric Inference")
```

This can take some time, and we can change the default MCMC specifications of `bayesparams()`, e.g.

```r
my.par <- bayesparams(maxit=11000, burn=1000, adapt=1000, thin=5)
```

We also have access to prior hyperparameters, proposal standard deviations, verbosity, etc.

Package Structure
-----------------

The package can be divided in 3 parts, corresponding to the 3 estimators above, each providing different routines that have not necessarily their corresponding counterparts in the other approaches

* empirical
    + `thetaruns()` to compute the runs estimator for the extremal index $\theta$;
* stepwise
    + `dep2fit()` to fit the conditional tail model of Heffernan and Tawn using the standard stepwise approach;
    + `theta2fit()` to compute the sub-asymptotic extremal index $\theta(x,m)$; it can use the output of a previous call to `dep2fit()`;
* Bayesian
    + `depfit()` to fit the conditional tail model of Heffernan and Tawn using our Bayesian semi-parametric approach; additional structure can be imposed to fit first order Markov chains;
    + `thetafit()` to compute posterior samples of the sub-asymptotic extremal index $\theta(x,m)$; it can use the output of a previous call to `depfit()`;
    + `chifit()` to compute posterior samples of the sub-asymptotic measure of extremal $\chi_t(x)$ defined below; it can use the output of a previous call to `depfit()`.

`plot()` functions are available for viewing outputs from all routines listed above. The sub-asymptotic measure of extremal dependence at lag $t$ for a stationary time series $(X_t)$ is

$$\chi_t(x) = {\rm Pr}(X_t > x\mid X_0 > x).$$