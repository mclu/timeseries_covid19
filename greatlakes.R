# Libraries
library(dplyr)
library(pomp)
library(ggplot2)
theme_set(theme_bw())
stopifnot(packageVersion("pomp")>="2.0")

data = read.csv("data.csv")
#data = read.csv("time_series_covid19_confirmed_global.csv")
#data = data %>% 
#  select(-c(Province.State, Lat, Long)) %>%
#  filter(Country.Region == "Italy") %>%
#  pivot_longer(cols = -1, names_to = "date", values_to = "C") %>%
#  transmute(date = as.Date(str_remove_all(date, "X"), format = "%m.%d.%y"),
#            C, day = c(1:95))

#80: ---------------------------------------------------------------------------
covid_statenames = c("S", "E", "I", "R")
covid_paramnames = c("Beta", "mu_EI", "rho", "mu_IR", "N", "eta", "k") 
covid_obsnames = "C"

covid_dmeasure = "lik = dpois(C, rho*I + 1e-6, give_log);"

covid_rmeasure = "C = rnbinom(rho*I, k);"

covid_rprocess = "
double dN_SE = rbinom(S, 1-exp(-Beta*I/N*dt));
double dN_EI = rbinom(E, 1-exp(-mu_EI*dt));
double dN_IR = rbinom(I, 1-exp(-mu_IR*dt));
S -= dN_SE;
E += dN_SE - dN_EI;
I += dN_EI - dN_IR;
R += dN_IR;
"

covid_rinit = "
S = 2500;
E = 1;
I = 1;
R = 0;
"

covid <- pomp(data = data, times = "day", t0 = 0,
              rprocess = euler(step.fun = Csnippet(covid_rprocess), delta.t = 1/7),
              rmeasure = Csnippet(covid_rmeasure),
              dmeasure = Csnippet(covid_dmeasure),
              partrans = parameter_trans(log=c("Beta","mu_EI","mu_IR", "k", "rho")),
              obsnames = covid_obsnames,
              statenames = covid_statenames,
              paramnames = covid_paramnames,
              rinit = Csnippet(covid_rinit)
)

# Parallel Setup
library(doParallel)
#cluster_cores <- as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE'))
#CLUSTER <- !is.na(cluster_cores)
#if(CLUSTER) registerDoParallel(cluster_cores) else registerDoParallel()
#library(doRNG)
#registerDoRNG(1929384)
cores = 36
registerDoParallel(cores)
mcopts = list(set.seed = TRUE)

# Run level
run_level <- 2
covid_Np <-      switch(run_level, 100, 10000, 60000)
covid_Nmif <-    switch(run_level,  10,   100,   300)
covid_Neval <-   switch(run_level,  10,    10,    10)
covid_Nglobal <- switch(run_level,  10,    10,   100)
covid_Nlocal <-  switch(run_level,  10,    100,    20)

covid_fixed_params = c(eta = 0.5, N = 3000)

covid_box <- rbind(
  Beta = c(10, 100),
  mu_EI = c(.01, .1),
  mu_IR = c(.01, .1), 
  k = c(.1, 1),
  rho = c(0.1, 5)
)

covid_rw.sd <- 0.02
covid_cooling.fraction.50 <- 0.5
stew(file=sprintf("box_eval-%d.rda", run_level),{
  t_local <- system.time({
    covid_local <- foreach(i=1:covid_Nlocal, .options.multicore=mcopts,
                           .packages='pomp', .combine=c) %dopar% mif2(
                             covid,
                             params=c(apply(covid_box,1,function(x) runif(1, x[1], x[2])),
                                      covid_fixed_params),
                             Nmif = covid_Nmif,
                             Np = covid_Np,
                             rw.sd=rw.sd(
                               Beta=covid_rw.sd,
                               mu_EI=covid_rw.sd,
                               mu_IR=covid_rw.sd,
                               k=covid_rw.sd,
                               rho=covid_rw.sd,
                               N = ivp(covid_rw.sd),
                               eta = ivp(covid_rw.sd)
                             ),
                             cooling.fraction.50 = covid_cooling.fraction.50
                           )
  })
},seed=900242058,kind="L'Ecuyer")


# likelihood local
stew(file=sprintf("lik_local_eval-%d.rda",run_level),{
  t_local_eval <- system.time({
    liks_local <- foreach(i=1:covid_Nlocal, .options.multicore=mcopts,
                          .packages='pomp',.combine=rbind) %dopar% {
                            evals <- replicate(covid_Neval,
                                               logLik(pfilter(covid,params=coef(covid_local[[i]]),Np=covid_Np)))
                            logmeanexp(evals, se=TRUE)
                          }
  })
})
results_local <- data.frame(logLik=liks_local[,1],
                            logLik_se=liks_local[,2],t(sapply(covid_local,coef)))

# global
# Global Search
stew(file=sprintf("g_box_eval-%d.rda", run_level),{
  t_global <- system.time({
    covid_global <- foreach(i=1:covid_Nglobal, .options.multicore=mcopts,
                           .packages='pomp', .combine=c) %dopar% mif2(
                             covid_local[[1]],
                             params=c(apply(covid_box,1,function(x) runif(1, x[1], x[2])),
                                      covid_fixed_params)
                           )
  })
},seed=900242058,kind="L'Ecuyer")

# likelihood eval at each point estimate
stew(file=sprintf("lik_global_eval-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:covid_Nglobal, .options.multicore=mcopts,
                           .packages='pomp',.combine=rbind) %dopar% {
                             evals <- replicate(covid_Neval,
                                                logLik(pfilter(covid,params=coef(covid_global[[i]]),Np=covid_Np)))
                             logmeanexp(evals, se=TRUE)
                           }
  })
},seed=442141592,kind="L'Ecuyer")

results_global <- data.frame(
  logLik=liks_global[,1], 
  logLik_se=liks_global[,2],
  t(sapply(covid_global,coef)))
