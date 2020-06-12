# POMP Modeling for COVID-19 in Italy

## Overview
The repo explores the transmission dynamics of COVID-19 in Italy using the partially observed Markov process model via the POMP package in R.

## Navigation
- [time_series_covid19_confirmed_global.csv](https://github.com/mclu/timeseries_covid19/raw/master/time_series_covid19_confirmed_global.csv): Daily new confirmed cases in Italy extracted from the website of Johns Hopkins Coronavirus Resource Center.
- [final.R](https://github.com/mclu/timeseries_covid19/blob/master/final.R): The script of analysis.
- [greatlakes.sbat](https://github.com/mclu/timeseries_covid19/blob/master/greatlakes.sbat): A script that used to submit on the cluster.
- [pf-1.rda](https://github.com/mclu/timeseries_covid19/blob/master/pf-1.rda), [box_eval-2.rda](https://github.com/mclu/timeseries_covid19/blob/master/box_eval-2.rda), [g_box_eval-2.rda](https://github.com/mclu/timeseries_covid19/blob/master/g_box_eval-2.rda), [lik_local_eval-2.rda](https://github.com/mclu/timeseries_covid19/blob/master/lik_local_eval-2.rda), [lik_global_eval-2.rda](https://github.com/mclu/timeseries_covid19/blob/master/lik_global_eval-2.rda): Files storing the simulation results.
- [final.Rmd](https://github.com/mclu/timeseries_covid19/blob/master/final.Rmd), [final.html](https://github.com/mclu/timeseries_covid19/blob/master/final.html): Files documenting the results.

## Installation
To run the R script, the following packages should be pre-installed in the IDE.
```r
install.packages('tidyverse')
install.packages('pomp')
install.packages('doParallel')
```

## Remark
The repo is the final project of the course Stats531 Analysis of Time Series in winter 2019.
