## General Semiparametric Shared Frailty Model: Estimation and Simulation with frailtySurv

Authors: John V. Monaco, Malka Gorfine, Li Hsu

The script to generate results in the manuscript is [frailtySurv-jss.R](R/frailtySurv-jss.R). **Note**, this script is very computationally demanding: each section takes about a day. Some results are saved as RData files and then loaded during manuscript compilation. This includes the case studies and simulation results.

To run as a background script and log the results:

   $ nohup Rscript frailtySurv-jss.R > log 2> log.err < /dev/null

Data preprocessing scripts for the case studies are located in [code/python](python).

For more information, see [frailtySurv](https://github.com/vmonaco/frailtySurv).
