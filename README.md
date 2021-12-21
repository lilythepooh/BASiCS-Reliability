# Simulation-based-Evaluation-of-BHM-for-scRNA-Data
* To install modified BASiCS package, dowload BASiCS_2.2.1.tar.gz from this repository, in your R terminal, write "R CMD INSTALL BASiCS_2.2.1.tar.gz"
* To simulate data fron non-regression BASiCS, run 'simulate_data_code_epsilon=0.R'.
* To simulate data from regression BASiCS, run 'sim_data_BASiCS2018.R'.
* To run synthetic data with non-regression BASiCS, run 'sim_data_epsilon=0.R'
* To run synthetic data with regression BASiCS, run 'sim_NewBAsiCS_default.R'.
* To resimulate and rerun BASiCS MCMC for non-regression BASiCS for 100 replications, run 'resimulate_data_epsilon=0.R'.
* To resimulate and rerun BASiCS MCMC for regression BASiCS for 100 replications, run 'resimulate_NewBASiCS_default.R'.
* To run fixed synthetic data on modified non-regression BASiCS MCMC with epsilon=0,0.25,0.5,0.5,0.75,1, run 'sim_data_epsilon=0quarter.R', 'sim_data_epsilon=1quarter.R', 'sim_data_epsilon=2quarter.R','sim_data_epsilon=3quarter.R','sim_data_epsilon=4quarter.R'.
* To calculate effective sample size for SBC and rerun BASiCS MCMC for SBC when the effective smaple size is too small, run 'sbc_minNeff.R'
* To calculate rank statistics for SBC, run 'SBC_rank_statistics.R'
* To plot rank statistics in ecdf plots, run'rank_plot.R'
