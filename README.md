# Tuas-South-Dorm

This folder includes the code for the Tuas South Dorm project (Estimating transmission dynamics at different intraspatial levels in an institutional outbreak).

Model 1 is the model in which we did not account for impacts of the interventions while model 2 is the model in which we considered the impacts of the interventions.

For each model, there are 3 files: "mcmc" is the file which we used to obtain the posterior samples; "functions" and "functions_rcpp" are the files for the R and cpp functions we used to obtain the likelihoods for Metropolis Hasting sampling. 

"infection_simulation.r" is the file for simulating number of infections when we changed number of blocks or floors in the dormitory area, as well as visualizing the simulation results. 

"plot.r" is the file for visualizing the datasets and the outputs.
