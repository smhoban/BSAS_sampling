# BSAS_sampling
All code for recreating simulations and analyses in paper 
**New guidance for ex situ gene conservation- sampling realistic population systems and accounting for collection attrition**

Simulations are started using the file "do_sim_BSAS.R". This will loop over all parameter combinations, create the .par files for simcoal (using the function "simc.write"), and run the simulations (suing the function "simc.run").

This will create 1076 folders, each with 15 (or whatever number of reps you choose) simulation output files. I copied all these into one folder called Simulation6.

Everything else...
