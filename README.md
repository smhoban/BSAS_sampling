# BSAS_sampling
All code for recreating simulations and analyses in paper 
**New guidance for ex situ gene conservation- sampling realistic population systems and accounting for collection attrition**

As explained in the paper Methods, the process consists of three phases.  First, the population system is simulated.  Second, the simulated system is sampled from at different levels of effort (number of plants sampled).  Third, the sample is compared to the entire system to determine what proportion of alleles (genetic variants) are actually present (and thus conserved) in the sample.  This calculation is done for several types of alleles and for one copy and multiple copies of each allele.  The following is instructions for each of these components of the work.  After this will be explanations of the code for creating the Figures and other calculations.

Step 1: Simulations 
Simulations are started using the file "do_sim_BSAS.R". This will loop over all parameter combinations, create the .par files for simcoal (using the function "simc.write"), and run the simulations (using the function "simc.run").  These write and run functions are in the files in the src folder "write.run.simcoal...".

This will create 1076 folders, each with 15 (or whatever number of reps you choose) simulation output files. I have then copied all these into one folder called Simulation6.

Step 2 and 3: Sampling and Comparison
Sampling and comparison are performed using the file "BSAS_sampling.R".  This will loop over all files that exist in the folder Simulations6.  This is well commented but briefly it will first convert all files from .arp to .gen format, calculate the number and type of alleles that exist in the system, then loop over all sampling efforts, and for each, create a genind object that represents a sample of a given effort, and then make the comparison of the sample to the total.  The data which is then recorded is the amount of genetic diversity captured by each sample effort.  Lastly a higher level of data is recorded using the gt95() function, this is the Nt or sampling effort required to capture 95% of the alleles in each allele category.  This is the result reported in the paper.
This step utilizes functions in the "sample_funcs_BSAS.R" file (in the main folder), as well as "arp2gen_edit.R" (in the src folder).  The latter file is a modification of the conversion script in adegenet, which corrects a small error in the original function.  I have created another function called "conv_arp_gen.R" but have not yet incorporated this.

Figure creation and calculations
