#################################################
# Allele Frequency Shift and Number Alleles		#
#################################################

#This will look at possible explanations for results, in particular for bottleneck and population size

setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
library(adegenet); library(diveRsity)
source("sample_funcs_BSAS.R"); source("src/arp2gen_edit.R")

load("all_res_final.Rdata")
scenarios_ran<-list.dirs(path = "./Simulations6/", full.names = TRUE,recursive=F)
#set up the params list
get_params<- function(scenario) strsplit(strsplit(scenario,"/")[[1]][4],"_")
sim_params<-matrix(unlist(lapply(scenarios_ran,get_params)),ncol=7,byrow=T)
sim_params[,5]<-paste(".",sim_params[,5],sep="")
allele_categories<-c("all alleles (including rare)", "overall common f>0.05", "overall low fr 0.10>f>0.01", 
					 "overall rare f<0.01", "local1", "local2", "somwhere common f>0.05", "somwhere common f>0.05", "BM overall f~0.05")
	
results_allele_shift<-matrix(nrow=length(scenarios_ran),ncol=9)

#We'll do all scenarios and 3 alleles with 3 stats (# alleles, median, .95 quantile) each

for (scen in 1:length(scenarios_ran))
	{
	num_reps<-4
	concat_af<-list(0,0,0)
	for (nrep in 1:num_reps){
		reps_ran_gen<-list.files(scenarios_ran[[scen]], pattern="gen")
		temp_file_name<-file.path(scenarios_ran[[scen]],reps_ran_gen[nrep],sep="")
			if (.Platform$OS.type=="unix")  temp_file_name<-substr(temp_file_name,1,nchar(temp_file_name)-1)
		BSAS_genind<-read.genepop(temp_file_name,ncode=3)
		BSAS_genpop<-genind2genpop(BSAS_genind)
	
		n_total_indivs<- length(BSAS_genind@tab[,1])
		n_ind_p_pop<-table(BSAS_genind@pop)
		allele_freqs<-colSums(BSAS_genpop@tab)/(n_total_indivs*2)	
	
		allele_cat<-get.allele.cat(BSAS_genpop, n_ind_p_pop, local=T)	
		concat_af[[1]]<-c(concat_af[[1]],allele_freqs[allele_cat[[1]]])		#all alleles
		concat_af[[2]]<-c(concat_af[[2]],allele_freqs[allele_cat[[3]]])		#low frequency
		concat_af[[3]]<-c(concat_af[[3]],allele_freqs[allele_cat[[7]]])		#somewhere common
	}
results_allele_shift[scen,1:3]<-unlist(lapply(concat_af,length))
results_allele_shift[scen,4:6]<-unlist(lapply(concat_af,median))
results_allele_shift[scen,7:9]<-unlist(lapply(concat_af,quantile,.95))
}
write.csv(results_allele_shift,"results_allele_shift.csv")
	
	#explanation of columns of results table
	#first three columns are number of alleles for all, low freq, and some_com alleles
	#next three columns are median frequency for all, low freq, and some_com alleles
	#next three columns are the .95 quantile of allele frequency for all, low freq, and some_com alleles
	
	
results_allele_shift<-read.csv("results_allele_shift.csv")	
results_allele_shift<-results_allele_shift[,-1]
#########################
#	with bottleneck		#
#########################
get_params<- function(scenario) strsplit(strsplit(scenario,"/")[[1]][4],"_") 
sim_params<-matrix(unlist(lapply(scenarios_ran,get_params)),ncol=7,byrow=T)
for (i in 1:9){
	print(c(mean(results_allele_shift[sim_params[,1]=="b1LBSAS",i]),
	mean(results_allele_shift[sim_params[,1]=="b5LBSAS",i]),
	mean(results_allele_shift[sim_params[,1]=="b25LBSAS",i])))
	}
#for i=1, all alleles ... 157486/201183= 0.783
#Results for i=4 (all alleles, allele frequency)...	.0384/.032= 1.2

#################################
#	With population size		#
#################################
#(pop size is "100"  "1000" "150"  "200"  "300"  "50"   "500"  "75" ).. double checked this 1/4/2019
num_mark<-c(10000,	1500,	9000,	8000,	5000,	15000,	4000,	10000)
boxplot((results_allele_shift[,1]/num_mark)~sim_params[,7])
for (i in 1:3){
	print(mean(results_allele_shift[sim_params[,7]=="1000",i])/1500)
	print(mean(results_allele_shift[sim_params[,7]=="50",i])/15000)
	}
#Results for i=1 (all alleles, number of alleles per locus)... 21.01/26.96=0.779
for (i in 4:9){
	print(mean(results_allele_shift[sim_params[,7]=="1000",i]))
	print(mean(results_allele_shift[sim_params[,7]=="50",i]))
	}
#Results for i=4 (all alleles, allele frequency)... .0308/.0403=0.764



