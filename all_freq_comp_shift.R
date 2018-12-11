setwd("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")

library(adegenet); library(diveRsity)
library(parallel); library(doParallel) #will load foreach
source("sample_funcs_BSAS.R"); source("src/arp2gen_edit.R")
colMax <- function(data) sapply(data, max, na.rm = TRUE)

load("all_res_final.Rdata")

scenarios_ran<-list.dirs(path = "./Simulations6/", full.names = TRUE,recursive=F)
#set up the params list
get_params<- function(scenario) strsplit(strsplit(scenario,"/")[[1]][4],"_")
sim_params<-matrix(unlist(lapply(scenarios_ran,get_params)),ncol=7,byrow=T)
sim_params[,5]<-paste(".",sim_params[,5],sep="")
allele_categories<-c("all alleles (including rare)", "overall common f>0.05", "overall low fr 0.10>f>0.01", 
					 "overall rare f<0.01", "local1", "local2", "somwhere common f>0.05", "somwhere common f>0.05", "BM overall f~0.05",
					"all alleles (including rare)", "overall common f>0.05", "overall low fr 0.10>f>0.01", 
					 "overall rare f<0.01", "local1", "local2", "somwhere common f>0.05", "somwhere common f>0.05", "BM overall f~0.05")

	#50,100,150,300,1000 for 10 populations 1 bneck, 3 populations 1 bneck, and 10 populations 25 bneck
	
results_allele_shift<-matrix(nrow=length(scenarios_ran),ncol=9)
	#We'll do 15 scenarios and 3 alleles with 3 stats (median, .9, # alleles) each
	#sets_to_do<-c(6,1,3,5,2,230,225,227,229,226,398,393,395,397,394)

for (scen in 1:length(scenarios_ran))	#for (scen in c(343,344,338,339,340,341,342,337))	for (scen in c(247,248,242,243,244,245,246,241))
	{
	num_reps<-4
	concat_af<-list(0,0,0)
	for (nrep in 1:num_reps){
		reps_ran_gen<-list.files(scenarios_ran[[scen]], pattern="gen")
		temp_file_name<-file.path(scenarios_ran[[scen]],reps_ran_gen[nrep],sep="")
			if (.Platform$OS.type=="unix")  temp_file_name<-substr(temp_file_name,1,nchar(temp_file_name)-1)
		BSAS_genind<-read.genepop(temp_file_name,ncode=3)
		BSAS_genpop<-genind2genpop(BSAS_genind)
	
		#--NUMBER OF POPULATIONS, INDIVIDUALS, REGIONS, REGIONAL MAKEUP--#
		n_pops<-length(levels(BSAS_genind@pop))
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
	
#with bottleneck
get_params<- function(scenario) strsplit(strsplit(scenario,"/")[[1]][4],"_") 
sim_params<-matrix(unlist(lapply(scenarios_ran,get_params)),ncol=7,byrow=T)
for (i in 1:9){
	print(c(mean(results_allele_shift[sim_params[,1]=="b1LBSAS",i]),
	mean(results_allele_shift[sim_params[,1]=="b5LBSAS",i]),
	mean(results_allele_shift[sim_params[,1]=="b25LBSAS",i])))
	}

#With population size
num_mark<-c(10000,	1000,	9000,	8000,	5000,	15000,	4000,	10000)
boxplot((results_allele_shift[,1]/num_mark)~sim_params[,7])
for (i in 1:3){
	print(mean(results_allele_shift[sim_params[,7]=="1000",i])/1000)
	print(mean(results_allele_shift[sim_params[,7]=="50",i])/15000)
	}
for (i in 4:9){
	print(mean(results_allele_shift[sim_params[,7]=="1000",i]))
	print(mean(results_allele_shift[sim_params[,7]=="50",i]))
	}

