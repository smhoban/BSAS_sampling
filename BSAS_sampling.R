
setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")

library(adegenet); library(diveRsity)
library(parallel); library(doParallel) #will load foreach
source("sample_funcs_BSAS.R"); source("src/arp2gen_edit.R")
colMax <- function(data) sapply(data, max, na.rm = TRUE)

#############################################################
#	PREP GENPOP FILES AND EMPTY RESULTS MATRICES			#
#############################################################

#Switch the comment to the other line if doing the stepping stone migration model 
#scenarios_ran<-list.dirs(path = "./mss_migr_model/", full.names = TRUE,recursive=F)
scenarios_ran<-list.dirs(path = "./Simulations6/", full.names = TRUE,recursive=F)

#This will convert all arp files to .gen (genepop) files
#If re-doing this you want to uncomment the first two lines which will wipe all .gen files
#If you want to wipe all .arp files after the conversion (to save disk space), uncomment the last line
for (scen in 1:length(scenarios_ran)){	
#	gen_files<-dir(path=scenarios_ran[scen],pattern="gen")
#	if (length(gen_files)!=0)  file.remove(file.path(scenarios_ran[scen],gen_files))
	reps_ran_arp<-list.files(scenarios_ran[scen], pattern="arp")
	arp_file_list<-file.path(scenarios_ran[scen],reps_ran_arp,sep="")
	if (.Platform$OS.type=="unix") arp_file_list<-substr(arp_file_list,1,nchar(arp_file_list)-1)
	mclapply(arp_file_list,arp2gen,mc.cores=16)
	#this deletes the arp files! DANGER!
	#file.remove(arp_file_list)
}
reps_ran_gen<-list.files(scenarios_ran[scen], pattern="gen")
	
#This is the range of sampling efforts to try, from 1 to 250
n_trees<-1:250
#Matrix to store results 
tree_results_temp<-matrix(nrow=length(n_trees),ncol=45)
all_res<-array(dim=c(length(scenarios_ran),45,length(reps_ran_gen)))
tree_results_all<-array(dim=c(length(n_trees),45,length(scenarios_ran),length(reps_ran_gen)))

##############################################
#	SAMPLING LOOP OVER ALL SIMULATION FILES	 #
##############################################

for (scen in 1:length(scenarios_ran)){	
	reps_ran_gen<-list.files(scenarios_ran[scen], pattern="gen")
	num_reps<-length(reps_ran_gen)
	
	cl <- makeCluster(28) # create a cluster with X cores
	registerDoParallel(cl) # register the cluster
	tree_results_list <- foreach (nrep = 1:num_reps, .packages="adegenet")	%dopar% {
	#If not using parallel computing, uncomment the below line for normal loop and comment out the above 3 lines
	#for (nrep in 1:num_reps){
		print(scenarios_ran[scen])

		#make a genind (by individual) and genpop (by population)
		temp_file_name<-file.path(scenarios_ran[scen],reps_ran_gen[nrep],sep="")
			if (.Platform$OS.type=="unix")  temp_file_name<-substr(temp_file_name,1,nchar(temp_file_name)-1)
		BSAS_genind<-read.genepop(temp_file_name,ncode=3)
		BSAS_genpop<-genind2genpop(BSAS_genind)


		#--NUMBER OF POPULATIONS, INDIVIDUALS, REGIONS, REGIONAL MAKEUP--#
		n_pops<-length(levels(BSAS_genind@pop))
		n_total_indivs<- length(BSAS_genind@tab[,1])
		n_ind_p_pop<-table(BSAS_genind@pop)
		allele_freqs<-colSums(BSAS_genpop@tab)/(n_total_indivs*2)	

		allele_cat<-get.allele.cat(BSAS_genpop,  n_ind_p_pop,local=F)	
		#global;			glob_com<-allele_cat;	glob_lowfr<-allele_cat;		glob_rare<-allele_cat
		#loc_com_int		loc_rare<-allele_cat	somew_com<-allele_cat;		somew_com<-allele_cat;			BM

		BSAS_genind_sep<-seppop(BSAS_genind)

		for (trees_samp in n_trees){
			all_caught<-matrix(nrow=n_pops,ncol=length(allele_freqs))
			for (p in 1:n_pops){
				do_replace<-table(BSAS_genind@pop)[p]<trees_samp
				if (trees_samp==1) 	all_caught[p,]<-BSAS_genind_sep[[p]]@tab[sample(1:nrow(BSAS_genind_sep[[p]]@tab), trees_samp),]
				else	all_caught[p,]<-colSums(BSAS_genind_sep[[p]]@tab[sample(1:nrow(BSAS_genind_sep[[p]]@tab), trees_samp,replace=do_replace),])
			}

			#as long as its not missed (0) it counts as being caught (also try 5, 10, 25, 50)
			for (l in 1:length(allele_cat)) tree_results_temp[trees_samp,l]<-sum(colSums(all_caught,na.rm=T)[allele_cat[[l]]]>0)
			for (l in 1:length(allele_cat))	
					tree_results_temp[trees_samp,l+length(allele_cat)]<-sum(colSums(all_caught,na.rm=T)[allele_cat[[l]]]>=5)
			for (l in 1:length(allele_cat))	
					tree_results_temp[trees_samp,l+length(allele_cat)*2]<-sum(colSums(all_caught,na.rm=T)[allele_cat[[l]]]>=10)
			for (l in 1:length(allele_cat))	
					tree_results_temp[trees_samp,l+length(allele_cat)*3]<-sum(colSums(all_caught,na.rm=T)[allele_cat[[l]]]>=25)
			for (l in 1:length(allele_cat))	
					tree_results_temp[trees_samp,l+length(allele_cat)*4]<-sum(colSums(all_caught,na.rm=T)[allele_cat[[l]]]>=50)
					
		}
			tree_results_temp<-t(t(tree_results_temp)/unlist(lapply(allele_cat,length)))
			tree_results_temp
	}
		gt95<- function(props) min(which(props>0.95))
			
		for (nrep in 1: num_reps) all_res[scen,,nrep]<-apply(tree_results_list[[nrep]],2,gt95)
		for (nrep in 1: num_reps) tree_results_all[,,scen,nrep]<-tree_results_list[[nrep]]
		stopCluster(cl)
}
#Again if running stepping stone model, switch comment out on next lines
#save(all_res,file="all_res_mss.Rdata") ;   save(tree_results_all,file="tree_results_mss.Rdata")
save(all_res,file="all_res_mult.Rdata") ;   save(tree_results_all,file="tree_results_mult.Rdata")
