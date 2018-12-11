
#######################
#---#DETERMINE WHAT ALLELES FALL IN WHAT CATEGORIES---#
#######################

get.allele.cat<-function(UK_genpop, n_ind_p_pop,local=TRUE){

	#--Set up categories for GLOBAL ALLELES
	allele_freqs<-colSums(UK_genpop@tab)/(n_total_indivs*2)	
	glob<-1:length(allele_freqs)
	glob_com<-as.vector(which(allele_freqs>0.05,arr.ind=T))
	glob_lowfr<-as.vector(which(allele_freqs<0.10&allele_freqs>0.01))
	glob_rare<-as.vector(which(allele_freqs<0.01,arr.ind=T))

	#--LOCAL ALLELES--#
	#get the populations having alleles which occur at frequency greater than a given frequency in focal pop'n
	#and not more than another frequency in more than two other populations
	if (local==TRUE) {
		local_freqs<-UK_genpop@tab/as.vector(n_ind_p_pop*2)
		#identify columns (alleles) for which at least one population exceeds min threshold
		loc_com_int<-as.vector(which((colMax(as.data.frame((local_freqs)))>=0.20)&(apply(local_freqs,2,sort)[(n_pops*.95),]<0.10),arr.ind=F))
		loc_rare<-as.vector(which(colMax(as.data.frame((local_freqs)))<0.10))
	}	else {loc_com_int<-NA; loc_rare<-NA}
	
	somew_comm<-as.vector(which(colMax(as.data.frame((UK_genpop@tab)/as.vector(n_ind_p_pop*2)))>0.05))
	somew_comm2<-as.vector(which(apply(as.data.frame((UK_genpop@tab)/as.vector(n_ind_p_pop*2)),2,function(col) sum(col>0.05)==1)))

	BM<-as.vector(which((allele_freqs>0.04&allele_freqs<0.06)==TRUE))

	list(glob,glob_com, glob_lowfr, glob_rare, loc_com_int,loc_rare,somew_comm, somew_comm2,BM)
}		


#####################
#--GENERAL SAMPLING FUNCTION--#
#####################

#we pass it a vector of populations and a vector of sample sizes, as well as the genind object

sample.pop<-function(genind_obj,vect_pop_ID,vect_samp_sizes){
	all_caught<-matrix(nrow=n_pops,ncol=length(allele_freqs))
	for (p in 1:length(vect_pop_ID)){
		if (vect_samp_sizes[p]==1) 	all_caught[p,]<-genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),]
		else	all_caught[p,]<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),])
	all_caught
	}
}
#all_caught[p,]<-BSAS_genind_sep[[1]]@tab[sample(1:50, 1),]
