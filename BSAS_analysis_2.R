

 
#####################################
#			ANALYSIS 
#####################################
		
setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
load("all_res_mult.Rdata")
#Remember the composition of the results matrix- this is the minumum sampling needed to capture 95% of the alleles 
#with columns as the minimum needed in a given category
#For example, once (cols 1:9), five times (cols 10:18), or 10 or 25 or 50 times 
#	There are nine allele categories able to calculate (see below) even though only four are used in the paper
#Rows are the scenarios (combinations of parameters)
			
#Some allele type/ parameter combinations have Inf/ NAs for many/ most replicates, which means that category not caught
#This replaces all Inf values (never reached 0.95) with NAs
is.na(all_res) <- sapply(all_res, is.infinite)
means_caught<-apply(all_res[,,],c(1,2),mean,na.rm=T)
#This will count the NAs- important to report!!!
apply(is.na(all_res[,,]),c(1,2),sum)
#run this and pull the params and then do table on that to see how many large/ small populations, bottleneck, migr etc.
#i.e. which situations is it most likely
apply(is.na(all_res[,,]),c(1,2),sum)[,3]>=15  



#####################################
#	MAIN PLOTS BY PARAMETER			#
#####################################  
		
scenarios_ran<-list.dirs(path = "./Simulations6/", full.names = TRUE,recursive=F)
#set up the params list
get_params<- function(scenario) strsplit(strsplit(scenario,"/")[[1]][4],"_")
lapply(scenarios_ran,get_params)
sim_params<-matrix(unlist(lapply(scenarios_ran,get_params)),ncol=7,byrow=T)
sim_params[,5]<-paste(".",sim_params[,5],sep="")
allele_categories<-c("all alleles (including rare)", "overall common f>0.05", "overall low fr 0.10>f>0.01", 
					 "overall rare f<0.01", "local1", "local2", "somwhere common f>0.05", "somwhere common f>0.05", "BM overall f~0.05")

summary(lm(sim_params[,3],means_caught[,2]))
plot(sim_params[,5],means_caught[,2])
plot(sim_params[,7],means_caught[,2])
	
#Or
#sim_params<-read.csv("sim_params.csv")[,-1]	
	
#####################################################
#	FIGURE 1 AND SUPP PLOTS- GRID ALL PARAMS		#
#####################################################
	
#PDFs showing each factor for the four kinds of alleles- will reproduce Figure 1

#Can do separately for the three kinds of bottlenecks by looking at rows for means_caught
	#in the inner loop just paste in 1:392 for rows of means_cause and sim_params etc.
	#into the below boxplot code to subset by bottleneck, otherwise the code calculates across all bottlenecks
	#and of course change the file output name to B1, B5, or B25 for bottleneck times
	#These are the different bottleneck scenarios simulated... 
	#B1 [1:392,]
	#B25 [393:784,]
	#B5 [785:1176,]

#Can also do separately for the multiple copies (for example for five copies) 
	#Just change the AT loop i.e. 1,3,8,9 sub for 10,12,17,18
	#and of course change the file name

pdf(file="BSAS_grid_final_1cop.pdf", height=12, width=14)
par(mfcol=c(4,3), mar=c(2,2,2,2), oma=c(6,5,4,4), cex.axis=2)
for (P in c(3,5,7)){
	do_yaxis="s"; do_xaxis="n"
	for (AT in c(1,3,8,9)){
		if (P>3) do_yaxis="n"; if(AT==9|AT==18) do_xaxis="s"
		boxplot(means_caught[,AT]~as.numeric(sim_params[,P]),xaxt=do_xaxis,yaxt=do_yaxis); if (P==7) axis(4,labels=F)
}	}
mtext("number of individuals needed", side=2, line=2,outer=T,cex=1.75)
mtext("Cat 4- species >0.05           Cat 3- local >0.05            Cat 2- 0.01 to 0.10          Cat 1- all alleles       ", side=4, line=2,outer=T,cex=1.4)
mtext("Number of populations                                 Migration rate                                     Population size", side=1, line=2,outer=T,cex=1.4)

dev.off()
	
	

	
	
##############################################################
#	POP SIZE AND MIGRATION COMBINED- FIGURE 2				 #
##############################################################

#The relatively flat boxplots for migration rate and population size in Figure 1 suggest that these factors 
#are dominated by the influence of number of populations which makes sense because the samples are 
#divided up among the different populations 
#So we need to parcel out these factors, to see if they have some influence
#The following will go through every combination of migration rate and population size
	
sim_params[sim_params[,7]==50,7]<-"  50";			sim_params[sim_params[,7]==75,7]<-"  75";		sim_params[sim_params[,7]==150,7]<-" 150"
sim_params[sim_params[,7]==100,7]<-" 100";		sim_params[sim_params[,7]==200,7]<-" 200";	sim_params[sim_params[,7]==300,7]<-" 300"
sim_params[sim_params[,7]==400,7]<-" 400";		sim_params[sim_params[,7]==500,7]<-" 500";		

#This one for publication- focuses on local alleles, has grey for separating each population size
psize_mig<-paste("p=",sim_params[,7],"  m=",sim_params[,5],sep="")
pdf("vary_migr_by_psize.pdf",height=6,width=12)
par(mar=c(10,5,3,2))
AT<-8; boxplot(means_caught[,AT]~psize_mig,las=2,col=c(rep("light grey",7),rep("white",7)),ylab="Ni, number of individuals to sample per population")
dev.off()
	
	#This one NOT for publication- has separate page for each allele category- just for examining
sim_params_temp<-cbind(sim_params,paste(sim_params[,5],sim_params[,7],sep="-"))
pdf("vary_psize_by_migr.pdf",height=7,width=15)
for (AT in c(1,2,3,4,8,9)) boxplot(means_caught[,AT]~sim_params_temp[,8],las=2,main=allele_categories[AT])
dev.off()

#This one NOT for publication-  focus on number of populations (i.e. there are separate graphs for each number of populations)
#To show large variation is due of course to number of populations 
#but not a very helpful set of plots
#Change the 7 to a 5 to vary by migration rate
#Change the AT for Allele Type (could look at any alleles or catch multiple)

pdf(file="BSAS_grid_final_psize.pdf", height=12, width=18)
par(mfcol=c(4,7), mar=c(2,2,2,2), oma=c(4,4,4,4))
num_pops_char<-unique(sim_params[,3])

#go through pop numbers, then allele types
for (nump in num_pops_char)		for (AT in c(1:4)){
		boxplot(means_caught[sim_params[,3]==nump,AT]~as.numeric(sim_params[sim_params[,3]==nump,5]),las=2,main=nump)
}	
dev.off()
	
	
#################################################
#	MULT FACTOR	CATCH 1,5,10,25,50- FIGURE 3	#
#################################################

#What is the multiplication factor for catching different numbers of alleles?
#In other words what is the ratio of number of samples needed to catch 1 vs. 5 (or vs. 10 or 25 or 50)
#across the types of alleles and all scenarios
pdf(file="mult_factor_boxplots.pdf",height=5,width=15)
par(oma=c(5,3,3,3),mfrow=c(1,4),cex.axis=1.5,cex.lab=1.7,cex.main=1.7)
boxplot(means_caught[,c(10,12,17,18)]/means_caught[,c(1,3,7,9)], names=1:4,main="5 copies",col="light grey"); abline(h=5,lty=2,col="salmon")
boxplot(means_caught[,9+c(10,12,17,18)]/means_caught[,c(1,3,7,9)],main="10 copies", names=1:4,col="light grey"); abline(h=10,lty=2,col="salmon")
boxplot(means_caught[,18+c(10,12,17,18)]/means_caught[,c(1,3,7,9)],main="25 copies", names=1:4,col="light grey"); abline(h=25,lty=2,col="salmon")
boxplot(means_caught[,27+c(10,12,17,18)]/means_caught[,c(1,3,7,9)],main="50 copies", names=1:4,col="light grey"); abline(h=50,lty=2,col="salmon")
mtext("Category of allele", outer=T,side=1,line=-0.5,cex=1.2)
	mtext("1-all alleles, 2-low frequency (0.01-0.10), 3-locally >0.05, 4-species wide >0.05", outer=T,side=1,line=2,cex=1.2); 
mtext("size by which collection must be increased",outer=T,side=2,cex=1.2)
 mtext("number of allele copies desired in the collection",outer=T,side=3,cex=1.2)
dev.off()

#Calculations for reporting
mean(means_caught[,c(10,12,17,18)]/means_caught[,c(1,3,7,9)],na.rm=T)	#4.16
mean(means_caught[,9+c(10,12,17,18)]/means_caught[,c(1,3,7,9)],na.rm=T)		#7.19
mean(means_caught[,18+c(10,12,17,18)]/means_caught[,c(1,3,7,9)],na.rm=T)	#15.8
mean(means_caught[,27+c(10,12,17,18)]/means_caught[,c(1,3,7,9)],na.rm=T)	#27.4



#################################################################################################

#############################################################
# EVERYTHING ELSE IS "CHECKS" TO MAKE SURE IT ALL WORKED	#
#############################################################

	
#####################################
#	MAKE SURE BOTTLENECK WORKED		#
#####################################
#At first it seems surprising that the bottleneck length had so little effect on sampling 
#In spite of the known effect of bottlenecks on number of alleles 	
#To make sure they worked, look at bottlenecks and migration examine genepop files 
#Calculate NUMB ALLELES and OBS HET over 10 reps, compare among bnecks, for 12 random scenarios

 setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
 load("all_res_final.Rdata")
 library(adegenet); library(diveRsity)
 scenarios_ran<-list.dirs(path = "./Simulations6/", full.names = TRUE,recursive=F)

#index 1,393,785 is the first scenario for the three bottlenecks
bneck_nall<-array(dim=c(3,12,10)); bneck_het<-array(dim=c(3,12,10))
bneck_base<-c(1,393,785); bneck_add<-seq(1,300,25)
#The 12 scenarios are:
# scenarios_ran[bneck_add]
# [1] "./Simulations6//b1LBSAS_Np_10_mim_000625_Ps_100" 	[2] "./Simulations6//b1LBSAS_Np_10_mim_005_Ps_1000"
# [3] "./Simulations6//b1LBSAS_Np_10_mim_04_Ps_150"		 	[4] "./Simulations6//b1LBSAS_Np_14_mim_0025_Ps_200"
# [5] "./Simulations6//b1LBSAS_Np_14_mim_02_Ps_300"			[6] "./Simulations6//b1LBSAS_Np_20_mim_00125_Ps_50"
# [7] "./Simulations6//b1LBSAS_Np_20_mim_01_Ps_500"			[8] "./Simulations6//b1LBSAS_Np_2_mim_000625_Ps_75"
# [9] "./Simulations6//b1LBSAS_Np_2_mim_01_Ps_100"			[10] "./Simulations6//b1LBSAS_Np_3_mim_000625_Ps_1000"
# [11] "./Simulations6//b1LBSAS_Np_3_mim_005_Ps_150"		[12] "./Simulations6//b1LBSAS_Np_3_mim_04_Ps_200"


for (b in 1:3){
	for (add in 1:12){
		this_b<-bneck_base[b]+bneck_add[add]
		setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
		reps_ran_gen<-list.files(scenarios_ran[this_b], pattern="gen")
		setwd(scenarios_ran[this_b])
		for (r in 1:10){
			BSAS_genind<-read.genepop(reps_ran_gen[r],ncode=3)
			#get alleles and heterozygosity, save
			bneck_nall[b,add,r]<-mean(unlist(summary(BSAS_genind)[4]))
			bneck_het[b,add,r]<-mean(unlist(summary(BSAS_genind)[6]))
}	}	}
#RESULTS: These plots show the ratio of non bottleneck to bottleneck- we can see it ranges 
#Thus the bottlenecks are having an effect on genetic diversity 
#Explanation is that bottlenecks only cause the loss of (mostly) rare alleles and those are likely
#below the threshold we are counting
setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
het_mean<-apply(bneck_het,c(1,2),mean,na.rm=T); nall_mean<-apply(bneck_nall,c(1,2),mean,na.rm=T)
pdf(file="bneck_nall.pdf");  boxplot(t(bneck_nall[1,,]/bneck_nall[2,,])); dev.off()
pdf(file="bneck_het.pdf");  boxplot(t(bneck_het[1,,]/bneck_het[2,,])); dev.off()
save(bneck_het,file="bneck_het_comp.Rdat"); save(bneck_nall,file="bneck_nall_comp.Rdat")
#het_mean= 1.03, max= 1.13	;	nall_mean= 1.26, max= 1.86


	#########################
	#COMPARE MSS AND MIM
	###########################
#Most simulations regard the island model.  This will analyze simulations under the stepping stone model (mss) and compare sampling needed
	
setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/Additional_checks/")
get_params<- function(scenario) strsplit(strsplit(scenario,"/")[[1]][4],"_") 

#First get the MSS data, then the MIM data
scenarios_ran<-list.dirs(path = "./mss_migr_model/", full.names = TRUE,recursive=F) 
load("all_res_mss.Rdata")
#Remember what is in these results- this is the minumum sampling needed to capture 95% of the alleles in a given category
#Either once (cols 1:9) or five times (cols 10:18) ... etc.
is.na(all_res) <- sapply(all_res, is.infinite)
means_caught_mss<-round(apply(all_res[,,],c(1,2),mean,na.rm=T),2)
sd_caught_mss<-apply(all_res[1:16,8,],1,sd,na.rm=T)	
sim_params_mss<-matrix(unlist(lapply(scenarios_ran,get_params)),ncol=7,byrow=T)
	
setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
scenarios_ran<-list.dirs(path = "./Simulations6/", full.names = TRUE,recursive=F)
load("all_res_final.Rdata")
is.na(all_res) <- sapply(all_res, is.infinite)
means_caught_mim<-round(apply(all_res[,,],c(1,2),mean,na.rm=T),2)
sim_params_mim<-matrix(unlist(lapply(scenarios_ran,get_params)),ncol=7,byrow=T)

#Recall that the MSS model was only for a subset of parameters, so narrow the data to parameters in common with MIM
shared_bw<-which(do.call("paste",as.data.frame(sim_params_mim[,c(1,3,5,7)])) %in% do.call("paste",as.data.frame(sim_params_mss[,c(1,3,5,7)])))
sd_caught_mim<-apply(all_res[shared_bw[1:16],8,],1,sd,na.rm=T)
	
#Analysis
diff_mss<-means_caught_mim[shared_bw,c(1,3,8,9)]/means_caught_mss[,c(1,3,8,9)]

colMeans(diff_mss[1:16,]); 
#[1] 1.0026127 1.0257867 0.7466511 1.0027878
colMeans(diff_mss[17:32,])
#[1] 1.0017947 0.9992222 1.0003993 0.9986894	#as expected there is no difference for the two population system
colMeans(diff_mss[33:48,])
#[1] 1.0038296 1.0010115 0.9471347 1.0083351
boxplot(diff_mss[,3]~as.numeric(paste(".",sim_params_mss[,5],sep="")))	#difference is highest for high gene flow- makes sense- stepping stone model is very restrictive in this case. 0.4 is much different than 0 while 0.000625 is not much different than 0
boxplot(diff_mss[,3]~as.numeric(sim_params_mss[,7]))					#difference is highest for larger population sizes- gene flow has most effect at larger population sizes, whereas at small sizes all populations experience high drift even if gene flow is high
#Which scenarios are especially different?
	(means_caught_mim[shared_bw,c(1,3,8,9)]/means_caught_mss[,c(1,3,8,9)])<.95
	(means_caught_mim[shared_bw,c(1,3,8,9)]/means_caught_mss[,c(1,3,8,9)])>1.05
#Create a summary table
	a<-matrix(cbind(
		round(means_caught_mss[,17]/means_caught_mim[shared_bw,17],2)[1:16],
		round(means_caught_mss[,17]/means_caught_mim[shared_bw,17],2)[33:48],
		sim_params_mss[1:16,c(5,7)]),ncol=4)
	colnames(a)<-c("14 pop","5 pop", "migr", "pop size") 
	mean(as.numeric(a[,1])); mean(as.numeric(a[,2]))
	
	
#############################################################
#	TO COUNT NUMBER OF SINGLETONS FOR DISCUSSION			#
	#########################################################

#To count the number of alleles occurring only once or twice in a sample and thus likely to be lost

library(adegenet); library(diveRsity)
source("sample_funcs_BSAS.R"); source("src/arp2gen_edit.R")
colMax <- function(data) sapply(data, max, na.rm = TRUE)

setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/Additional_checks/onep_no_bn/nob5LBSAS_Np_1_Ps_1000/")

BSAS_genind<-read.genepop("nob5LBSAS_Np_1_Ps_500_1_3.gen",ncode=3)

sum(colSums(BSAS_genind@tab[sample(1:nrow(BSAS_genind@tab), 500),])==1)
sum(colSums(BSAS_genind@tab[sample(1:nrow(BSAS_genind@tab), 500),])==0)
length(colSums(BSAS_genind@tab[sample(1:nrow(BSAS_genind@tab), 500),]))
hist(colSums(BSAS_genind@tab[sample(1:nrow(BSAS_genind@tab), 500),]),breaks=c(0,seq(5,100,by=5)))

	
#####################################
#	FST CHECKING MIGR RATE			#
#####################################
#Ok now do FSTs based on migration rate
#It seems odd that the migration rate has little affect on how to sample for allele capture
#So we want to make sure the migration rate is actually affecting the FST 
#We will look at several population sizes, for several numbers of populations, across the range of mig rates
#RESULTS: migration does influence FST, but less than population size, i.e. allele frequencies don't change
#a lot across the range of migration rates, as the populations are already quite small
#In small populations the population size will have more influence on drift than migr!
#(Also, again, some rare alleles might be lost or shared due to migration but we're ignoring rare alleles)
#To report migration rate in the results also

 library(hierfstat)
 pops_focus<-c("50","100","300","500")
 migr_comp<-array(dim=c(4,7,10))
 
for (pnum in c("2","3","5","10")){
	for (psize in 1:4){
		for (migr in 1:7){
				setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems")
				setwd(scenarios_ran[sim_params[,3]==pnum&sim_params[,1]=="b5LBSAS"&sim_params[,7]==pops_focus[psize]][migr])
				reps_ran_gen<-list.files(pattern="gen")
				for (r in 1:10){
				BSAS_genind<-read.genepop(reps_ran_gen[r],ncode=3)
				migr_comp[psize,migr,r]<-mean(pairwise.fst(BSAS_genind))
	}	}	}
	fst_mean<-apply(migr_comp,c(1,2),mean,na.rm=T)
	colnames(fst_mean)<- unique(sim_params[,5]) ; rownames(fst_mean)<-pops_focus
	setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
	write.csv(fst_mean,file=paste("fst_w_mig",pnum,".csv"))
}



#################################
#	ANALYZING ONE POPS			#
#################################
setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/Additional_checks")

#set up the params list
pdf("one_population.pdf", height=6, width=9);	par(mfrow=c(1,2))
load("all_res_onep_bn.Rdata")	
is.na(all_res) <- sapply(all_res, is.infinite)
means_caught<-apply(all_res[,,],c(1,2),mean,na.rm=T)
boxplot(means_caught[,c(1,3,9)],ylim=c(0,250),ylab="number plants to sample", names=c("all","low fr","0.05"),main="bottleneck, 1 population")
axis(side=4,at=c(0,50,100,150,200,250),labels=F); 	mtext("allele category",side=1,line=3,cex=1.3)
load("all_res_onep_nobn.Rdata")
is.na(all_res) <- sapply(all_res, is.infinite)
means_caught<-apply(all_res[,,],c(1,2),mean,na.rm=T)
boxplot(means_caught[,c(1,3,9)],ylim=c(0,250),ylab="number plants to sample", names=c("all","low fr","0.05"),main="no bottleneck, 1 population")
axis(side=4,at=c(0,50,100,150,200,250),labels=F); 	mtext("allele category",side=1,line=3,cex=1.3)
dev.off()

#For reporting, one population, for BM, compare to theoretical expecations
summary(c(all_res[,9,])) #mean of 27.39 median of 28;	mean of 27.9 median of 29	
sd(c(all_res[,9,]))	#sd 2.67	2.47
#good, it does conform to expectations of approximately 28 individuals

#and the ratio- the "big result" of the paper in terms of capturing multiple alleles- in this case five copies
mean(means_caught[,10:18]/means_caught[,1:9],na.rm=T)
#no bn then bn		#[1] 4.824	[1] 4.136
boxplot(means_caught[,10:18]/means_caught[,1:9])



#################################################
# ALLELE FREQUENCY SHIFT and NUMBER ALLELES		#
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
#	With bottleneck		#
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



	