

setwd("C:/Users/shoban/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")
setwd("/home/user/Dropbox/Projects/IN_PROGRESS/Best_sampling_across_systems/")

source("src/write.run.simcoal.bneck.R")
source("src/array.mat.funcs.R")
#source("fsc.code.R")

#I think this is where I want to start... 
#I want to read in a map and run a simcoal
 
base_filename<-"BSAS"
Np<-c(2,3,5,7,10,14,20)
m<-c(0.04,0.02,0.01,0.005,0.0025,0.00125,0.000625); migchar<-gsub("0.0","0",as.character(m))
Ps<-2*c(50,75,100,150,200,300,500,1000)
PbPs<-2*12500; PbPschar<-"L"	#Pre bottleneck Population size
gen_bneck<-c(25,5,1)				#generations ago when bottleneck occurred
mig_mod<-"mim"
parfilenames<-list(); parfoldernames<-list()
#Rather than have nested loops, lets make a reference table for all combinations
#and pull from each row of that table
scen_vars<-expand.grid(Np,m,Ps)

for (bn in gen_bneck){
for (i in 1:nrow(scen_vars))
	
parfilenames[i]<-paste("b",bn,PbPschar,base_filename,"_Np_",scen_vars[i,1],"_",mig_mod,"_",gsub("0.0","0",as.character(scen_vars[i,2])),"_Ps_",scen_vars[i,3]/2,".par",sep="")
#for (i in 1:nrow(scen_vars)) parfoldernames[i]<-paste(base_filename,"_Np_",scen_vars[i,1],"_",mig_mod,"_",gsub("0.0","0",as.character(scen_vars[i,2])),"_Ps_",scen_vars[i,3]/2,"_b",gen_bneck,PbPschar,sep="")


num_mark<-rep(c(15000,10000,10000,9000,8000,5000,4000,1500,500,450),each=49)
num_rep<-3
m_core<-T
for (i in 1:nrow(scen_vars)){
	this_migMat<-matrix(rep(scen_vars[i,2],(scen_vars[i,1])^2),nrow=scen_vars[i,1])
	this_migMat[ row(this_migMat) == col(this_migMat) ] <- 0
	time_coalesc<-50000
	num<-1
	filename<-parfilenames[[i]]
	simc.write(filename,num_rep,round((num_mark[i]/scen_vars[i,1])*1.5),this_migMat,rep(scen_vars[i,3],scen_vars[i,1]),time_coalesc,bn,PbPs)
	simc.run(filename,num_rep,m_core)
}
}

