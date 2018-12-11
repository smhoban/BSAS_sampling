######################################
#---FILE CHECK AND FILE CONVERSION---#
######################################

for (scen in 1:length(scenarios_ran)){	
	#Check for and remove genind
	gen_files<-dir(path=scenarios_ran[scen],pattern="gen")
	if (length(gen_files)!=0)  file.remove(file.path(scenarios_ran[scen],gen_files))
	#convert to genind	
	reps_ran_arp<-list.files(scenarios_ran[scen], pattern="arp")
	arp_file_list<-file.path(scenarios_ran[scen],reps_ran_arp,sep="")
	if (.Platform$OS.type=="unix") arp_file_list<-substr(arp_file_list,1,nchar(arp_file_list)-1)
	mclapply(arp_file_list,arp2gen,mc.cores=16)
	#foreach (nrep = 1:length(reps_ran_arp)) %dopar% {	arp2gen(arp_file_list[nrep])	}		#alternative MC loop
}
