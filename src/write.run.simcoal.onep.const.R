#NOTE THE LAST POPULATION TO GO BACK TO COALESCING IS HARD CODED
simc.write<-function(basefile,m_core,num_rep,num_mark,migMat,pop_sizes,time_coalesc, time_b,size_bb){
	filename<-basefile
	write("//Parameters for the coalescence simulation program : simcoal.exe",file=filename)
	write(paste(length(pop_sizes),"samples to simulate : mini uk"),file=filename,append=T)
	write("//Population effective sizes (number of genes 2*diploids)",file=filename,append=T)
	for (p in 1:length(pop_sizes)) write(pop_sizes[p],file=filename,append=T)
	write("//Samples sizes (number of genes 2*diploids)",file=filename,append=T)
	for (p in 1:length(pop_sizes)) write(pop_sizes[p],file=filename,append=T)
	write("//Growth rates	: negative growth implies population expansion",file=filename,append=T)
	for (p in 1:length(pop_sizes)) write(0,file=filename,append=T)
	write("//Number of migration matrices : 0 implies no migration between demes",file=filename,append=T)
	write("0",file=filename,append=T)
	write("//historical event: coalesc time, source, sink, migrants, new deme size, new growth rate, migration matrix index",file=filename,append=T)
	write("0 historical event",file=filename,append=T)
	write(paste("//Number of independent loci [chromosome]\n",num_mark," 1",sep=""),file=filename,append=T)
	for (l in 1:num_mark) write(paste("//Number of contiguous linkage blocks in chromosome ",(l-1),"\n1\n//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters\nMICROSAT 1  0.0005 5.0E-4 1.0 50",sep=""),file=filename,append=T)
}

simc.run<-function(basefile,num_rep){
	filename<-basefile
	if (.Platform$OS.type=="windows"){
		m_core=F
		if (m_core==T)
			err <- system(paste("fastsimcoal.exe -i ", filename, " -c4 -B4 -n ", num_rep, " -g -q",sep=""))
		if (m_core==F)
			err <- system(paste("fastsimcoal.exe -i ", filename, " -n ", num_rep, " -g -q",sep=""))
	}
	if (.Platform$OS.type=="unix"){		
		if (m_core==T)
			err <- system(paste("./fsc -i ", filename, " -c4 -B4 -n ", num_rep, " -g -q",sep=""))
		if (m_core==F)
			err <- system(paste("./fsc -i ", filename, " -n ", num_rep, " -g -q",sep=""))
	}
	if(err == 0) { cat("fastsimcoal exited normally\n")  } else { stop("fastsimcoal exited with error ", err, "\n") }
}


