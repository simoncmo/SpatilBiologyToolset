# load qs or RDS
library(qs)
library(stringr)

message("LoadObject( object_path, nthreads_to_use (default 50))")
LoadObject = function(path, nthreads = 50){
	message("Loading ", path)

	if(str_detect(path,'.qs')){
		message("Loading using qread with n threads = ", nthreads)
		
		object = qread(path, nthreads = nthreads)
	}else{
		message("Loading using readRDS")
		
		object = readRDS(path)
	}
}

