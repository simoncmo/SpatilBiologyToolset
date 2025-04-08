# 2024-04-23
library(dplyr)
library(stringr)
library(readr)
library(purrr)


# Process name
message("SplitMergedID(merged_id)")
message("merged_id format: {Sample_ID}_{ID in [A-Z,a-z]-1}")
SplitMergedID = function(merged_id){
	# stop if null
	if(is.null(merged_id)){
		stop("merged_id is NULL")
	}
	split_list = merged_id %>% str_split('_(?=[a-z,A-Z]+-1)')

	if(length(split_list[[1]]) == 1){
		# split failed
		print(split_list %>% map_chr(1))
		stop("No pattern found. Please check input are {Sample_ID}_{ID in [A-Z,a-z]-1}")
	}
	data.frame(
		Sample_ID = split_list %>% map_chr(1),
		Cell_ID = split_list %>% map_chr(2)
	)
}

message("merged_id format: {Sample_ID}__{ID in [A-Z,a-z]-1}")
SplitMergedID2Dash <- function(merged_id){
	# stop if null
	if(is.null(merged_id)){
		stop("merged_id is NULL")
	}
	split_list = merged_id %>% str_split('__(?=[a-z,A-Z]+-[0-9]+)')

	if(length(split_list[[1]]) == 1){
		# split failed
		print(split_list %>% map_chr(1))
		stop("No pattern found. Please check input are {Sample_ID}__{ID in [A-Z,a-z]-1}")
	}
	data.frame(
		Sample_ID = split_list %>% map_chr(1),
		Cell_ID = split_list %>% map_chr(2)
	)
}

# Try first run first, if failed, try second run
# if success, return the result
SplitMergedID2 <- function(merged_id){
	tryCatch({
		SplitMergedID(merged_id)
	}, error = function(e){
		message("merged_id format: {Sample_ID}__{ID in [A-Z,a-z]-1}")
		SplitMergedID2Dash(merged_id)
	})
}



# Test 
# meta = read_tsv('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/5_integration/2_Harmony_v2/Run_S18-15142-B17/out/4_assign_celltype/2_celltype_assignment.tsv')

# SplitMergedID(meta$merged_id %>% head)

