# Add missing meta
# This function take meta table, assume first column to be barcode, 
# add any missing column that's not in the seurat object@meta.data to the meta
library(Seurat)

message("AddMissingMeta( obj, df )")
message("Assuming first column of the df is the barcode/cell-id")
AddMissingMeta = function(obj, df){
	# 1. make first column is the barcode/cell-id
	df = MakeIDinRowname(df)
	
	# 2. keep only columns not in the obj@meta.data
	columns_keep = setdiff(names(df), names(obj@meta.data))
	print(columns_keep)
	# if no columns to add, return obj
	if(length(columns_keep) == 0){
		message("No columns to add")
		return(obj)
	}

	# 3. AddMetaData to obj
	obj = AddMetaData(obj, df[,columns_keep, drop = F])

	return(obj)
}

# Make sure cell id in rownames. If not, default make the first column the rowname
MakeIDinRowname <- function(df, id_column_idx = 1){
	# Make sure is data.frame
	df = as.data.frame(df)

	# If have rownames. return as is
	if(HaveRownames(df)){
		return(df)
	}

	# Make the designated column the rowname
	rownames(df) = df[,id_column_idx]
	df = df[,-id_column_idx,drop = F]
	return(df)
}

# Check if rownames is generate 
# https://stackoverflow.com/questions/28777073/test-for-existing-row-names-and-col-names-in-data-frame
HaveRownames <- function(df){
	# a negative sign indicates the row names were generated automatically.
	return(.row_names_info(df) > 0)
}
