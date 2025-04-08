# Description: Assign cell types based on RCTD scores
message("Assigning cell types based on RCTD scores")
message('assign_RCTD_celltypes(obj, combine_mode = "double", single_comp_cutoff = 0.5, double_comp_cutoff = 0.15, ignore_order = TRUE)')
### FUNCTION VERSION
# Function to assign cell type based on cutoffs
assign_cell_type <- function(row, 
	single_comp_cutoff=0.5, double_comp_cutoff=0.2,
	combine_mode = c('double','single','all'),
	ignore_order = TRUE
	) {
	combine_mode = match.arg(combine_mode)
	# If single model just return the top one
	if(combine_mode == 'single') return(names(row)[which.max(row)])
	
  if (sum(row > single_comp_cutoff) >= 1) {
    return(names(row)[which.max(row)])
  } else if (sum(row > double_comp_cutoff) >= 1) {
	# COMBINE:
	if(combine_mode == 'all'){ # Get all cell type passed the cutoff
		celltype_vector = names(row)[row > double_comp_cutoff]
	}else{
		# combine only top 2
		celltype_vector = names(row)[order(row, decreasing = TRUE)[1:2]]
	}
	if(ignore_order) celltype_vector = sort(celltype_vector)
	return(paste(celltype_vector, collapse = " & "))
  } else {
	return("Mixture")
  }
}

#library(Seurat)

assign_RCTD_celltypes <- function(obj, 
	combine_mode = c("double","single","all"), 
	single_comp_cutoff = 0.5, double_comp_cutoff = 0.15,
	ignore_order = TRUE,
	confident_cutoff = 0.55 # for sinlge mode only
	) {
  # Get RCTD celltype
  rctd_df <- obj@meta.data %>% select(starts_with('rctd_', ignore.case = FALSE))
  # Remove rctd_
  colnames(rctd_df) <- gsub('rctd_', '', colnames(rctd_df))

  # combine_mode
  combine_mode = match.arg(combine_mode)
  if(combine_mode == 'single') message("Warning, currently doesn't filter but return all the cell with max values")
  
  # Apply the function row-wise
  assigned_df <- rctd_df
  assigned_df[[str_c('RCTD_celltype_',combine_mode)]] <- apply(rctd_df, 1, 
  	assign_cell_type, 
	combine_mode = combine_mode,
	ignore_order = ignore_order,
  	single_comp_cutoff = single_comp_cutoff, 
	double_comp_cutoff = double_comp_cutoff)

  # If mode is Single, also assign the proportion
  # And a column to keep only the 'confident spot'
  if(combine_mode == 'single'){
	assigned_df$RCTD_celltype_proportion = apply(rctd_df, 1, max, na.rm=T)
	print(head(assigned_df$RCTD_celltype_proportion))
	# Assign confident celltype
	assigned_df = assigned_df %>% mutate(
		RCTD_celltype_confident = ifelse(RCTD_celltype_proportion >= confident_cutoff, RCTD_celltype_single, 'Mixed')
	)
  }
  
  # Add cell_id column
  assigned_df$cell_id <- rownames(assigned_df)
  
  message("save result to column: ", str_c('RCTD_celltype_',combine_mode))
  # Return the assigned dataframe
  return(assigned_df %>% select(cell_id, contains('RCTD_celltype')) )
  # return(assigned_df[, c("cell_id", paste0("RCTD_celltype_", combine_mode))])
}


# Example usage:
# obj <- LoadYourSeuratObject()
# assigned_celltypes <- assign_RCTD_celltypes(obj)
# assigned_celltypes <- assign_RCTD_celltypes(obj, combine_mode = "double")
# assigned_celltypes %>% count(RCTD_celltype_double, sort = T)
