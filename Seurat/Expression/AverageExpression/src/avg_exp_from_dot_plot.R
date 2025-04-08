# Get Average expression with DotPlot function
AvgExpressionFromDotPlot = function(obj, group.by = 'seurat_clusters', features){
	pdot = DotPlot(obj, group.by = group.by, features = features)
	df = pdot$data %>% mutate(id = as.character(id))
	# Add n cell per cluster
	df = left_join(df,
		obj@meta.data %>% count(.data[[group.by]]) %>% setNames(c('id','n_cell')) %>% mutate(id = as.character(id)),
		by = c('id')
	)
	
	return(df)

}
message("AvgExpressionFromDotPlot(obj, group.by = 'seurat_clusters', features)")
