FetchDataFov = function(obj, values , fov_use = NULL){
	if(is.null(fov_use)) fov_use = Images(obj)[[1]]
	fov_obj = obj@images[[fov_use]]

	# v1. Use GetTissueCoordinates to get the cell id
	#fov_df = GetTissueCoordinates(fov_obj)
	#fov_cellid = fov_df$cell

	# Update with new method: directly get the cell id from the fov object
	# v2
	fov_cellid = fov_obj@boundaries$centroids@cells
	fov_df = as.data.frame(fov_obj@boundaries$centroids@coords) %>% mutate(cell = fov_cellid)

	value_df = FetchData(obj, values, cells= fov_cellid) %>% rownames_to_column('cell')
	left_join(fov_df, value_df, by = 'cell')
}
message("FetchDataFov(obj, values , fov_use = 'fov') ")
