message("MergeBPCellObjectList(obj_list, assay_keep_name = 'Xenium', ...)")
MergeBPCellObjectList = function(obj_list, assay_keep_name = 'Xenium', ...){
	# Remove any assay
	obj_list_clean = map(obj_list, function(obj){
		assay_keep = obj@assays[[assay_keep_name]]
		obj@assays = list()
		obj@assays[[assay_keep_name]] = assay_keep
		return(obj)
	})

	# TEST4: merge BPCells directly 
	obj_merge = merge(obj_list_clean[[1]], obj_list_clean[-1], add.cell.ids = names(obj_list_clean), ...)
	return(obj_merge)
}