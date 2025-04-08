LoadXeniumhMulti = function(){
	df = read.table('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Xenium_panel/table/xenium_hMulti_v1.tsv',
		header = T, sep = '\t')
	return(df)
}

LoadXeniumhMultiGenes = function(){
	df = LoadXeniumhMulti()
	return(unique(df[,'Gene']))
}