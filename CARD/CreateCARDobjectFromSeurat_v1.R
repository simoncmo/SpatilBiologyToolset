# 2024-05-21 Simon Mo
# Create CARD object seurat wrapper
# tutorial: https://yma-lab.github.io/CARD/documentation/04_CARD_Example.html
library(tidyverse)
library(Seurat)
library(CARD)

# function
RemoveSpotNoCoordinates = function(ST, image = NULL){
	# Use first one if NULL
	image = image %||% names(ST@images)[[1]] 

	# If expression matrics ID the same number as location ID. return as is
	exp_cells = Cells(ST)
	loc_cells = Cells(ST@images[[image]])
	if(length(exp_cells) == length(loc_cells) ) return(ST) 

	# else, subset to only contains the intersect
	message("Expression Cells count: ", length(exp_cells))
	message("Coordinate Cells count: ", length(loc_cells))
	shared_cells = intersect(exp_cells, loc_cells)
	message("Shared Cell count:", length(shared_cells))
	return(subset(ST, cells = shared_cells))
	
}

ExtractSpatialCoordinate = function(obj, ST_slice = NULL){
	if(is.null(ST_slice)) ST_slice = Images(obj)[[1]]
	# check ST type
	# a. Visium
	tissue_df = GetTissueCoordinates(obj)
	if('imagerow' %in% names(tissue_df)){
		message('Visium')
		spatial_location = obj@images[[ST_slice]]@coordinates
		spatial_location = spatial_location[,c('imagecol','imagerow')] %>% setNames(c('x','y'))
	}else if('cell' %in% names(tissue_df)){
	# b. Xenium
		message('Xenium')
		spatial_location = tissue_df %>% column_to_rownames('cell') %>% .[,c('x','y')]
	}else{
		stop('ST SeuratObj type not supported!')
	}
	return(spatial_location)
}

createCARDObjectFromSeuratObject = function(
	sc_object,
	ST_object,
	sc_meta = NULL,
	sn_assay = 'RNA',
	ST_assay = 'Spatial',
	ST_slice = NULL,
	ct.varname = "", # Cell type column 
	sample.varname = NULL, # sample info column
	minCountGene = 100,
	minCountSpot = 5,
	...
	){

	# 1. Extract sndata
	sc_count = GetAssayData(sc_object@assays[[sn_assay]])

	if(is.null(sc_meta)) sc_meta = sc_object@meta.data

	# 1a. Check if ct.varname in metadata
	if(!ct.varname %in% colnames(sc_meta)){
		message('ct.varname ', ct.varname, ' not found!')
		message("available columns: ")
		print(colnames(sc_meta))
		stop('Please double check ct.varname value')
	}

	# 1b. check sample info column
	if(is.null(sample.varname)){
		message('sample.varname not provided. Assume single sample in the reference')
		message("Setting sample.varname as 'sampleInfo' and assign values as 'sample1'")
		sample.varname = "sampleInfo"
		sc_meta$sampleInfo = "sample1"
	}

	# 2a. Make sure ST data and ST coordinate has same row
	ST_object = RemoveSpotNoCoordinates(ST_object)

	# 2. Extract ST data
	if(!ST_assay %in% names(ST_object@assays)){
		message('ST assay ', ST_assay, ' not found!')
		message("available assays: ", colnames(ST_object@assays))
		stop('Please double check ST_assay value')
	}
	spatial_count = GetAssayData(ST_object@assays[[ST_assay]])

	# 3. Extract ST coordinate 
	spatial_location = ExtractSpatialCoordinate(ST_object, ST_slice)

	# DEBUG:
	# print(head(sc_count[1:3,1:3]))
	# print(head(sc_meta[1:3,1:3]))
	# print(head(spatial_count[1:3,1:3]))
	# print(head(spatial_location[1:3,]))
	
	# Create object 
	CARD_obj = createCARDObject(
		sc_count = sc_count,
		sc_meta = sc_meta,
		spatial_count = spatial_count,
		spatial_location = spatial_location,
		ct.varname = ct.varname,
		ct.select = unique(sc_meta[[ct.varname]]),
		sample.varname = sample.varname,
		minCountGene = minCountGene,
		minCountSpot = minCountSpot,
		...
		) 
	# QC on scRNASeq dataset! ...
	# QC on spatially-resolved dataset! ..

}

createCARDandRunDeconvFromSeurat = function(
	sc_object,
	ST_object,
	sc_meta = NULL,
	sn_assay = 'RNA',
	ST_assay = 'Spatial',
	ST_slice = NULL,
	ct.varname = "", # Cell type column 
	sample.varname = NULL, # sample info column
	minCountGene = 100,
	minCountSpot = 5,
	...
){
	# Create CARD object
	CARD_obj = createCARDObjectFromSeuratObject(
		sc_object = sc_object,
		ST_object = ST_object,
		sc_meta = sc_meta,
		sn_assay = sn_assay,
		ST_assay = ST_assay,
		ST_slice = ST_slice,
		ct.varname = ct.varname,
		sample.varname = sample.varname,
		minCountGene = minCountGene,
		minCountSpot = minCountSpot,
		...
	)

	# Run Deconv
	CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
	return(CARD_obj)
}

