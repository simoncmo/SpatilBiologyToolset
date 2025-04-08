# Make seurat object with coordinate from card
library(CARD)
library(Seurat)
# FUNCTION
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/CARD/ExtractCoordiate.R')

message("CreateSeuratFromCARD(card_obj) is loaded.")
# Function
CreateSeuratFromCARD <- function(card_obj){
	# 1. Get Coordinate and barcdoe
	coordinate_df = ExtractCoordinateFromCARDobject(card_obj)

	# 2. Use complicated colnames as rownames for coordinate
	coordinate_df_complex = coordinate_df
	rownames(coordinate_df_complex) = colnames(card_obj)

	# Creat object using count and coordiinate in metadata
	sn_new_obj = CreateSeuratObject(
		counts = card_obj@assays@data@listData$counts,
		meta.data = coordinate_df_complex
	)
	message("Cooridnate store in meta.data slot")
	return(sn_new_obj)
}
