# get average expression from object
library(optparse)

option_list = list(
	# sample name
	make_option(c("-s", "--sample"), type="character", default=NULL, help="Sample name"),
	# input path
	make_option(c("-i", "--input"), type="character", default=NULL, help="Input Seurat object"),
	# identity to take average expression
	make_option(c("-g", "--groupbys"), type="character", default='seurat_clusters', help="Identity to take average expression. common seperated for all idents"),
	# output file
	make_option(c("-o", "--output"), type="character", default=NULL, help="Output file"),
	# Assay use
	make_option(c("-a", "--assays"), type="character", default='RNA', help="Assays to use. seperate by comma. Default is SCT"),
	# optional input - metadata table path. if provided will be added to the exisiting obj@meta.data
	make_option(
		c("-t", "--metadata"),
		type = "character",
		default = NULL,
		help = "Path to metadata table. If provided, will be added to the existing obj@meta.data"
	)
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(opt)

library(tidyverse)
library(qs)
library(Seurat)
library(qs)

# # TEST
# opt = list(
# 	sample = 'S18-15142Fp1Us1_1',
# 	input = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/3_austin_obj/out/S18-15142Fp1Us1_1/S18-15142Fp1Us1_1_processed.qs',
# 	groupbys = 'seurat_clusters,Celltype',
# 	output = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/8_Epi_marker/3_avg_expression/out/'
# )
# obj_test = subset(obj, downsample = 10)

# Load object
obj = qread(opt$input, nthreads = 50)

# Add meta.data if file exists
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Metadata/add_missing_meta_v1.R')
if( file.exists(opt$metadata) ){
    message('Adding metadata from ', opt$metadata)
    meta = read_tsv(opt$metadata)
    obj = AddMissingMeta(obj, meta)
}

# Set assay
assays_all = strsplit(opt$assays, ',')[[1]]
assays_all = intersect(assays_all, Assays(obj)); print(assays_all)

# Process groupby
group_by_all = strsplit(opt$groupbys, ',')[[1]]
group_by_all = intersect(group_by_all, names(obj@meta.data)); print(group_by_all)


# Function
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Expression/AverageExpression/src/avg_exp_from_dot_plot.R')

# features defaults to all features
features_use = rownames(obj)

# Run Per assay 
for(assay_use in assays_all){
	message("Processing assay: ", assay_use)
	DefaultAssay(obj) = assay_use

	# Extract data
	exp_list = map(group_by_all, function(group_by){
		message("Extracting average expression for ", group_by, ", from N features: ", length(features_use))
		avg_exp = AvgExpressionFromDotPlot(obj, group.by = group_by, features = features_use)
		return(avg_exp)
	}) %>% setNames(group_by_all)

	# save to file
	iwalk(exp_list, function(df, ident){
		out_sample = str_glue('{opt$output}/{opt$sample}/{assay_use}/')
		dir.create(out_sample, showWarnings = FALSE, recursive = TRUE)

		message("Saving average expression for ", ident)
		write_tsv(df, file.path(out_sample, paste0('1_Avg_exp_',opt$sample, '__', ident, '.tsv')))
	})

}