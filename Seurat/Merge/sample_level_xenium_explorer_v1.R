# Make sample level assignment and xenium explorer input
library(tidyverse)
library(optparse)

option_list = list(
	make_option(c("-m", "--meta"), type="character", default=NULL, 
							help="Input meta data"),
	# merged id colum name
	make_option(c("-i", "--id_column"), type="character", default=NULL, 
							help="Column name in meta data that contains barcode information"),
	# # column to keep as group in xenium explorer input
	# make_option(c("-c", "--cluster_column"), type="character", default=NULL, 
	# 						help="Column name in meta data that contains cluster information"),
	# output directory
	make_option(c("-o", "--output"), type="character", default=NULL,
							help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

print(opt)


# # TEST
# opt = list(
# 	meta = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/5_integration/2_Harmony_v2/Run_S18-15142-B17/out/4_assign_celltype/2_celltype_assignment.tsv',
# 	id_column = 'merged_id',
# 	output = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/5_integration/2_Harmony_v2/Run_S18-15142-B17/out/5_sample_level/'
# )

# Function to parse merged id
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Merge/split_merged_id_v1.R')

# Load 
meta = if(str_detect(opt$meta, 'tsv')) read_tsv(opt$meta) else read_csv(opt$meta)

# Process
# Note SplitMergedID2 will try to split by '_' first, if failed, it will try to split by '__'
meta_id = bind_cols(SplitMergedID2(meta[[opt$id_column]]), meta)
print(head(meta_id))

# Save
dir.create(opt$output, showWarnings = FALSE)
write_tsv(meta_id, file.path(opt$output, 'meta_sample_level_ID.tsv'))


# 2. Save sample level output and uses xenium_explorer format
meta_id_list = meta_id %>% split(f = .$Sample_ID)
map(meta_id_list, head)

idents_use = setdiff(names(meta), opt$id_column); print(idents_use)
for(ident in idents_use){
	iwalk(meta_id_list, function(df, name){
		# meta.data format
		df %>% select('Cell_ID', everything()) %>% 
		write_tsv(str_glue('{opt$output}/{name}_meta.tsv'))

		# Xenium explorer format
		df %>% 
			select('Cell_ID', all_of(ident)) %>% 
			setNames(c('cell_id','group')) %>% 
			write_csv(str_glue('{opt$output}/{name}_xenium_explorer_{ident}.csv'))
	})
}

