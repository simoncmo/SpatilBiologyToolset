# Assign cell type to barcode
library(readr)
library(dplyr)
library(Seurat)
library(optparse)

# Take meta data, celltype assignment file, id_column, cluster_column, and output path
options <- list(
    make_option(c("-m", "--meta"), type="character", default=NULL,
        help="Input meta data"),
    make_option(c("-i", "--id_column"), type="character", default=NULL,
        help="Column name in meta data that contains barcode information"),
    make_option(c("-c", "--cluster_column"), type="character", default=NULL,
        help="Column name in meta data that contains cluster information"),
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="Output directory"),
    make_option(c("-a", "--assignment"), type="character", default=NULL,
        help="Cell type assignment file")
)   

opt_parser <- OptionParser(usage = "Usage: %prog -m meta -c cluster_column -o output -a assignment", option_list = options)
opt <- parse_args(opt_parser)

# TESTING
# opt = list(
#     meta = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/5_integration/2_Harmony_v2/Run_S18-5591-C8Us2/out/S18-5591-C8Us2_xenium_merged_UMAP_meta.tsv',
#     id_column = 'merged_id',
#     cluster_column = 'seurat_clusters',
#     output = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/5_integration/2_Harmony_v2/Run_S18-5591-C8Us2/out/4_assign_celltype/',
#     assignment = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/5_integration/2_Harmony_v2/Run_S18-5591-C8Us2/out/4_assign_celltype/1_celltype_seurat_clusters.tsv'
# )

print(opt)

# Parameters 
id_column = opt$id_column
cluster_column = opt$cluster_column
output_dir = opt$output

# Read meta data
# Load
meta = read_tsv(opt$meta) %>% select(all_of(c(id_column, cluster_column)))
assignment = read_tsv(opt$assignment) %>% .[,1:2] %>% setNames(c(cluster_column,'Celltype'))

# Report 
print(paste('Meta data:', dim(meta)))
print(paste('Assignment:', dim(assignment)))

print(paste('Unique clusters in meta:', length(unique(meta$seurat_clusters))))
print(paste('Unique clusters in assignment:', length(unique(assignment$seurat_clusters))))
if(length(unique(meta$seurat_clusters)) != length(unique(assignment$seurat_clusters))){
    warning('Number of clusters in meta and assignment do not match! Please check if using the correct files')
}

# Left join
meta_joined = left_join(meta, assignment, by = cluster_column); print(head(meta_joined))

# report NA after assignment
print(paste('NA after assignment:', sum(is.na(meta_joined$Celltype))))

# save to output
write_tsv(meta_joined, paste0(output_dir, '/2_celltype_assignment.tsv'))

# Do some stats
# 3.1 Count number of cell per CellType
cluster_count = assignment %>% group_by(Celltype) %>% 
    summarise(clusters = toString(.data[[cluster_column]]), n_cluster = n()) %>% 
    arrange(desc(n_cluster))
cell_count = meta_joined %>% group_by(Celltype) %>% summarise(n_cell = n()) %>% arrange(desc(n_cell)) 
both_count = left_join(cluster_count, cell_count, by = 'Celltype')
write_tsv(both_count, paste0(output_dir, '/3_celltype_count.tsv'))

