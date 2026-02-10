# Assign cell type to barcode
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(Seurat)
library(optparse)
library(qs)


# Take meta data, celltype assignment file, id_column, cluster_column, and output path
options <- list(
    make_option(c("-m", "--meta"), type="character", default="",
        help="Input seurat meta. If provided will not use the object since much faster "),
    make_option(c("-j", "--object"), type="character", default=NULL,
        help="Input seurat object"),
    make_option(c("-i", "--id_column"), type="character", default=NULL,
        help="(Deprecating) Column name in meta data that contains barcode information"),
    make_option(c("-c", "--cluster_column"), type="character", default=NULL,
        help="Column name in meta data that contains cluster information"),
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="Output directory"),
    make_option(c("-a", "--assignment"), type="character", default=NULL,
        help="Cell type assignment file")
)   

opt_parser <- OptionParser(usage = "Usage: %prog -j object -c cluster_column -o output -a assignment", option_list = options)
opt <- parse_args(opt_parser)



# # TESTING
# opt = list(
#     object = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/5_integration/2_Harmony_v2/Run_S18-15142-B17/out/S18-15142-B17_xenium_merged_UMAP.qs',
#     id_column = 'merged_id',
#     cluster_column = 'seurat_clusters',
#     output = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/5_integration/2_Harmony_v2/Run_S18-15142-B17/out/4_assign_celltype/',
#     assignment = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/Xenium/5_integration/2_Harmony_v2/Run_S18-15142-B17/out/4_assign_celltype/1_celltype_seurat_clusters.tsv'
# )

print(dput(opt))

# Parameters 
# id_column = opt$id_column # Deprecating not useing this
cluster_column = opt$cluster_column
output_dir = opt$output

# Read meta data
if(file.exists(opt$meta)){
    message('Reading meta')
    meta = read_tsv(opt$meta) 

    # turn first fow into cell_id
    meta = if(!'cell_id' %in% names(meta)){
        meta <- meta %>% mutate(cell_id = meta[[1]])
    }

    meta = meta %>% 
        select(all_of(c('cell_id', cluster_column))) %>%
        mutate(across(all_of(cluster_column), as.character))
    
} else if(!is.null(opt$object)){
    message('Loading object')
    # Load
    source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Loading/load_object_v1.R')
    obj = LoadObject(opt$object, nthreads = 50)
    # Get meta data
    meta = if(!'cell_id' %in% names(obj@meta.data)) obj@meta.data %>% rownames_to_column('cell_id') else obj@meta.data
    meta = meta %>% select(all_of(c('cell_id', cluster_column))) %>% 
        mutate(across(all_of(cluster_column), as.character))

} else {
    stop('Please provide meta or object ')
}


# Load assignment
assignment = read_tsv(opt$assignment) %>% .[,1:2] %>% setNames(c(cluster_column,'Celltype')) %>% 
    mutate(across(all_of(cluster_column), as.character))


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


# 4. make a version that contains sample id and individual barcode per sample
meta_joined_individual <- meta_joined %>% 
    separate(cell_id, into = c('sample_id', 'sample_cell_id'), sep = '__', remove = F) 
write_tsv(meta_joined_individual, paste0(output_dir, '/4_celltype_assignment_individual.tsv'))