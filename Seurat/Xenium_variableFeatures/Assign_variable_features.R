message("AssignVariableFeatures(obj, assay_use = NULL, features)")
message("Assign377GenesAsVariableFeatures(obj, assay_use = NULL)")

AssignVariableFeatures = function(obj, assay_use = NULL, features){
  if(is.null(assay_use)) assay_use = DefaultAssays(obj)
  message(assay_use)
  # Assay3 method
  if(class(obj@assays[[assay_use]]) == 'Assay'){
    message(assay_use, "is Assay3")
    obj@assays[[assay_use]]@var.features = features
  }

  # Assay5 method:
  if(class(obj@assays[[assay_use]]) == 'Assay5'){
    message(assay_use, "is Assay5")
    obj@assays[[assay_use]][[]]$var.features = features
  }
  return(obj)
}

source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Xenium_panel/LoadPanel.R')
Assign377GenesAsVariableFeatures = function(obj, assay_use = NULL){
  if(is.null(assay_use)) assay_use = DefaultAssays(obj)
  message("Skip FindVariableFeatures. Assign the shared 377 genes as variable features")
  message("Loading shared genes");print(Sys.time())
  shared_genes = LoadXeniumhMultiGenes(); print(length(shared_genes))
  
  # Even faster way to set variable features
  message("Setting variable features");print(Sys.time())
  shared_genes = intersect(rownames(obj), shared_genes); print(length(shared_genes))
  
   # Assign
  obj = AssignVariableFeatures(obj, assay_use = assay_use, features = shared_genes)
  return(obj)
}