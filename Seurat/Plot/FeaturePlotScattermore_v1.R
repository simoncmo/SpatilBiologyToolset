library(scattermore)
message("FeaturePlotScattermore(obj, features, assay = 'RNA', slot = 'data', sort = T, reduction = 'umap', pt_size = 4, pix = 700, ...) is loaded.")
FeaturePlotScattermore = function(obj, features, assay = 'RNA', slot = 'data', sort = T, reduction = 'umap', pt_size = 4, pix = 700, ...){
  DefaultAssay(obj) <- assay
  # find reduction
  if(!reduction %in% Reductions(obj)){
    message('Reduction not found in object. Available reduction: ', Reductions(obj))
    reduction = str_detect(Reductions(obj), 'umap')[[1]]
    print(reduction)
  }
  # Get umap position 
  position = Embeddings(obj, reduction = reduction)
  x_col = colnames(position)[[1]]#; print(x_col)
  y_col = colnames(position)[[2]]#; print(y_col)

  # get Data
  df = FetchData(obj, features, layer = slot) %>% rownames_to_column('ids')
  # Add embedding
  df = cbind(df, position)
  
  # loop throught all the features
  p_all = map(features, function(feature_use){
    # sort 
    if(sort){
      df = df %>% arrange(desc(.data[[feature_use]]))
    }
    # ggplot 
    p = ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[feature_use]])) + 
      geom_scattermore(pointsize=pt_size,
    pixels=c(pix,pix)) + 
      cowplot::theme_cowplot() + 
      labs(title = feature_use) +
      theme(aspect.ratio = 1) + 
      #scale_color_viridis_c(option = 'F', begin = 0.3, end = 1, direction = -1) 
      scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, 'BuPu'))
      #scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, 'Reds'))
    return(p)
  }) %>% wrap_plots()
  return(p_all)
}