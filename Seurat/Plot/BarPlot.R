source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Color/colorspace/ExtendedColor.R')
library(Seurat)
library(patchwork)
message("BarPlot(obj, group.by = 'seurat_clusters', fill.by = 'sample_id', plot = c('both','fill','stack'))")
BarPlot = function(obj, group.by = 'seurat_clusters', fill.by = 'sample_id', plot = c('both','fill','stack'), meta_df = NULL, filter_fillby_na = FALSE){
    # meta df
    if(is.null(meta_df)) meta_df = obj@meta.data

    # filter filter_fillby_na
    if(filter_fillby_na){
        message('Filtering out NA in fill.by')
        meta_df = meta_df %>% filter(!is.na(.data[[fill.by]]))
    }
    
    # Count group.by to sort
    group_by_order = meta_df %>% count(.data[[group.by]]) %>% arrange(desc(n)) %>% pull(.data[[group.by]]) %>% as.character
    
    # Reoder 
    df = meta_df
    df[[group.by]] = factor(df[[group.by]], levels = group_by_order)
    
    # Get Color
    n_ident_fillby = meta_df[[fill.by]] %>% unique %>% length
    color_fillby =  ExtendedColor(n_ident_fillby)

    # Make fill.by factor as well
    df[[fill.by]] = factor(df[[fill.by]], levels = unique(df[[fill.by]]))
    
    # Choose plot type
    plot = match.arg(plot)
    # Make fill plot
    if(plot %in% c('both', 'fill')){
        p_fill = df %>% 
            ggplot(aes(x = .data[[group.by]], fill = .data[[fill.by]])) + 
            geom_bar(position = 'fill') +
            cowplot::theme_cowplot() + 
            labs(x = group.by, y = 'Proportion') +
            scale_fill_manual(values = color_fillby) + 
            theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
            ) 
    }
    if(plot %in% c('both', 'stack')){
        p_stack = df %>% 
            ggplot(aes(x = .data[[group.by]], fill = .data[[fill.by]])) + 
            geom_bar(position = 'stack') +
            cowplot::theme_cowplot() + 
            labs(x = group.by, y = 'Count') +
            scale_fill_manual(values = color_fillby) + 
            theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
            ) 
    }
    # Return plots
    if(plot == 'both'){
        return(wrap_plots(p_stack, p_fill, ncol = 1, guides = 'collect') + plot_annotation(title = str_glue("Barplot by {group.by} and {fill.by}")))
    }else if(plot == 'fill'){
        return(p_fill + labs(title = str_glue("Barplot by {group.by} and {fill.by}")))
    }else if(plot == 'stack'){
        return(p_stack + labs(title = str_glue("Barplot by {group.by} and {fill.by}")))
    }
    
}   
