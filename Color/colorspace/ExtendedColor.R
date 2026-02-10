source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Color/RcolorBrewer/GetAllQualColors.R')
library(colorspace)
# 2024-05-03 
# Check if a color is gray
is_gray = function(color_hex){
	rgb = col2rgb(color_hex)
	if(rgb[1] == rgb[2] & rgb[2] == rgb[3]) return(TRUE) else return(FALSE)
}
#is_gray('#0F0F0F')

# - add exclude_color argument to multiple functions
# - modifed DistinctColorListFromTable function to prevent color reuse
ExtendedColor = function(n = NULL, col_orig = GetAllQualColors(), exclude_color = NULL, no_gray = TRUE){
    col_light = lighten(col_orig, 0.2)
    col_darker = darken(col_orig, 0.2)
    col_use = c(col_orig, col_light, col_darker)
    # exclude gray
    if(no_gray) col_use = col_use[!sapply(col_use, is_gray)]
    # exclude
    col_use = setdiff(col_use, exclude_color)
    if(!is.null(n)) col_use = col_use[1:n]
    return(col_use)
}

AssignDistinctColor = function(ident_vector, exclude_color = NULL){
    col = ExtendedColor(length(unique(ident_vector)), exclude_color = exclude_color) 
    setNames(col, unique(ident_vector))
}

# V1 original version
# DistinctColorListFromTable = function(df){
#     columns_use = colnames(df)
#     map(columns_use, ~AssignDistinctColor(df[[.x]])) %>% setNames(columns_use)
# }

# 2024-05-03
# V2 assign color for all ident first then distribute to avoid same color
DistinctColorListFromTable <- function(df, assigned_colors = NULL, exclude_color = NULL) {
  columns_use <- colnames(df)
  ident_list <- purrr::map(columns_use, ~as.character(unique(df[[.x]])))
  flat_ident <- purrr::reduce(ident_list, union)
  ident_need_colors <- setdiff(flat_ident, names(assigned_colors))
  color_map <- c(AssignDistinctColor(ident_need_colors, exclude_color = c(unlist(assigned_colors), exclude_color)), assigned_colors)
  
  result_list <- purrr::map(ident_list, ~{
    colors <- sapply(.x, function(cell) if (cell %in% names(color_map)) color_map[cell] else NA)
    names(colors) <- .x  # Assign names to the resulting vector
    return(colors)
  })
  
  names(result_list) <- columns_use  # Assign names to the list elements
  return(result_list)
}

# Test the function
#result_list <- DistinctColorListFromTable(df, assigned_colors)



# # TEST
# test_df = tumor_gene_labeled_df %>% 
#     distinct(cluster, genelabel) %>% 
#     mutate(text = sample(c("",'A','B'), n(), replace = T))
# test_df
# DistinctColorListFromTable(test_df)


# Add color
# Given some palette and new identity, Add color to ones doesn't have color
FillColors = function(ident_vector, existing_palette = NULL, remove_not_exists = T){
    missing_ident = setdiff(ident_vector, names(existing_palette))
    # Keep color that used yet
    unused_color = setdiff(ExtendedColor(), existing_palette)
    missing_ident_color = setNames(unused_color[1:length(missing_ident)], missing_ident)
    # make final
    final_color = c(existing_palette, missing_ident_color)
    # Remove idents not in target ident
    if(remove_not_exists) final_color = final_color[names(final_color) %in% ident_vector]
    return(final_color)
}
# 2024-05-03 
# Check if a color is gray
is_gray = function(color_hex){
	rgb = col2rgb(color_hex)
	if(rgb[1] == rgb[2] & rgb[2] == rgb[3]) return(TRUE) else return(FALSE)
}
#is_gray('#0F0F0F')


# 2024-05-17 
# scPalete function from Cellchat
# From: https://github.com/sqjin/CellChat/blob/master/R/visualization.R
#' ggplot theme in CellChat
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 theme_classic element_rect theme element_blank element_line element_text
CellChat_theme_opts <- function() {
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme_classic() +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(color = "black")) +
    theme(axis.line.y = element_line(color = "black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank()) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}


## Other palettes


#' Generate ggplot2 colors
#'
#' @param n number of colors to generate
#' @importFrom grDevices hcl
#' @export
#'
ggPalette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Generate colors from a customed color palette
#'
#' @param n number of colors
#'
#' @return A color palette for plotting
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}

AssignDistinctScPalette = function(ident_vector, exclude_color = NULL){
  # generate a color palette
    col = scPalette(length(unique(ident_vector)))
  # exclude color, rerun with colorRampPalette
  if(!is.null(exclude_color)){
    col_remain = setdiff(col, exclude_color)
    col = grDevices::colorRampPalette(col_remain)(length(unique(ident_vector)))
  }
  # set names
    setNames(col, unique(ident_vector))
}

message("Loaded ExtendedColor(n = NULL, col_orig = GetAllQualColors())")
message("Loaded AssignDistinctColor(ident_vector, exclude_color = NULL)")
message("Loaded DistinctColorListFromTable(df)")
message("FillColors(ident_vector, existing_palette)")
message("is_gray(color_hex)")
message("scPalette(n)")
message("CellChat_theme_opts()")
message("ggPalette(n)")
message("scPalette(n)")
message("AssignDistinctScPalette(ident_vector, exclude_color = NULL)")

