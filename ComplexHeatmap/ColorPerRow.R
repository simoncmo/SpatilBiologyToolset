library(ComplexHeatmap)
library(circlize)
# 2025-09-08 Simon
# Function to create layer_fun for binary heatmap with custom colors
# 
# This function creates a layer function for ComplexHeatmap that colors binary matrix cells
# based on row names (e.g., cell types) using a custom color mapping. Cells with value 1
# are colored according to the mapping, while cells with value 0 remain white.
#
# @param binary_matrix A binary matrix (0/1) with named rows
# @param color_mapping A named vector where names are row names and values are colors
# @return A layer function for use in ComplexHeatmap::Heatmap()
#
# Example usage:
# binary_mat <- matrix(sample(0:1, replace = T, size = 25), nrow=5, dimnames=list(c("TypeA","TypeB","TypeC","TypeD","TypeE"), NULL))
# colors <- c("TypeA"="#E74C3C", "TypeB"="#3498DB", "TypeC"="#2ECC71", "TypeD"="#F39C12", "TypeE"="#9B59B6")
# layer_fun <- ColorPerRow(binary_mat, colors)
# Heatmap(binary_mat, layer_fun = layer_fun, show_heatmap_legend = FALSE)
ColorPerRow <- function(binary_matrix, color_mapping) {
	function(j, i, x, y, w, h, fill) {
		library(grid)
		# Get values and row indices
		v = pindex(binary_matrix, i, j)
		row_names_vec = rownames(binary_matrix)
		
		# Color based on row (luminal type)
		col = sapply(seq_along(v), function(idx) {
			row_name = row_names_vec[i[idx]]
			if(v[idx] == 1) {
				return(color_mapping[row_name])
			} else {
				return("white")
			}
		})
		
		grid.rect(x, y, w, h, gp = gpar(fill = col, col = "gray80"))
	}
}