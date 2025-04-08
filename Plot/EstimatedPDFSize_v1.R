message("EstimateHeatmapPDFsize(idents, n_gene_per_ident = 5, 
	width_base = 5, height_base = 4,
	width_ratio = 3, height_ratio = 5, n_idents = NULL, n_gene = NULL) is loaded.")
message("EstimateHeatmapPDFsizeFromMatrix(matrix)")
message("calc_ht_size(ht, unit = 'inch')")

# Function to estimate pdf size
EstimateHeatmapPDFsize = function(idents, n_gene_per_ident = 5, 
	width_base = 5, height_base = 4,
	width_ratio = 3, height_ratio = 5, n_idents = NULL, n_gene = NULL){
	# Assume width for ident, height for gene
	# heatmap annotation on right
	if(is.null(n_idents)) n_idents = length(unique(idents)) 
	# if n_gene provided directly use that to calculate the height. esle calculate based on idents * n_gene_per_ident
	if(is.null(n_gene)) n_gene = n_idents * n_gene_per_ident
	return(
		list(
			width = width_base + n_idents / width_ratio,
			height = height_base + n_gene / height_ratio
		)
	)
}


# Function to estimate pdf size
EstimateHeatmapPDFsizeFromMatrix = function(matrix, 
	width_base = 5, height_base = 2,
	width_ratio = 6, height_ratio = 6, n_idents = NULL, n_gene = NULL){
	# Assume width for ident, height for gene
	# heatmap annotation on right
	if(is.null(n_idents)) n_idents = nrow(matrix)
	# if n_gene provided directly use that to calculate the height. esle calculate based on idents * n_gene_per_ident
	if(is.null(n_gene)) n_gene = ncol(matrix)
	return(
		list(
			width = width_base + n_idents / width_ratio,
			height = height_base + n_gene / height_ratio
		)
	)
}


# From the author
# https://jokergoo.github.io/2020/05/11/set-cell-width/height-in-the-heatmap/
calc_ht_size = function(ht, unit = "inch") {
    pdf(NULL)
    ht = draw(ht)
    w = ComplexHeatmap:::width(ht)
    w = convertX(w, unit, valueOnly = TRUE)
    h = ComplexHeatmap:::height(ht)
    h = convertY(h, unit, valueOnly = TRUE)
    dev.off()

    c(w, h)
}