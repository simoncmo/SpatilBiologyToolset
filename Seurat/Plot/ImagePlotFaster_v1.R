# conda activate seurat5_v2
# 2024-05-02 Simon Mo
library(tidyverse)
require(scattermore)
library(patchwork)
library(ggplot2)
library(Seurat)
# library(ggrastr)
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Plot/FetchDataFov.R')
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Color/colorspace/ExtendedColor.R')

# Message for loaded functions
message("theme_black")
message("ImagePlotFaster(obj, group.by = NULL, features = NULL, max.cutoff = 'q90', min.cutoff = 'q0', feature_palette = NULL, dimension = 600, fov = 'fov', black_theme = TRUE, groupby_palette = list()")
message("ImagePlotFasterAllFOV(obj, group.by = NULL, features = NULL, max.cutoff = 'q90', min.cutoff = 'q0', feature_palette = NULL, dimension = 600, nrow = NULL, black_theme = TRUE, groupby_palette = list())")
message("GetFovNames(obj, ident_column = 'orig.ident',  fov = 'fov')")
message("GetAllFovNames(obj, ident_column = 'orig.ident')")
message("ReplaceFovNames(obj, ident_column = 'orig.ident')")
message("internal_image_plot(data, column_use, plot_type = c('ident','feature'), dimension = 2000, colors = NULL)")
message("internal_cut_values(values, cuts_max = 'q100', cuts_min = 'q0')")


# internal
# Set background black
theme_black = theme(
		plot.background = element_rect(fill = 'black'),
		panel.background = element_rect(fill = 'black'),
		# All text white
		plot.title = element_text(color = 'white'),
		plot.subtitle = element_text(color = 'white'),
		legend.text = element_text(color = 'white'),
		legend.title = element_text(color = 'white'),
	)

internal_image_plot = function(data, column_use, plot_type = c('ident','feature'), 
	plot_method = c('scattermore', 'ggrastr'),
	# scattermore parameter
	dimension = 2000, 
	pointsize = 2, 
	# ggrastr parameter
	raster_ptsize = 0.1,
	raster_dpi = 300,
	raster_scale = 1,
	# other parameters
	colors = NULL, label = F, repel = T, subtitle = NULL,
	split.by = NULL
){
	plot_type = match.arg(plot_type)
	# Get plot ratio
	yx_ratio = diff(range(data$y)) / diff(range(data$x))

	# Set up base data plot
	p = ggplot(data = data, aes(x = x, y = y, color = .data[[column_use]]))

	# Determin which method to plot
	plot_method = match.arg(plot_method)
	print(paste0("Plotting method: ", plot_method))
	if(plot_method == 'scattermore'){
		message("Using scattermore")
		# Get pixels ratio
		dimension_x = dimension
		dimension_y = dimension * yx_ratio
		#message("dimension_x: ", dimension_x, " dimension_y: ", dimension_y)
	
	 	p <- p + geom_scattermore(pointsize = pointsize, pixels = c(dimension_x,dimension_y), interpolate = TRUE) 
	} else {
		message("Using ggrastr")
		# Use ggrastr to rasterize the plot
		p <- p + ggrastr::rasterize(geom_point(size = raster_ptsize), dpi = raster_dpi, scale = raster_scale)
	} 
	# Set theme
	p <- p + theme_void() +
			labs(title = column_use, subtitle = subtitle) + 
			# set ratio
			#theme(aspect.ratio = yx_ratio) +
			coord_fixed()
			# set title bold and center
			theme(
				plot.title = element_text(face = "bold", hjust = 0.5),
				plot.subtitle = element_text(hjust = 0.5)
				) 
	
	# label with text 
	if(label){
		# Calculate center per identity
		center_position_data = data %>% group_by(.data[[column_use]]) %>% summarize(x = mean(x), y = mean(y))
		if(repel){
			message("Repel")
			# p = p + ggrepel::geom_label_repel(data = center_position_data, aes(label = .data[[column_use]], fill = .data[[column_use]]), size = 3, color = 'gray10', max.overlaps = Inf)
		 	# No fill just tex
			p = p + ggrepel::geom_text_repel(data = center_position_data, aes(label = .data[[column_use]]), size = 3, color = 'gray10', max.overlaps = Inf)
		}else{
			p = p + geom_text(data = center_position_data, aes(label = .data[[column_use]]), size = 3)
		}
	}
	
	# split.by:
	if(!is.null(split.by)){
		message("Split by: ", split.by)
		p = p + facet_wrap(as.formula(paste0('~', split.by)))
	}

	# Set color
	p = if(plot_type == 'ident'){
		p + scale_color_manual(values = colors)
	}else{
		p + scale_color_gradientn(colors = colors)
	}
}

internal_cut_values = function(values, cuts_max = 'q100', cuts_min = 'q0', remove_zero = T){
	message("Cutting values with max: ", cuts_max, " min: ", cuts_min)
	cut_quantile_max = str_remove(cuts_max, 'q') %>% as.numeric() %>% `/`(100)
	cut_quantile_min = str_remove(cuts_min, 'q') %>% as.numeric() %>% `/`(100)
	# remove zero if remove_zero is TRUE
	values_calculate = if(remove_zero) values[values != 0] else values
	message("values_calculate: ", length(values_calculate))
	# if no values, return values
	if(length(values_calculate) == 0){ message("All value are zero. No Cutting");return(values)}
	max_value = quantile(values_calculate, cut_quantile_max, na.rm = T)
	min_value = quantile(values_calculate, cut_quantile_min, na.rm = T)
	message("max_value: ", max_value, " min_value: ", min_value)
	if(max_value > 0) values[values > max_value] <- max_value else message("max_value is 0. Not cutting max")
	values[values < min_value] <- min_value
	return(values)
}


# Faster image plot
ImagePlotFaster = function(obj, group.by = NULL, features = NULL, split.by = NULL,
	max.cutoff = "q90", min.cutoff = "q0", 
	feature_palette = NULL,  fov = NULL,
	black_theme = TRUE,
	groupby_palette = list(),
	nrow = NULL,
	plot_method = c('scattermore', 'ggrastr'),
	# scattermore parameter
	dimension = 2000, 
	pointsize = 2, 
	# ggrastr parameter
	raster_ptsize = 0.05,
	raster_dpi = 300,
	raster_scale = 1,
	# other parameters
	label = F, repel = T,
	sort = T,
	flip_y = FALSE
	){

	# Get data
	fov = fov %||% Images(obj)[[1]]
	message(fov)
	meta_use = FetchDataFov(obj, c(group.by, features, split.by), fov_use = fov)
	
	# Check if plot group.by
	p_dim_list = list()
	if(any(group.by %in% colnames(meta_use))){
		message("group.by: ", toString(group.by))
		group.by = intersect(group.by, colnames(meta_use))
		message("plotting: ", toString(group.by))
		p_dim_list = map(group.by, function(group_by_use){
			# set color
			cols_use = FillColors(obj@meta.data[[group_by_use]], groupby_palette[[group_by_use]])
			# filter out NA values
			meta_use = meta_use %>% filter(!is.na(.data[[group_by_use]]))
			# make sure all groupby value are discrete
			meta_use = meta_use %>% mutate({{group_by_use}} := as.factor(.data[[group_by_use]]))
			internal_image_plot(meta_use, group_by_use, plot_type = 'ident',
				plot_method = plot_method,
				dimension = dimension, 
				pointsize = pointsize, 
				raster_ptsize = raster_ptsize,
				raster_dpi = raster_dpi,
				raster_scale = raster_scale,
				colors = cols_use, 
				label = label, repel = repel, subtitle = fov, split.by=split.by)
		})
	}

	# CHeck if plot features
	if(is.null(feature_palette) & black_theme) feature_palette =  viridis::inferno(n = 10, begin = 0.3, end = 1)
	if(is.null(feature_palette) & !black_theme) feature_palette =  RColorBrewer::brewer.pal(9, "BuPu")
	p_feature_list = list()
	if(any(features %in% colnames(meta_use))){
		message("features: ", toString(features))
		message("plotting: ", toString(features))
		p_feature_list = map(features, function(feature_use){
			# Cut balues 
			meta_use = meta_use %>% mutate({{feature_use}} := internal_cut_values(.data[[feature_use]], cuts_max = max.cutoff, cuts_min = min.cutoff))
			# filter out NA values
			meta_use = meta_use %>% filter(!is.na(.data[[feature_use]]))
			# sort values
			if(sort) meta_use = meta_use %>% arrange(desc({{feature_use}}))
			internal_image_plot(meta_use, feature_use, plot_type = 'feature', 
				plot_method = plot_method,
				dimension = dimension, 
				pointsize = pointsize, 
				raster_ptsize = raster_ptsize,
				raster_dpi = raster_dpi,
				raster_scale = raster_scale,
				colors = feature_palette, 
				subtitle = fov, split.by=split.by)
		})
	}
	
	p_all = c(p_dim_list, p_feature_list)
	
	# Combine plot
	p = wrap_plots(p_all, nrow = nrow)
	if(black_theme) p = p & theme_black
	# flip y
	if(flip_y) p <- p & scale_y_reverse()
	return(p)
}

# Plot All Fov
ImagePlotFasterAllFOV = function(obj, group.by = NULL, features = NULL, split.by = NULL,
	max.cutoff = "q90", min.cutoff = "q0", 
	feature_palette = NULL, dimension = 2000, 
	combine = T,
	pointsize = 2, label = FALSE,
	nrow = NULL, black_theme = TRUE,
	groupby_palette = list()
	){
	p_out = map(Images(obj), function(fov_use){
		p = ImagePlotFaster(obj, group.by = group.by, 
			features = features, split.by = split.by,
			max.cutoff = max.cutoff, 
			min.cutoff = min.cutoff, feature_palette = feature_palette, 
			dimension = dimension, fov = fov_use, 
			black_theme = black_theme, groupby_palette = groupby_palette,
			nrow = 1, # Always 1 row per sample
			pointsize = pointsize, label = label
			) 
		if(black_theme) p = p & theme_black
		return(p)
	})
	if(combine){
		p_out = wrap_plots(p_out, nrow = nrow)
		if(black_theme) p_out = p_out & theme_black
	}
	return(p_out)
}

# Plot All Features with same scale
ImagePlotFasterAllFOVSameScale = function(obj, feature = NULL, split.by = NULL,
	max.cutoff = "q90", min.cutoff = "q0", 
	feature_palette = NULL, dimension = 2000, 
	pointsize = 2, label = FALSE,
	nrow = NULL, black_theme = TRUE,
	groupby_palette = list()
	){
	# Only Support 1 feature
	if(length(feature) > 1){feature = feature[[1]]; message("Only support 1 feature. Using first feature: ", feature)}
	p_list = map(Images(obj), function(fov_use){
		ImagePlotFaster(obj, group.by = NULL, 
			features = feature,  split.by = split.by,
			max.cutoff = 'q100', min.cutoff = 'q0', # no cutting
			feature_palette = feature_palette, 
			dimension = dimension, fov = fov_use, 
			black_theme = black_theme, groupby_palette = groupby_palette,
			nrow = 1, # Always 1 row per sample
			pointsize = pointsize, label = label
			) 
	})
	# Get expension range with cut factor
	value_cut <- internal_cut_values(FetchData(obj, feature), cuts_max = max.cutoff, cuts_min = min.cutoff)
	p_all = wrap_plots(p_list, nrow = nrow) & scale_color_gradientn(limits = range(value_cut), colors = viridis::inferno(n = 10, begin = 0.3, end = 1))
	if(black_theme) p_all = p_all & theme_black
	return(p_all)
}


## FOV related functions

# Fov rename functions
GetFovNames = function(obj_func, ident_column = 'orig.ident',  fov = 'fov'){
	message("ident_column use : ", ident_column)
	message("fov:", fov)
	fov_obj = obj_func@images[[fov]]
	fov_cells = Cells(fov_obj)
	# Get meta
	obj_func@meta.data[fov_cells, ident_column] %>% unique 
}

GetAllFovNames = function(obj_func, ident_column = 'orig.ident'){
	map(Images(obj_func), ~GetFovNames(obj_func, ident_column = ident_column, fov = .x)) %>% 
		setNames(Images(obj_func)) %>% unlist
}


ReplaceFovNames = function(obj_func, ident_column = 'orig.ident'){
	all_fov_names = GetAllFovNames(obj_func, ident_column = ident_column)
	print(all_fov_names)
	obj_func@images = obj_func@images %>% 
		setNames(
			all_fov_names[names(.)]
		)
	return(obj_func)
}



### FOR UMAP
FetchDim <- function(obj_func, reduction_use = NULL){
	message("Available Reductions: ", toString(Reductions(obj_func)))
	reduction_use = reduction_use %||% str_subset(Reductions(obj_func), 'umap')[[1]]
	message("Using: ", reduction_use)
	# Get data
	coord_df <- obj_func@reductions[[reduction_use]]@cell.embeddings[,1:2] %>% as.data.frame %>% rownames_to_column('cell') %>% setNames(c('cell','x','y'))
	return(coord_df)
}
FetchDataDim <- function(obj_func, values = NULL, reduction_use = NULL){
	coord_df = FetchDim(obj_func, reduction_use = reduction_use)
	cells = coord_df$cell
	value_df = FetchData(obj_func, values, cells= cells) %>% rownames_to_column('cell')
	left_join(coord_df, value_df, by = 'cell')
}
# Faster image plot
DimPlotFaster = function(obj_func, group.by = NULL, features = NULL, split.by = NULL,
	max.cutoff = "q90", min.cutoff = "q0", 
	feature_palette = NULL, dimension = 2000, fov = NULL,
	black_theme = TRUE,
	groupby_palette = list(),
	nrow = NULL,
	pointsize = 2,
	label = F, repel = T,
	sort = T,
	reduction_use = NULL
	){
	# Get data
	meta_use = FetchDataDim(obj_func, c(group.by, features, split.by), reduction_use = reduction_use)
	
	# Check if plot group.by
	p_dim_list = list()
	if(any(group.by %in% colnames(meta_use))){
		message("group.by: ", toString(group.by))
		group.by = intersect(group.by, colnames(meta_use))
		message("plotting: ", toString(group.by))
		p_dim_list = map(group.by, function(group_by_use){
			# Make ident discrete
			meta_use = meta_use %>% mutate({{group_by_use}} := as.factor(.data[[group_by_use]]))
			# set color
			
			cols_use = FillColors(obj_func@meta.data[[group_by_use]], groupby_palette[[group_by_use]])
			print(cols_use)
			# filter out NA values
			meta_use = meta_use %>% filter(!is.na(.data[[group_by_use]]))
			internal_image_plot(meta_use, group_by_use, plot_type = 'ident',dimension = dimension, colors = cols_use, pointsize = pointsize, label = label, repel = repel, subtitle = fov, split.by=split.by)
		})
	}

	# CHeck if plot features
	if(is.null(feature_palette) & black_theme) feature_palette =  viridis::inferno(n = 10, begin = 0.3, end = 1)
	if(is.null(feature_palette) & !black_theme) feature_palette =  RColorBrewer::brewer.pal(9, "BuPu")
	p_feature_list = list()
	if(any(features %in% colnames(meta_use))){
		message("features: ", toString(features))
		message("plotting: ", toString(features))
		p_feature_list = map(features, function(feature_use){
			# Cut values 
			meta_use = meta_use %>% mutate({{feature_use}} := internal_cut_values(.data[[feature_use]], cuts_max = max.cutoff, cuts_min = min.cutoff))
			# filter out NA values
			meta_use = meta_use %>% filter(!is.na(.data[[feature_use]]))
			# sort values
			if(sort) meta_use = meta_use %>% arrange(desc({{feature_use}}))
			internal_image_plot(meta_use, feature_use, plot_type = 'feature', dimension = dimension, colors = feature_palette, pointsize = pointsize, subtitle = fov, split.by=split.by)
		})
	}
	
	p_all = c(p_dim_list, p_feature_list)
	
	# Combine plot
	p = wrap_plots(p_all, nrow = nrow)
	if(black_theme) p = p & theme_black
	return(p)
}


## For spatial And DimPlot
ImageDimComboPlotFaster <- function(obj_func, group.by = NULL, features = NULL, split.by = NULL,
	nrow = NULL, ncol = NULL, black_theme = TRUE, flip_y = FALSE, reduction_use = NULL, ...
	){
	# if nrow not set, decided by number of values to plot. >1 values, nrow = 2
	if(is.null(nrow)) nrow = ifelse(length(group.by) + length(features) > 1, 2, 1); message("nrow: ", nrow)
	p_dim = DimPlotFaster(obj_func, group.by = group.by, 
			features = features, split.by = split.by, nrow = 1, reduction_use=reduction_use,...
			)
	p_image = ImagePlotFaster(obj_func, group.by = group.by,
			features = features, split.by = split.by, nrow = 1, flip_y=flip_y, ...
			)
	p_all = wrap_plots(list(p_dim, p_image), nrow = nrow, ncol = ncol)
	if(black_theme) p_all = p_all & theme_black
	return(p_all)
}

