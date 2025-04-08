# 2024-05-03
# FUnction to get data.frame to heatmap
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Color/colorspace/ExtendedColor.R')

# Report function
message("summarizeByColumn(df, by_col = NULL)")
message("getHierarchicalOrder(mtx, k = 3, fill_na = T)")
message("TableToHeatmap(df, x='', y='', value_column = '', scale_by = c('None','row','column'), cluster_rows = TRUE, cluster_columns = TRUE, heatmap_title = NULL, top_columns=NULL, top_color_list = NULL, left_color_list = NULL, left_columns = NULL)")
message("AverageByGroup(meta_df, group = '', make_long = T, names_to ='name', values_to = 'value')")
message("MakeColorRamp(matrix, colors = c('#0066ff','white','#da5c2e'))")
message("TableToDummyHeatmap(ident_df, id_col = 'gene', ident_col = 'cluster', row_order = NULL)")
message("AverageExpressionHeatmap(obj, id_col = 'seruat_clusters', heatmap_colors = c('white','#ff3300'), id_order = NULL, meta_column = NULL, gene_use = NULL, gene_label_df = NULL, idents_only = NULL)")
message("AverageExpressionHeatmapFromDEG(obj, id_col = 'seurat_clusters', heatmap_colors = c('white','#ff3300'), meta_column = NULL, deg = NULL, gene_use = NULL, fc_cutoff = 1, p_val_cutoff = 0.01, top_n = 10, idents_only= NULL, row_title_rot = 0, cluster_columns = F, cluster_rows = F)")


library(ComplexHeatmap)

# summarize without dropping
summarizeByColumn = function(df, by_col = NULL, id_to_rownames = TRUE){
	# use first if by_col is NULL
	if(is.null(by_col)){
		by_col = colnames(df)[[1]]
	}
	df = df %>% group_by(.data[[by_col]]) %>% 
			summarize(across(everything(), ~toString(sort(unique(.x))))) 

	if(id_to_rownames) df = column_to_rownames(df, by_col) 
	return(df)
}

# FactorEveryColumn
FactorEveryColumn = function(df){
	rownames_of_df = rownames(df)
	df = df %>% mutate(across(everything(), ~factor(., levels = unique(.))))
	df = as.data.frame(df) 
	rownames(df) = rownames_of_df
	return(df)
}

# Get cluster
getHierarchicalOrder = function(mtx,  fill_na = T){
	# 0. fill na 
	if(fill_na){
		mtx[is.na(mtx)] = 0
	}
	# 0. run hc
	hc = hclust(dist(mtx))
	return(hc$labels[hc$order])
}

# # Test
# random_mtx = matrix(rnorm(100), nrow = 10)
# rownames(random_mtx) = paste0('row', 1:10)
# colnames(random_mtx) = paste0('col', 1:10)
# hc = hclust(dist(random_mtx))
# hc$labels[hc$order]

TableToHeatmap = function(df, x="", y="", 
	value_column = "", 
	heatmap_colors = c('#0066ff','white','#da5c2e'),
	scale_by = c('None',"row",'column'),
	cluster_rows = TRUE, cluster_columns = TRUE,
	heatmap_title = NULL,
	top_columns=NULL, 
	annotation_colors = NULL, 
	left_columns = NULL,
	column_order = NULL,
	row_order = NULL,
	split_row = F, split_col = F,
	...
	){
	# By default, use first column as x and second column as y, 3rd as value
	if(x == "") x = colnames(df)[[1]]
	if(y == "") y = colnames(df)[[2]]
	if(value_column == "") value_column = colnames(df)[[3]]

	# make sure values in df
	if(any(!c(x,y,value_column,top_columns,left_columns) %in% colnames(df))){
		missing_column = setdiff(c(x,y,value_column,top_columns,left_columns), colnames(df))
		message(paste("Missing columns: ", missing_column))
		stop("Not all columns are in the dataframe. Please double check!")
	}
	# 1. create matrix
	mtx = df %>% 
		pivot_wider(id_cols= {{y}}, names_from = {{x}} ,values_from = {{value_column}}) %>%
		column_to_rownames({{y}}) %>% as.matrix

	# 3a. title
	if(is.null(heatmap_title)){
		heatmap_title = str_glue('Heatmap of {value_column} by {x} and {y}')
	}
	# 3b. scale
	scale_by = match.arg(scale_by)
	if(scale_by == 'row'){
		mtx = t(scale(t(mtx)))
	} else if(scale_by == 'column'){
		mtx = scale(mtx)
	}
	# 3c. cluster
	# if has na in mtx and want to cluster, fill NA with 0 and pre calculate the order
	if(cluster_rows){
		if(any(is.na(mtx))){
			message("Found NA in matrix. Use external clustering function")
			cluster_rows = F
		}
		row_order = getHierarchicalOrder(mtx)
		mtx = mtx[row_order,]
	}

	if(cluster_columns){
		if(any(is.na(mtx))){
			message("Found NA in matrix. Use external clustering function")
			cluster_columns = F
		}
		col_order = getHierarchicalOrder(t(mtx))
		mtx = mtx[,col_order]
	}

	# 3cc. matrix order 
	if(!is.null(row_order)){
		print('change row order')
		mtx = mtx[row_order,]
		cluster_rows = F
	}
	if(!is.null(column_order)){
		print('change column order')
		mtx = mtx[,column_order]
		cluster_columns = F
	}
	print(rownames(mtx))
	print('mtx-2')
	#print(mtx[1:3,1:3])
	# note, annotation comes after mtx to prevent order difference
	# 3d. get top annotation
	color_used = c()
	if(!is.null(top_columns)){
		top_df = df[,c(x, top_columns)] %>% summarizeByColumn() %>% .[colnames(mtx), ,drop = F]
		# Make the top_columns factors so the heatmap doesn't change order on it own when adding annotation
		top_df = FactorEveryColumn(top_df)
		# 2a0. fill in any color list that's not provided 
		top_color_list = DistinctColorListFromTable(top_df, annotation_colors)
		# remove na in color_list
		top_color_list = map(top_color_list, ~.x[!is.na(.x)])
		# update color used
		color_used = c(color_used, reduce(top_color_list, union))
		# annotation
		message('hERERE')
		top_annot = HeatmapAnnotation(df = top_df, which = 'column', col = top_color_list)
	}
	# 3e. get left annotation
	if(!is.null(left_columns)){
		left_df = df[,c(y, left_columns)] %>% summarizeByColumn() %>% .[rownames(mtx), ,drop = F]
		# Make the left_columns factors so the heatmap doesn't change order on it own when adding annotation 
		left_df = FactorEveryColumn(left_df)
		# 2b0. fill in any color list that's not provided 
		print(left_df)
		left_color_list = DistinctColorListFromTable(left_df, c(annotation_colors, color_used))
		# remove na in color_list
		left_color_list = map(left_color_list, ~.x[!is.na(.x)])
		# annotation
		left_annot = HeatmapAnnotation(df = left_df, which = 'row', col = left_color_list)
	}
	# 4. make heatmap
	
	p_hm = Heatmap(
		mtx,
		name = {{value_column}},
		col = MakeColorRamp(mtx, heatmap_colors), 
		cluster_rows = cluster_rows,
		cluster_columns = cluster_columns,
		column_title=heatmap_title,
		# make bold and bigger
		column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
		top_annotation = if(!is.null(top_columns)) top_annot else NULL,
		left_annotation = if(!is.null(left_columns)) left_annot else NULL,
		# splits
		row_split = if(split_row) left_df else NULL, # Use first column of left annotation as split
		column_split = if(split_col) top_df else NULL, # Use first column as top annotation as split
		# remove space betwee slices
		row_gap = unit(0, "mm"),
		column_gap = unit(0, "mm"),
		# color
		na_col = 'gray90',
		...
	)
	# 6. return
	return(p_hm)

	#return(p_hm)
}

# For x, y not unique. summarize them first
# return columns: group, names_to, values_to
AverageByGroup = function(meta_df, group = "", make_long = T, names_to ='name', values_to = 'value'){
	if(group == "") group = colnames(meta_df)[[1]]
	df = meta_df %>% 
		group_by(.data[[group]]) %>% 
		mutate(across(where(is.factor), as.character)) %>%
		summarize(
			across(where(is.numeric), ~mean(.x, na.rm = T)),
			across(where(is.character), ~toString(sort(unique(.x))))
		) %>% 
		ungroup
	if(make_long){
		df = df %>% 
			pivot_longer(-{{group}}, names_to = names_to, values_to = values_to)
	}
	return(df)
}

# Get gradient color for heatmap
library(circlize)
MakeColorRamp = function(matrix, colors = c('#0066ff','white','#da5c2e')){
	# get breaks based on number of colors
	n_breaks = length(colors)
	quantiles_use = seq(0,1,length.out = n_breaks)
	breaks = quantile(matrix, quantiles_use, na.rm=T)
	color_fun = colorRamp2(breaks, colors)
	return(color_fun)
}

# With ID and multiple-idents, plot the dummy heatmap
# e.g. ID: DEG, idents = Celtype1, Celltyp2 ....
TableToDummyHeatmap = function(ident_df, 
	id_col = 'gene', ident_col = 'cluster', row_order = NULL){
	# 1. make dummy matrix
	dummy_mtx = table(ident_df[[id_col]], ident_df[[ident_col]])
	# 2. reorder
	if(!is.null(row_order)) dummy_mtx = dummy_mtx[row_order,]
	# 3. make heatmap
	p_hm = Heatmap(
		dummy_mtx,
		col = colorRamp2(c(0, max(dummy_mtx)), c('white', '#23ade9')),
		name = 'count',
		cluster_rows = T,
		cluster_columns = T,
		column_title='Membership Heatmap',
		# make bold and bigger
		column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
		# color
		na_col = 'gray90'
	)
	return(p_hm)
}

# # TEST sampple
# n_test_row = 20
# test_df = data.frame(
# 	gene = sample(str_c('gene', 1:5), n_test_row, replace = T),
# 	cluster = sample(str_c('celltype', 1:5), n_test_row, replace = T)
# )

# TableToDummyHeatmap(test_df)




AverageExpressionHeatmap = function(obj, 
	id_col = 'seruat_clusters', 
	heatmap_colors = c('white','#ff3300'),
	id_order = NULL, # if ident have predecided order, added here
	meta_column = NULL, gene_use = NULL, gene_label_df = NULL,
	idents_only = NULL,
	...
	){
		if(is.null(gene_use)){
			gene_use = sample(rownames(obj),20)
		}
		# Get data
		gene_use = unique(gene_use)
		message("Fetching data ...")
		data_df = FetchData(obj, vars = c(id_col,gene_use))
		if(!is.null(idents_only)) data_df = data_df %>% filter(.data[[id_col]] %in% idents_only)

		# get average expression
		message("Calculating expression ...")
		exp_df = AverageByGroup(data_df, group = id_col, names_to ='gene', values_to = 'Expression')
		# return columns: id_col, gene, Expression
		#print(names(exp_df))
		
		# reorder ident
		message('id_order')
		print(id_order)
		if(!is.null(id_order)){
			exp_df = exp_df %>% 
				mutate({{id_col}}:= factor(.data[[id_col]], levels = id_order)) %>% 
				arrange({{id_col}})
			print(unique(exp_df[[id_col]]))
		}

		# add meta
		message("Fetching meta ...")
		#print(c(id_col,meta_column))
		meta_distinct = obj@meta.data[,c(id_col,meta_column), drop = F] %>% 
			AverageByGroup(names_to='Expression', make_long = F)
		
		# Add gene label
		if(!is.null(gene_label_df)){
			message("Adding gene label ...")
			#print(gene_label_df %>% head())
			exp_df = left_join(exp_df, gene_label_df, by = 'gene')
		}

		# plot
		message("Plotting ...")
		plt_df = left_join(exp_df, meta_distinct, by = id_col)  
		#return(plt_df)
		TableToHeatmap(plt_df, scale_by='row', 
			top_columns = meta_column, 
			left_columns = if(is.null(gene_label_df)) NULL else setdiff(names(gene_label_df), 'gene'),
			heatmap_colors = heatmap_colors,
			...
			)
	}

AverageExpressionHeatmapFromDEG = function(obj, 
	id_col = 'seurat_clusters', 
	heatmap_colors = c('white','#ff3300'),
	meta_column = NULL, deg = NULL,
	gene_use = NULL,
	fc_cutoff = 1, p_val_cutoff = 0.01, top_n = 10,
	idents_only= NULL,
	row_title_rot = 0,
	cluster_columns = F, cluster_rows = F,
	 ...
	){
		if(is.null(deg) | !all(c('cluster','gene') %in% colnames(deg))){
			stop("Please provide deg table from FindMarkerAll with cluster and gene column")
		}
		# filter deg
		filter_df = deg %>% 
			filter(avg_log2FC > fc_cutoff, p_val_adj < p_val_cutoff) %>% 
			group_by(cluster) %>% 
			slice_max(n = top_n, order_by = avg_log2FC) %>%
			select(cluster,gene) %>% 
			arrange(cluster)
		
		if(!is.null(idents_only)) filter_df = filter_df %>% filter(cluster %in% idents_only)
		message("After filtering, there are ", nrow(filter_df), " genes left")

		# Get cluster order to use in next function
		# THIS WILL CAUSE ISSUE SINCE cluster will be removed from filtering 
		# cluster_order = filter_df$cluster %>% unique
		#return(filter_df)
		# Aggregate
		deg_label_df = summarizeByColumn(filter_df[,c('gene','cluster')], by_col = 'gene', id_to_rownames = F)
		#return(deg_label_df)
		# select gene
		if(is.null(gene_use)) gene_use = filter_df$gene %>% unique
		#return(filter_df)
		# # 4a. deg membership heatmap
		# deg_hm = TableToDummyHeatmap(filter_df, id_col = 'gene', ident_col = 'cluster')
		# deg_mtx = deg_hm@matrix

		# Plot
		exp_hm = AverageExpressionHeatmap(obj, 
			id_col = id_col, 
			id_order = NULL, # cluster_order,
			meta_column = meta_column, 
			heatmap_colors = heatmap_colors,
			gene_use = gene_use, 
			gene_label_df = deg_label_df,
			split_row = T,
			idents_only = idents_only, 
			row_title_rot = row_title_rot,
			cluster_columns = cluster_columns, cluster_rows= cluster_rows,
			#column_order = colnames(deg_mtx),
			#row_order = rownames(deg_mtx),
			...
			)
		
		return(exp_hm)
	}

