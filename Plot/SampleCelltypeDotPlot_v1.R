library(patchwork)
library(ggplot2)

message("SampleCelltypeDotPlot(meta_df, ident=NULL, feature = NULL, color_by = NULL)")
SampleCelltypeDotPlot = function(meta_df, ident=NULL, feature = NULL, color_by = NULL, 
	color_cell_fill = "#78cfa1", color_sample_fill = "#ec9f5c"){
	if(any(! c(ident, feature, color_by) %in% colnames(meta_df))){
		stop('ident, feature, color_by must be in count_df')
	}
	# If no color provided, created dummy one
	if(is.null(color_by)){
		meta_df = meta_df %>% mutate(color_column = 'dummy')
		color_by = 'color_column'
	}
	# Count 
	count_df = meta_df %>% count(.data[[ident]], .data[[feature]], .data[[color_by]])

	p_dot2 = count_df %>% ggplot(aes(x = .data[[ident]], y = .data[[feature]], size = n, color = .data[[color_by]])) + 
	geom_point(stat = 'identity', position = 'dodge', alpha = 0.6) + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
	labs(title = str_c(feature,' count per ', ident), x = 'Sample ID', y = 'Cell type', fill = 'Count') + 
	theme(legend.position = 'bottom')
	
	# part2a: Cell count bar
	p_bar1 = count_df %>% group_by(.data[[feature]]) %>% summarize(n = sum(n)) %>%
		mutate(n_show = sprintf("%.2g", n)) %>% 
		ggplot(aes(x = n, y = .data[[feature]], label = n_show)) +
		geom_bar(stat = 'identity', fill = color_cell_fill) +
		#geom_text(position = position_stack(vjust = 0))+
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
		theme(legend.position = 'bottom') +
		labs(x = 'Total Cells') + 
		# remove y axis
		theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
	# part2 : Add side bar plot count how many sample have given cell type. Fill by .data[[color_by]]
	p_bar2 = count_df %>% count(.data[[feature]], .data[[color_by]]) %>%
		ggplot(aes(y = .data[[feature]], x = n, fill = .data[[color_by]], label = n)) +
		geom_bar(stat = 'identity') +
		geom_text(position = position_stack(vjust = 0.5))+
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
		theme(legend.position = 'bottom') +
		labs(x = str_c("Per ", color_by)) + 
		# remove y axis
		theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

	# part3 : total count
	p_bar_total = count_df %>% count(.data[[feature]]) %>%
		ggplot(aes(x = n, y = .data[[feature]], label = n)) +
		geom_bar(stat = 'identity', fill = color_sample_fill) +
		geom_text(position = position_stack(vjust = 0.5))+
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
		theme(legend.position = 'bottom') +
		labs(x = str_c('Total ', ident)) + 
		# remove y axis
		theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

	p_final = wrap_plots(p_dot2, p_bar1, p_bar2, p_bar_total, widths = c(5,1, 1,1))
	return(p_final)
}