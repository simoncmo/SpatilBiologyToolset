library(tidyverse)
library(ggplot2)

## Create synthetic data for testing 

# samples = str_c('Sample', 1:10)
# celltypes = c('Tumor','Immune','Stromal', str_c('Celltype', 1:10))
# genes = c('EPCAM','AR','KRT17','ACTA2','CD3')


# # Make table
# library(dplyr)

# exp_df <- expand.grid(samples, celltypes, genes) %>% 
#   setNames(c('Sample','Celltype','Gene')) %>% 
#   group_by(Celltype, Gene) %>% 
#   mutate(
#     Exp = case_when(
#       Celltype == 'Tumor' & Gene %in% c('EPCAM','KRT17','AR') ~ rnorm(n(), mean = 10),
#       Celltype != 'Tumor' & Gene %in% c('EPCAM','KRT17','AR') ~ abs(rnorm(n(), mean = 0.5)),
#       Celltype == 'Immune' & Gene %in% c('CD3') ~ rnorm(n(), mean = 10),
#       Celltype != 'Immune' & Gene %in% c('CD3') ~ abs(rnorm(n(), mean = 0.6)),
#       Celltype == 'Stromal' & Gene %in% c('ACTA2') ~ rnorm(n(), mean = 10),
#       Celltype != 'Stromal' & Gene %in% c('ACTA2') ~ abs(rnorm(n(), mean = 0.6))
#     )
#   ) %>% 
#   mutate(pct1 = rnorm(n = n(), mean = 5))

message("DotViolinPlot(exp_df, sample_column = 'Sample', expresssion_column = 'Exp', celltype_column = 'Celltype', gene_column = 'Gene')")
library(patchwork)
DotViolinPlot = function(exp_df, sample_column = 'Sample', expresssion_column = 'Exp', 
	celltype_column = 'Celltype', gene_column = 'Gene',
	sample_facet = "."
	){

	# Make plot
	p_dot = ggplot(exp_df, aes(x = .data[[sample_column]], y= .data[[expresssion_column]], color = .data[[celltype_column]])) + 
		facet_grid(.data[[gene_column]]~., space = 'fixed', scales = 'free') + 
		geom_point() +
		theme_bw()

	# Violin
	p_vio = ggplot(exp_df, aes(x = .data[[celltype_column]], y = .data[[expresssion_column]], fill = .data[[celltype_column]])) + 
		facet_grid(.data[[gene_column]] ~., space = 'fixed', scales = 'free') + 
		geom_violin(scale = 'width' ) + 
		theme_bw() + 
		theme(
			axis.text = element_blank(),
			axis.ticks.y = element_blank()
		) + 
		labs(y = "")


	# determine width
	n_sample = length(unique(exp_df[[sample_column]])) 
	n_celltype = length(unique(exp_df[[celltype_column]]))

	wrap_plots(p_dot, p_vio, widths = c(n_sample, n_celltype)) & 
		theme(legend.position = 'bottom') &
		theme(
			strip.text.y = element_text(angle = 0),
			panel.grid.minor.y=element_blank(),
			panel.grid.major.y=element_blank(),
			panel.grid.major.x=element_line(linewidth=  0.1), 
			axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)
			)
}


# DotViolinPlot(exp_df)
