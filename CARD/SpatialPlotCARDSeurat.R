library(Seurat)
library(scattermore)

message("SpatialPlotCARDSeurat(seurat_card_obj, features = feature_plt, pt_size = 4, max_cutoff = 0.75, pixels = 1000) is loaded.")
SpatialPlotCARDSeurat = function(seurat_card_obj, features = feature_plt, slot = 'data',
	pt_size = 4,
	max_cutoff = 0.85, pixels = 1000){
	exp_df = FetchData(seurat_card_obj, c(feature_plt,'sn_x','sn_y'), layer = 'data') 
	exp_max = exp_df[[feature_plt]] %>% max %>% `*`(., max_cutoff); print(exp_max)
	p_sn_spatial = ggplot(exp_df, aes(x = sn_x, y = -sn_y, color = .data[[feature_plt]])) + 
		geom_scattermore(pointsize=pt_size,
		pixels=c(pixels,pixels)
		) + 
		scale_color_viridis_c(
			option = 'A', 
			limits = c(0, exp_max), 
			oob = scales::squish) +
		coord_fixed() + 
		theme_void()
	return(p_sn_spatial)
}