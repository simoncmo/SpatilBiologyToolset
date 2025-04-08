# Set heatmapcolor
library(circlize)
library(pals)
# GetColorRampZero: Function to get color ramp that have 0 in the middle
GetColorRampZero <- function(matrix, max_scale = 0.8, palette = pals::parula(4)){
    mtx_range = range(matrix) %>% `*`(max_scale)
    mtx_range = c(mtx_range[[1]], mtx_range[[1]]/2, 0, mtx_range[[2]]/2, mtx_range[[2]])
    palette_use = c(palette[1:2], 'white', palette[c(length(palette)-1,length(palette))])
    ramp_color = colorRamp2(mtx_range, palette_use)
    return(ramp_color)
}

# GetColorRamp: Function to get color ramp that does not have 0 in the middle
GetColorRamp <- function(matrix, max_scale = 0.8, palette = pals::parula(4)){
    mtx_range = range(matrix) %>% `*`(max_scale)
    mtx_range_seq = mtx_range %>% {seq(from = .[[1]], to = .[[2]], length.out = length(palette))}
    ramp_color = colorRamp2(mtx_range_seq, palette)
    return(ramp_color)
}