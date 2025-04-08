message("Unwrap to plot to list or add new plots to wrapped plot")
message('WrappedPlotToList(plot)')
message('AddPlotToWrappedPlot(wrapped_plot, plot_add, rewrap = T)')
library(patchwork)
WrappedPlotToList = function(plot){
  if(class(plot)[[1]] == 'patchwork'){
    n_plots = length(plot)
    # return list(plot[[1]], plot[[2]], ..., plot[[n_plots]])
    return(purrr::map(1:n_plots, ~plot[[.x]]))
  }else{
    return(plot)
  }
}

AddPlotToWrappedPlot = function(wrapped_plot, plot_add, placement = c('back','front'), rewrap = T, ...){
  unwrapped_list = WrappedPlotToList(wrapped_plot)
  # check if plot is list
  if(class(plot_add)[[1]] != 'list'){
    plot_add = list(plot_add)
  }
  # combine based on placement
  placement = match.arg(placement)
  if(placement == 'back'){
	unwrapped_list = c(unwrapped_list, plot_add)
	  }else{
	unwrapped_list = c(plot_add, unwrapped_list)
	}	
  
  # if rewrap, return as wrapplot
  if(rewrap){
    return(wrap_plots(unwrapped_list, ...))
  }else{
    return(unwrapped_list)
  }
}
