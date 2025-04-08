library(RColorBrewer)
GetAllQualColors = function(category_use = 'qual'){
    brewer.pal.info %>% 
    rownames_to_column('Pal') %>%
    filter(category == category_use) %>% 
    pmap(function(Pal, maxcolors, ...){
        RColorBrewer::brewer.pal(n = maxcolors, name = Pal)
    }) %>% unlist
}
