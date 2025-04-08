message("ExtractCoordinateFromCARDobject(CARDobject) is loaded.")

ExtractCoordinateFromCARDobject = function(SCEobject){
	# Note, assume coordinate store in colnames as follow:
	# "AGAACCGCACCCTCAC-1:2159x7464:2152.82024847367x7423.57504707975"
	str_split(colnames(SCEobject), ":", 3, simplify = T) %>% 
		as.data.frame %>% 
		setNames(c('sn_Barcode','ST_coordinates', 'sn_coordinates')) %>% 
		separate(ST_coordinates, c('ST_x','ST_y'), sep = 'x', convert = T) %>%
		separate(sn_coordinates, c('sn_x','sn_y'), sep = 'x', convert = T)
}