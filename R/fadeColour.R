# Function to fade a colour by specified alpha
	
fade.colour <- function(colours,alpha=1){
	rgb.dat <- col2rgb(colours)
	return(rgb(rgb.dat[1,],rgb.dat[2,],rgb.dat[3,],alpha*255,maxColorValue=255))
	}	