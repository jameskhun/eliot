# Function to subtract mean and divide by standard deviation

standardize <- function(x) return((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))

# Apply standardize within subgroups defined by "by"

standardize.by <- function(x,by){
	by <- as.factor(by)
	by.levels <- levels(by)
	out <- rep(NA,length(x))
	for (i in 1:length(by.levels)){
		out[as.numeric(by) == i] <- standardize(x[as.numeric(by) == i])
	}
	return(out)	
	}
