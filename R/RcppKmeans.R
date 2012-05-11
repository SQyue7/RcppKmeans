## TODO validate inputs
Kmeans <- function(x,y,z,w){
	.Call( "Kmeans", x, y, z, w, PACKAGE = "RcppKmeans" )
}

