## TODO validate inputs
## points clusterSize iter method
Kmeans <- function(x,y,z,w){
	.Call( "Kmeans", x, y, z, w, PACKAGE = "RcppKmeans" )
}
