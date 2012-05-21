## TODO validate inputs
## points clusterSize iter method
Kmeans <- function(x,y,z,w,e){
    if(missing(e))
        e <- 0;
	.Call( "Kmeans", x, y, z, w, e, PACKAGE = "RcppKmeans" )
}
