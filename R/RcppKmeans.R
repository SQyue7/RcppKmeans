## TODO validate inputs
## points clusterSize iter method
Kmeans <- function(x,y,z,e){
    if(missing(e))
        e <- 0;
    if (missing(z)) z <- 100;
    if (missing(y)) stop("Kmeans(points, cluster.no,iter, epsilon)")
    ## ##
    ## TODO: 1 is cosine, no other has been implemented
    ## ##
	.Call( "Kmeans", x, y, z, 1, e, PACKAGE = "RcppKmeans" )
}
