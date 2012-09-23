## points clusterSize iter method
Kmeans <- function(strings, cluster.no=10, max.iter=100, epsilon=0.001){
    if (missing(strings)) stop("Kmeans(points, cluster.no,iter, epsilon)")
    if (cluster.no < 2) stop("cluster.no must be greater than 1.")
    if (cluster.no < 2) stop("max.iter must be greater than 1.")
    ## ##
    ## TODO: 1 is cosine, no other has been implemented
    ## ##
	.Call( "Kmeans", strings, cluster.no, max.iter, 1, epsilon, PACKAGE = "RcppKmeans" )
}
