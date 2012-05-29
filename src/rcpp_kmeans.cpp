#include "rcpp_kmeans.h"
#ifdef VERBOSE
#include <sys/time.h>
#endif
#include <omp.h>

//@{
// Tip
// Avoid writing to global data that is accessed from multiple threads.
// Tip
// Align shared global data to cache line boundaries.
// Tip
// Don't store temporary, thread specific data in an array indexed by the thread id or rank.
// Tip
// When parallelizing an algorithm, partition data sets along cache lines, not across cache lines.
//@}

using namespace Rcpp ;
/**
 * Types, Constants and functions
 */
static double eps = 0.0005;
typedef struct {
    int center_index; // one of dusts
    std::list<int> dusts; // we use remove
} cloud;
typedef std::map<std::string, int> Point;
typedef std::map<int, Point> Map;
static double (*dist) (Point* xv, Point* yv);
/**
 * Distances
 */
// TODO add more methods to calculate distances
//@{
//FIXME NOT CORRECT
// Euclidean distance
/**
 * @brief Euclidian distance
 *
 * @param xv ...
 * @param yv ...
 * @return double
 **/
double euclidean (Point* xv_in, Point* yv_in)
{
    volatile double sum = 0.0;
    Point::iterator xv_it;
    Point::iterator yv_it;
    // match xs against ys
    for (xv_it = xv_in->begin(); xv_it != xv_in->end(); ++xv_it) {
        yv_it = yv_in->find (xv_it->first);
        if (yv_it != yv_in->end()) {
            sum += (xv_it->second - yv_it->second)
                   * (xv_it->second - yv_it->second);
        } else {
            sum += xv_it->second * xv_it->second;
        }
    }
    for (yv_it = yv_in->begin(); yv_it != yv_in->end(); ++yv_it) {
        if (xv_in->find (yv_it->first) == xv_in->end()) {
            sum += yv_it->second * yv_it->second;
        }
    }
    return sum;
}
//@}

inline double norm (Point* xv_in)
{
    double n = 0.0;
    Point::iterator xv_it = xv_in->begin();
    for (xv_it = xv_in->begin(); xv_it != xv_in->end(); ++xv_it) {
        n += xv_it->second * xv_it->second;
    }
    return sqrt (n);
}

//@{
// Cosine distance
/**
 * @brief Cosine dissimilarity (1 - cos<x,y>)
 *
 * @param xv ...
 * @param yv ...
 * @return double
 **/
double cosine (Point* xv_in, Point* yv_in)
{
    if (*xv_in == *yv_in) return 0.0;
    if (xv_in->size() == 0  || yv_in->size() == 0) return 0.0;
    volatile double sum = 0.0, xs = 0.0, ys = 0.0;
    Point::iterator xv_it = xv_in->begin();
    Point::iterator zv_it = yv_in->begin();
    Point::iterator yv_it;
    // inner product, <x,y>, <x,x>
    for (; xv_it != xv_in->end(); ++xv_it) {
        for (yv_it = zv_it; yv_it != yv_in->end(); ++yv_it) {
            if (yv_it->first == xv_it->first) {
                sum += yv_it->second * xv_it->second;
                zv_it = yv_it;
                break;
            }
        }
    }
    xs = norm (xv_in);
    ys = norm (yv_in);
    return 1.0 - sum / (xs * ys);
}
//@}
// <x,y> = cos<x,y> |x|*|y|
//1-cos<x,y>= 1 - <x,y> /(|x|*|y|)
// <x,y> /|x|
double projection (Point* xv_in, Point* yv_in)
{
    if (*xv_in == *yv_in) return 1.0;
    if (xv_in->size() == 0  || yv_in->size() == 0) return 0.0;
    volatile double sum = 0.0;
    Point::iterator xv_it = xv_in->begin();
    Point::iterator zv_it = yv_in->begin();
    Point::iterator yv_it;
    // inner product, <x,y>
    for (; xv_it != xv_in->end(); ++xv_it) {
        for (yv_it = zv_it; yv_it != yv_in->end(); ++yv_it) {
            if (yv_it->first == xv_it->first) {
                sum += yv_it->second * xv_it->second;
                zv_it = yv_it;
                break;
            }
        }
    }
    return sum / norm (xv_in);
}
//@{
/**
 * @brief Interface to dist for caller
 *
 * @param freqTable_in ...
 * @param x ...
 * @param y ...
 * @return double
 **/
double distance (Map* freqTable_in, int x_in, int y_in)
{
    if (x_in == y_in) return 0.0;
    return dist (& (*freqTable_in) [x_in], & (*freqTable_in) [y_in]);
}
//@}



/**
 * @brief find one center_index that is nearest to @a index_in
 *
 * @param freqTable_in ...
 * @param index_in ...
 * @param centerIndex_in ...
 * @param centers_in ...
 * @return int
 **/
int whichClosest (Map* freqTable_in,
                  int index_in,
                  int centerIndex_in, // pre-center = centers_in(centerIndex_in)
                  std::map<int, cloud>  &clouds,
                  int centersCount_in)
{

    double dist_t = 0.0;
    double dist = eps * (2 - eps);// eps + eps(1-eps)
    int which = centerIndex_in;
    for (int i = 0 ; i < centersCount_in; ++i) {
        dist_t = distance (freqTable_in, index_in, (clouds [i]).center_index);
        if (dist_t < dist) {
            dist = dist_t;
            which = i;
        }
    }
    return which;
}


//@{
/**
 * @brief remove point @a value from center @a key
 *
 * @param clouds clusters
 * @param key the center index
 * @param value the dust that was included
 * @return void
 **/
inline bool removeFromCenters (std::map<int, cloud> & clouds_io,
                               const int key_in, const int value_in)
{
    if (key_in == -1) return true;
    if (clouds_io.find (key_in) != clouds_io.end() &&
            clouds_io[key_in].dusts.size() > 1) {
        clouds_io[key_in].dusts.remove (value_in);
        return true;
    }
    return false;
}
//@}

//@{
/**
 * @brief add @a value to the list centers[keys]
 *
 * @param centers Pointer to centers table
 * @param key
 * @param value ...
 * @return int
 **/
inline void addToCenters (std::map<int, cloud> & clouds,
                          int key_in, const int value_in)
{
    clouds [key_in].dusts.push_back (value_in);
}
//@}

// NOTE performance battle-neck
int centroidX (std::list<int> & vec_in, Map* freqTable_in)
{
    int s = vec_in.size();
    if (s == 0) {
        std::cout << "!!_BUG_!!" << std::endl;
        return 0;
    }
#ifdef VERBOSE
    //@{
    timeval tv0;
    timeval tv1;
    gettimeofday (&tv0, NULL);
    std::cout << "centroidX "  <<  omp_get_thread_num()
              << " begins" << std::endl;
    std::cout << "Size " << s << std::endl;
    //@}
#endif
    Point p;
    Point::iterator p_it, ep_it;
    std::list<int>::iterator vec_it;
    for (vec_it = vec_in.begin(); vec_it != vec_in.end(); ++vec_it) {
        Point& ep = (*freqTable_in) [*vec_it];
        for (ep_it = ep.begin(); ep_it != ep.end(); ++ep_it) {
            p_it = p.find (ep_it->first);
            if (p_it != p.end()) {
                p_it->second += ep_it->second;
            } else {
                p.insert (*ep_it);
            }
        }
    }
    std::vector<int> vec (vec_in.begin(), vec_in.end());
    std::vector<double> es (s);
    std::fill (es.begin(), es.end(), 0.0);
    // get a point nearest to centroid p
    // the larger the projection is, the closer the point gets to p
    {
        #pragma omp parallel for shared(p,vec,es)
        for (int i = 0; i < s; ++i) {
            #pragma omp flush(es)
            es[i] = projection (& (*freqTable_in) [vec[i]], &p);
            #pragma omp flush(es)
        }
    }

    volatile double large = 0.0;
    volatile int x = 0;
    for (int i = 0; i < s; ++i) {
        if (es[i] > large) {
            large = es[i];
            x = vec[i];
        }
    }
#ifdef VERBOSE
//@{
    gettimeofday (&tv1, NULL);
    std::cout << (tv1.tv_sec - tv0.tv_sec) * 1000000
              + (tv1.tv_usec - tv0.tv_usec)
              << std::endl << "centroidX " <<  omp_get_thread_num()
              << " ends" << std::endl;
//@}
#endif
    return x;
}


void updateCloudCentersX (std::map<int, cloud> & clouds,
                          int cloudSize, Map* freqTable_in)
{
    // NOTE clouds[cloudSize] is home for all homeless, replace center with next one
    for (int i = 0; i < cloudSize; ++i) {
        clouds [i].center_index = centroidX (clouds [i].dusts, freqTable_in);
    }
}
//@{
/**
 * @brief ...
 *
 * @param points_in ...
 * @param clusterSize_in ...
 * @param iter_in ...
 * @return SEXP
 **/
SEXP Kmeans (SEXP points_in, SEXP clusterSize_in,
             SEXP iter_in, SEXP method_in, SEXP epsilon_in)
{
    // input variables should be validated in R
//     BEGIN_RCPP
    double epsilon = as<double> (epsilon_in);
    if (epsilon < eps) {
        std::cout << "Given epsilon " << epsilon << " is less than " << eps
                  << std::endl << "Using default epsilon " << eps << std::endl;
    } else if (epsilon > 1.0) {
        std::cout << "Given epsilon " << epsilon << " must be less than 1.0"
                  << std::endl << "Using default epsilon " << eps << std::endl;
    } else {
        eps = epsilon;
    }
    List points (points_in);
    unsigned int pointsSize = points.size();
//     switch (as<int> (method_in)) {
    dist = cosine;

    int cloudSize = as<int> (clusterSize_in);
    /**
     * it's better to match exactly rather than hash
     */
    Map freqTable;
    std::vector<std::string>::iterator v_it;
#ifdef VERBOSE
    //@{
    timeval tv0;
    timeval tv1;
    gettimeofday (&tv0, NULL);
    //@}
#endif
    {
        for (int i = 0; i < pointsSize; ++i) {
            Point v;
            std::vector<std::string> s =
                as<std::vector<std::string> > (points[i]);
            for (v_it = s.begin(); v_it != s.end(); ++v_it) {
                if (v.find (*v_it) != v.end()) {
                    v[*v_it] ++;
                } else {
                    v[*v_it] = 1;
                }
            }
            freqTable[i] = v;
        }
    }
#ifdef VERBOSE
    //@{
    gettimeofday (&tv1, NULL);
    std::cout << "freqTable: "
              << (tv1.tv_sec - tv0.tv_sec) * 1000000
              + (tv1.tv_usec - tv0.tv_usec)
              << std::endl;
    //@}
#endif
    // initialize a cluster of size cluster
    // for std::random_shuffle
    srand ( (int) time (NULL));
    rand();
    // the value of cluster's element is the index of cloud that it belongs to
    std::vector<int> cluster (pointsSize);
    std::vector<int> clusterTmp (pointsSize);
    for (int i = 0; i < pointsSize; ++i) {
        cluster[i] = i;
    }
    std::random_shuffle (cluster.begin(), cluster.end());
    std::vector<int>::iterator clusterIt = cluster.begin();
#ifdef VERBOSE
    //@{
    gettimeofday (&tv0, NULL);
    //@}
#endif
    volatile bool skip = false;
    //TODO free
    int* tmp = (int*) alloca (cloudSize * sizeof (int));
    tmp[0]  = *clusterIt;
    for (int i = 1; i < cloudSize; ++i) {
        // attempt to get a distinct point
        do {
            skip = false;
            clusterIt++;
            if (clusterIt == cluster.end()) {
                std::cout << "Error on reach end when initialize center["
                          << i << "]" << std::endl;
                std::cout << "epsilon (" << eps << ") is too large" << std::endl;
                return R_NilValue;
            }
            for (int j = 0; j < i; ++j) {
                double d = distance (&freqTable, *clusterIt, tmp[j]);
                // if found one close enough
                if (d < eps) {
                    skip = true;
                    break;
                }
            }
        } while (skip);
        tmp[i] = *clusterIt;
    }
    std::cout << "skip " << (skip ? "true" : "flase") << std::endl;
#ifdef VERBOSE
//@{
    gettimeofday (&tv1, NULL);
    std::cout << "tmp: "
              << (tv1.tv_sec - tv0.tv_sec) * 1000000
              + (tv1.tv_usec - tv0.tv_usec)
              << std::endl;
    //@}
#endif
    std::fill (cluster.begin(), cluster.end(), -1);
    std::fill (clusterTmp.begin(), clusterTmp.end(), -1);
#ifdef VERBOSE
    //@{
    gettimeofday (&tv0, NULL);
    //@}
#endif
    // TODO free them
    // centers of clounds
    int* centers = (int*) malloc (cloudSize * sizeof (int));
    // key is a placeholder
    std::map<int, cloud> clouds;
#ifdef VERBOSE
    std::cout << "T ";
#endif
    for (int i = 0; i < cloudSize; ++i) {
        std::list<int> l;
        l.push_back (tmp[i]);
        clouds[i].dusts = l;
        clouds[i].center_index = tmp[i];
        centers[i] = tmp[i];
        cluster[tmp[i]] = i;
#ifdef VERBOSE
        std::cout << centers[i] << " ";
#endif
    }
#ifdef VERBOSE
    std::cout << std::endl;
#endif
#ifdef VERBOSE
    //@{
    gettimeofday (&tv1, NULL);
    std::cout << "clouds: "
              << (tv1.tv_sec - tv0.tv_sec) * 1000000
              + (tv1.tv_usec - tv0.tv_usec)
              << std::endl;
    //@}
#endif
    int iter = as<int> (iter_in); // assert positive in R level
    bool changed = false;
    while (--iter) {
        changed = false;
#ifdef VERBOSE
        //@{
        std::cout << "Loop Begin" << std::endl;
        gettimeofday (&tv0, NULL);
        //@}
#endif
        //@{
        // 1. use a temporary array to store updated centers
        //    and parallelize this part
        //@}
        //@{
        {
            #pragma omp parallel for
            for (int i = 0; i < pointsSize; ++i) {
                #pragma omp flush(clusterTmp,cluster)
                clusterTmp[i] = whichClosest (&freqTable, i,
                                              cluster[i],
                                              clouds, cloudSize);
                #pragma omp flush(clusterTmp)
            }
            #pragma omp barrier
        }
        //@}
        //@{
        // 2. for each point, if it moves to new cloud, update it,
        //    mark changed true, if one does
        //@}
        for (int i = 0; i < pointsSize; ++i) {
            // At initial, clusterTmp[i] == -1 means no center is close enough
            // to @a i and as a result, @a i may be isolated eventually
            // while @a eps is relatively large.
            if (cluster[i] != clusterTmp[i]
                    || (clusterTmp[i] != cloudSize
                        && distance (&freqTable,
                                     centers[cluster[i]],
                                     centers[clusterTmp[i]]) > eps)) {
                /** change cluster[i] state
                *   remove and add to
                **/
                if (removeFromCenters (clouds, cluster[i], i)) {
                    /* update */
                    cluster[i] = clusterTmp[i];
                    addToCenters (clouds, clusterTmp[i], i);
                    // record clusterTmp[i]
                    changed = true;
                }
            }
        }
#ifdef VERBOSE
        //@{
        gettimeofday (&tv1, NULL);
        std::cout << (tv1.tv_sec - tv0.tv_sec) * 1000000
                  + (tv1.tv_usec - tv0.tv_usec)
                  << std::endl << "Loop End" << std::endl;
        //@}
#endif
        //@{
        // 3. break if no change else update each centroids
        // 4. got 1.
        //@}
        if (not changed) break;
#ifdef VERBOSE
        //@{
        std::cout << "updateCloudCenters Begin" << std::endl;
        gettimeofday (&tv0, NULL);
        //@}
#endif
        //@{
        // IMPROVE ME performance battle-neck
        updateCloudCentersX (clouds, cloudSize, &freqTable);
        //@}
#ifdef VERBOSE
        //@{
        gettimeofday (&tv1, NULL);
        std::cout << (tv1.tv_sec - tv0.tv_sec) * 1000000
                  + (tv1.tv_usec - tv0.tv_usec)
                  << std::endl << "updateCloudCenters End" << std::endl;
        //@}
#endif
        // compare centers and clouds.index
        {
            #pragma omp parallel for
            for (int i = 0; i < cloudSize; ++i) {
                #pragma omp flush(centers,clouds)
                centers[i] = clouds[i].center_index;
                #pragma omp flush(centers,clouds)
            }
            #pragma omp barrier
        }
#ifdef VERBOSE
        //@{
        // TODO output current clusters, it's nice for visualization
        // Perhaps we could have a std::vector<clouds>
        //@}
#endif
    }
    /**
     * KMEANS BUG
     * suppose there are k+1 groups of vectors, whose spaces are Orthogonal
     * each other, hence it's possible that a tiny cloud might be clusted while
     * a big one not.
     */
    int ccsum = 0;
    for (int i = 0; i < cloudSize; ++i) {
        ccsum += clouds[i].dusts.size();
    }
    std::vector<int> divergent;
    if (ccsum < pointsSize) {
        for (int i = 0; i < pointsSize; ++i) {
            if (cluster[i] == -1) divergent.push_back (i);
        }
    }
    std::vector<IntegerVector> x;
    for (int i = 0; i < cloudSize; ++i) {
        x.push_back (IntegerVector (clouds[i].dusts.begin(),
                                    clouds[i].dusts.end()));
    }
    List z =  List::create (_["clusters"] = x,
                            _["iterations"] = as<int> (iter_in) - iter,
                            _["divergent"] = divergent);
    return z;
//     END_RCPP
}
//@}
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;




