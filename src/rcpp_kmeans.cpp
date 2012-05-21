#include "rcpp_kmeans.h"
#include <sys/time.h>
#include <tr1/functional>
using namespace Rcpp ;


// TODO add more methods to calculate distances

const double eps = 1e-6;
typedef struct {
    int center_index; // one of dusts
    std::list<int> dusts; // we use remove
} cloud;
// typedef std::map<std::string, int> Point;
typedef std::map<size_t, int> Point;
typedef std::map<int, Point> Map;
typedef struct {
    std::map<int, double> product;
    std::vector<int> zeroes;
} products;
// we use double to avoiding memory overflow, but this is not a good idea anyway
// static int square[ (int) 1e+8];
// <c ,<i ,<j, r> > >
static std::map<int, std::map<int, products> > table;

// _EXPERIMENTAL_
static std::map<std::pair<int,int>, double> fulltable;
static double (*dist) (Point* xv, Point* yv);
//NOTE keys in std::map are sorted
static bool operator== (Point& l, Point& r)
{
    if (l.size() != r.size()) {
        return false;
    }
    Point::iterator l_it;
    Point::iterator r_it;
    for (l_it = l.begin(), r_it = r.begin(); l_it != l.end(); ++l_it, ++r_it) {
        if (*l_it != *r_it) {
            return false;
        }
    }
    return true;
}

//@{
// Euclidean distance
double euclidean (Point* xv, Point* yv)
{
    volatile double sum = 0.0;
    Point::iterator xv_it;
    Point::iterator yv_it;
    // match xs against ys
    for (xv_it = xv->begin(); xv_it != xv->end(); ++xv_it) {
        yv_it = yv->find (xv_it->first);
        if (yv_it != yv->end()) {
            sum += (xv_it->second - yv_it->second)
                   * (xv_it->second - yv_it->second);
        } else {
            sum += xv_it->second * xv_it->second;
        }
    }
    // match ys against xs, only take non-matched account
    for (yv_it = yv->begin(); yv_it != yv->end(); ++yv_it) {
        if (xv->find (yv_it->first) == xv->end()) {
            sum += yv_it->second * yv_it->second;
        }
    }
    return sum;
}
//@}

//@{
// Cosine distance
double cosine (Point* xv, Point* yv)
{
    if (*xv == *yv) return 0.0;
    volatile double sum = 0.0, xs = 0.0, ys = 0.0;
    Point::iterator xv_it;
    Point::iterator yv_it;
    //  TODO parallelize by using openmp
    // inner product, <x,y>, <x,x>
    for (xv_it = xv->begin(); xv_it != xv->end(); ++xv_it) {
        yv_it = yv->find (xv_it->first);
        if (yv_it != yv->end()) {
            sum += (xv_it->second * yv_it->second);
        }
        xs += xv_it->second * xv_it->second;
    }
    // inner product, <y,y>
    for (yv_it = yv->begin(); yv_it != yv->end(); ++yv_it) {
        ys += yv_it->second * yv_it->second;
    }
    return 1.0 - (sum / xs) * (sum / ys);
}
//@}

//@{
double distance (Map* freqTable_in, int x, int y)
{
    if (x == y) return 0.0;
    int x_, y_;
    double z;
    if (x < y) {
        x_ = x, y_ = y;
    } else {
        x_ = y, y_ = x;
    }
    // TODO query table first
    // _IF_ _NOT_ FOUND FROM table

    //      _THEN_
    //@{
    z = dist (& (*freqTable_in) [x_], & (*freqTable_in) [y_]);
    //@}
    return z;
}

//@}

//@{

/**
 * @brief ...
 *
 * @param freqTable_in ...
 * @param index point index
 * @param clouds ...
 * @param cloudSize ...
 * @return int the center of that cloud
 **/
int whichClosest (Map* freqTable_in, int index,
                  std::map<int, cloud> &clouds, int cloudSize)
{
    // avoid empty cloud
    // TODO keep or move one from another cloud that has a long distance from
    // its center into this (*clouds)[index]
    double dist_t = 0.0, dist = 1.0e+63;
    int which = index;
    int p = 0, q = 0; //p<<q
    // Performance bottle kneck
    // TODO try to make use of openmp to parallelize
    std::map<int, products>::iterator it;
    std::map<int, double>::iterator itt;
    for (int i = 0 ; i < cloudSize; ++i) {
        if (index < clouds [i].center_index) {
            p = index;
            q = clouds [i].center_index;
        } else {
            q = index;
            p = clouds [i].center_index;
        }
        it = table[i].find (p);
        if (it != table[i].end()) {
            itt = it->second.product.find (q);
            if (itt != it->second.product.end()) {
                dist_t = it->second.product[q];
            }
            int zeroes = it->second.zeroes.size();
            for (int j = 0; j < zeroes; ++j) {
                if (it->second.zeroes[j] == q) {
                    dist_t = 0.0;
                    break;
                }
            }
        }
        dist_t = distance (freqTable_in, p, q);
        if (dist_t < dist) {
            dist = dist_t;
            which = i;
        }
    }
    return which;
}
//@}


//@{
/**
 * @brief remove point @a value from center @a key
 *
 * @param clouds clusters
 * @param key the center index
 * @param value the dust that was included
 * @return void
 **/
inline bool removeFromCenters (std::map<int, cloud> & clouds,
                               const int key, const int value)
{
    // when initial
    if (key == -1) return true;
    if (clouds[key].dusts.size() > 1) {
        clouds[key].dusts.remove (value);
        // remove value from table[key] to save memory?
        return true;
    }
    // FIXME we could take one from other cloud that is far from its center
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
inline void addToCenters (std::map<int, cloud> & clouds, Map* freqTable_in,
                          int key, const int value)
{
    // What about to puttable in external file
    // and read in on call
    clouds [key].dusts.push_back (value);
    std::vector<int> dusts (clouds [key].dusts.begin(),
                            clouds [key].dusts.end());
    int dustsSize = dusts.size();
    int p = 0, q = 0; //p << q
    // OMP
    {
//         #pragma omp parallel for shared(dustsSize,dusts,table)
        for (int j = 0; j < dustsSize; ++j) {
            // could opt p==q
            if (value > dusts[j]) {
                p = dusts[j];
                q = value;
            } else {
                p = value;
                q = dusts[j];
            }
            std::map<int, products>::iterator it = table[key].find (p);
            if (it == table[key].end()) {
                double z = distance (freqTable_in, p, q);
                products tkp;
                if (z > eps) {
                    std::map<int, double> e;
                    tkp.product.insert (std::pair<int, double> (q, z));
                } else {
                    tkp.zeroes.push_back (q);
                }
                table[key].insert (std::pair<int, products> (p, tkp));
            } else {
                // check if zeroes has p
                int zeroes = it->second.zeroes.size();
                bool inZeroes = false;
                for (int i = 0; i < zeroes; ++i) {
                    if (it->second.zeroes[i] == q) {
                        inZeroes = true;
                        break;
                    }
                }
                if (!inZeroes) {
                    // check if products has p
                    std::map<int, double>::iterator itt =
                        it->second.product.find (q);
                    if (itt == it->second.product.end()) {
                        double z = distance (freqTable_in, p, q);
                        if (z < eps) {
                            it->second.zeroes.push_back (q);
                        } else {
                            it->second.product.insert
                            (std::pair<int, double> (q, z));
                        }
                    }
                }
            }
        }
    }
}
//@}


//@{
// the median is one that gets closest to others
// \min_i\max_j S_{i,j}
// O(n^2)

int median (std::list<int> & vec, int cloudIndex_in)
{
    int x = 0/*, y = 0*/;
    double max = 0.0;
    double min = 1.0e+63;
    volatile double e = 0.0;
    std::map<int, products>  & indexer = table[cloudIndex_in];
    // deep copy
    std::vector<int> dusts (vec.begin(), vec.end());
    std::sort (dusts.begin(), dusts.end());
    std::vector<int>::iterator x_it;
    std::vector<int>::iterator y_it;
    std::map<int, double>::iterator z_it;
    for (x_it = dusts.begin(); x_it != dusts.end(); ++x_it) {
        for (y_it = x_it + 1; y_it != dusts.end(); ++y_it) {
            products& p = indexer[*x_it];
            z_it = p.product.find (*y_it);
            if (z_it != p.product.end()) {
                e = z_it->second;
                if (max < e) {
                    max = e;
                }
            }
        }
        if (min > max) {
            min = max;
            x = *x_it;
        }
    }
    std::cout << "Median " << x << std::endl;
    return x;
}
//@}

//@{
void updateCloudCenters (std::map<int, cloud> & clouds,
                         int cloudSize)
{
    // FIXME expect median executing time diff slightly for all dusts
    // NOTE parallel threads are limited by cloudSize, however, luckly,
    // cloudSize is above 10 usually.
    {
//         #pragma omp parallel for
        for (int i = 0; i < cloudSize; ++i) {
            // What about to put table in external file
            // and read in on call
            clouds [i].center_index = median (clouds [i].dusts, i);
        }
    }
}
//@}
/**
 * @brief ...
 *
 * @param points_in ...
 * @param clusterSize_in ...
 * @param iter_in ...
 * @param method_in ...
 * @return SEXP
 **/
SEXP Kmeans (SEXP points_in, SEXP clusterSize_in, SEXP iter_in, SEXP method_in)
{
    // inputs should be validated in R
    BEGIN_RCPP
    List points (points_in);
    unsigned int pointsSize = points.size();
    switch (as<int> (method_in)) {
    case 1:
        dist = cosine;
        break;
    default:
        dist = euclidean;
    }
    int cloudSize = as<int> (clusterSize_in);
    //@{

    // Require CXX0X
//     std::hash< std::string > hashFun;
    std::tr1::hash<std::string> hashFun;
    Map freqTable; // < index <<hash,count>,<hash,count>,.. > >
    std::vector<std::string>::iterator v_it;
    volatile std::size_t v_t = 0LU;
    timeval tv0;
    timeval tv1;
    gettimeofday (&tv0, NULL);
    {
//         #pragma omp parallel for
        for (int i = 0; i < pointsSize; ++i) {
            Point v;
            std::vector<std::string> s = as<std::vector<std::string> > (points[i]);
            for (v_it = s.begin(); v_it != s.end(); ++v_it) {
                v_t = hashFun (*v_it);
                if (v.find ( (int) v_t) != v.end()) {
                    v[ (int) v_t] ++;
                } else {
                    v[ (int) v_t] = 1;
                }
            }
            freqTable[i] = v;
        }
    }
    gettimeofday (&tv1, NULL);
    std::cout << "TIME: "
              << (tv1.tv_sec - tv0.tv_sec) * 1000000 + (tv1.tv_usec - tv0.tv_usec)
              << std::endl;
    //@}
    //@{
    // initialize a cluster of size cluster
    // for std::random_shuffle
    srand ( (int) time (NULL));
    rand();
    // for quickly find out which cluster a point belongs to
    // cluster to each points, -1 indicates not assigned yet
    // the value of @a cluster is the center of that cluster
    std::vector<int> cluster (pointsSize), tmp;
    std::fill (cluster.begin(), cluster.end(), -1);

    // make sure that no same points as centers
    // check first cloudSize points

    std::vector<int> shuffle (pointsSize);
    for (int i = 0; i < pointsSize; ++i) {
        shuffle[i] = i;
    }
    std::random_shuffle (shuffle.begin(), shuffle.end());
    std::vector<int>::iterator sh_it = shuffle.begin();
    gettimeofday (&tv0, NULL);
    volatile bool skip = false;
    tmp.push_back (*sh_it);
    for (int i = 1; i < cloudSize; ++i) {
        // attempt to get a distinct point
        do {
            skip = false;
            sh_it++;
            if (sh_it == shuffle.end()) {
                std::cout << "Error on reach end" << std::endl;
                break;
            }
            for (int j = 0; j < i; ++j) {
                // duplicated found
                if (freqTable[*sh_it] == freqTable[tmp[j]]) {
                    /* std::cout << ".. " << *sh_it << " .."<< std::endl; */
                    skip = true;
                    break;
                }
            }
        } while (skip);
        if (sh_it != shuffle.end())
            tmp.push_back (*sh_it);
    }
    gettimeofday (&tv1, NULL);
    std::cout << "TIME: "
              << (tv1.tv_sec - tv0.tv_sec) * 1000000 + (tv1.tv_usec - tv0.tv_usec)
              << std::endl;
    shuffle.clear();
    gettimeofday (&tv0, NULL);

    // key is a placeholder
    std::map<int, cloud> clouds;
    for (int i = 0; i < cloudSize; ++i) {
        std::list<int> l;
        l.push_back (tmp[i]);
        clouds[i].dusts = l;
        clouds[i].center_index = tmp[i];
        cluster[tmp[i]] = i;
//         std::cout << tmp[i] << std::endl;
    }
    gettimeofday (&tv1, NULL);
    std::cout << "TIME: "
              << (tv1.tv_sec - tv0.tv_sec) * 1000000 + (tv1.tv_usec - tv0.tv_usec)
              << std::endl;
    tmp.clear();
    //@}

    volatile int which = 0;
    int iter = as<int> (iter_in); // assert positive in R level
    volatile int changed = 0;
    // TO sort clouds by center_index and adjust cluster in each iteration
    while (iter--) {
        changed = 0;
        std::cout << "Loop Begin" << std::endl;
        gettimeofday (&tv0, NULL);
        for (int i = 0; i < pointsSize; ++i) {
            // skip centers as they ain't changed
            if (cluster[i] == -1 ||
                    clouds [cluster[i]].center_index != cluster[i]) {
                which = whichClosest (&freqTable, i, clouds, cloudSize);
                if (cluster[i] != which) {
                    /* change cluster[i] state */
                    /* remove */
                    if (removeFromCenters (clouds, cluster[i], i)) {
                        /* update */
                        cluster[i] = which;
                        addToCenters (clouds, &freqTable, which, i);
                        // record which
                        changed ++;
//                         std::cout << "Looping ..." << std::endl;
                    }
                }
            }
        }
        gettimeofday (&tv1, NULL);
        std::cout << (tv1.tv_sec - tv0.tv_sec) * 1000000
                  + (tv1.tv_usec - tv0.tv_usec)
                  << std::endl << "Loop End " << std::endl;
        if (changed == 0) break;
        // TODO what about to let table update on removeFromCenters and
        // addToCenters
        // TODO it's better that only update which and cluster[i]
        // use a vector to store (which and cluster[i])
        std::cout << "updateCloudCenters Begin" << std::endl;
        gettimeofday (&tv0, NULL);
        //@{
        updateCloudCenters (clouds, cloudSize);
        //@}
        gettimeofday (&tv1, NULL);
        std::cout << (tv1.tv_sec - tv0.tv_sec) * 1000000
                  + (tv1.tv_usec - tv0.tv_usec)
                  << std::endl << "updateCloudCenters End " << std::endl;
    }
    std::vector<IntegerVector> x;
    for (int i = 0; i < cloudSize; ++i) {
        x.push_back (IntegerVector (clouds[i].dusts.begin(),
                                    clouds[i].dusts.end()));
    }
    List z =  List::create (_["clusters"] = x,
                            _["iterations"] = as<int> (iter_in) - iter);
    // clean up
    clouds.clear();
    freqTable.clear();
    table.clear();
    return z;
    END_RCPP

}

// handle small samples and return centers
// merge all centers above and return final centers
// loop all to find out which centers they get nearest



// kate: indent-mode cstyle; indent-width 4; replace-tabs on;


