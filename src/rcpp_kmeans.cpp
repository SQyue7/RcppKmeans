#include "rcpp_kmeans.h"
#include <sys/time.h>
#include <tr1/functional>
using namespace Rcpp ;
// TODO add more methods to calculate distances
static double eps = 0.0001;
typedef struct {
    int center_index; // one of dusts
    std::list<int> dusts; // we use remove
} cloud;

typedef std::pair<int, int> iipair;
typedef std::map<iipair, double> fulltableType;
// typedef std::map<std::string, int> Point;
typedef std::map<size_t, int> Point;
typedef std::map<int, Point> Map;
typedef struct {
    std::map<int, double> product;
    std::vector<int> zeroes;
} products;
static double (*dist) (Point* xv, Point* yv);

//NOTE keys in std::map are sorted
static bool operator== (Point& l, Point& r)
{
    if (l.size() != r.size()) {
        return false;
    }
    Point::iterator l_it = l.begin();
    Point::iterator r_it = r.begin();
    for (; l_it != l.end() && r_it != r.end(); ++l_it, ++r_it) {
        if (*l_it != *r_it) {
            return false;
        }
    }
    return true;
}

//@{
//FIXME
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
    double z = 0.0;
    if (x < y) {
        z = dist (& (*freqTable_in) [x], & (*freqTable_in) [y]);
    } else {
        z = dist (& (*freqTable_in) [y], & (*freqTable_in) [x]);
    }
    return z;
}
//@}

// NOTE for cosine, 1 - eps as MAX
// NOTE Don't use it, it require a lot of RAM
// void generateFullTable (Map* freqTable_in)
// {
//     int dim = freqTable_in->size();
//     volatile double r = 0.0;
//     {
//         #pragma omp parallel for
//         for (int i = 0; i < dim; ++i) {
//             for (int j = i + 1; j < dim; ++j) {
//                 r = distance (freqTable_in, i, j);
//                 /* std::cout << r << " "; */
//                 if (r < 1 - eps)
//                     fulltable.insert (std::pair<iipair, double>
//                                       (iipair (i, j), (double) r));
//             }
// //             std::cout << i << std::endl;
//         }
//     }
// }

void addToFullTable (Map* freqTable_in, fulltableType& table, const int x_in, const int y_in)
{
    double r = 0.0;
    if (x_in == y_in) return;
    if (x_in < y_in) {
        #pragma omp flush(table)
        if (table.find (iipair (x_in, y_in)) == table.end()) {
            r = distance (freqTable_in, x_in, y_in);
            if (r < 1 - eps) {
                table.insert (std::pair<iipair, double>
                              (iipair (x_in, y_in), (double) r));
            }
        }
    } else {
        #pragma omp flush(table)
        if (table.find (iipair (y_in, x_in)) == table.end()) {
            r = distance (freqTable_in, y_in, x_in);
            if (r < 1 - eps) {
                table.insert (std::pair<iipair, double>
                              (iipair (y_in, x_in), (double) r));
            }
        }
    }
}

// value_in against list_in
void addToFullTable (Map* freqTable_in, fulltableType& table, std::list<int> & list_in, const int value_in)
{
    int dim = list_in.size();
    {
        #pragma omp parallel for
        for (int i = 0; i < dim; ++i) {
            addToFullTable (freqTable_in, table, i, value_in);
        }
    }
}

void getDistance (fulltableType& table, int x, int y, double& r)
{
    std::map<iipair, double>::iterator it;
    if (x > y) {
        it = table.find (iipair (y, x));
    } else {
        it = table.find (iipair (x, y));
    }
    if (it != table.end()) {
        r = it->second;
    } else {
        r = 1.0;
    }
}

int whichClosest (Map* freqTable_in, fulltableType& table, int index,
                  std::map<int, cloud> &clouds, int cloudSize)
{
    double dist_t = 0.0, dist = 1.0e+63;
    int which = 0;
    // TODO try to make use of openmp to parallelize
    std::map<iipair, double>::iterator  it;
    for (int i = 0 ; i < cloudSize; ++i) {
        addToFullTable (freqTable_in, table, index, clouds [i].center_index);
        getDistance (table, index, clouds [i].center_index, dist_t);
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
                          fulltableType& table, int key_in, const int value_in)
{
    addToFullTable (freqTable_in, table, clouds [key_in].dusts, value_in);
    clouds [key_in].dusts.push_back (value_in);
}
//@}

int median (std::list<int> & vec, fulltableType& table)
{
    int x = 0/*, y = 0*/;
    double max = 0.0;
    double min = 1.0e+63;
    double e = 1.0;
    vec.sort();
    std::list<int>::iterator vecItx, vecIty;
    for (vecItx = vec.begin(); vecItx != vec.end(); ++vecItx) {
        for (vecIty = vec.begin(); vecIty != vec.end(); ++vecIty) {
            getDistance (table, *vecItx, *vecIty, e);
            if (max < e) {
                max = e;
            }
        }
        if (min > max) {
            min = max;
            x = *vecItx;
        }
    }
//     std::cout << "MedianX " << x << std::endl;
    return x;
}


//@{
void updateCloudCenters (std::map<int, cloud> & clouds,
                         int cloudSize, fulltableType& table)
{
    {
        #pragma omp parallel for
        for (int i = 0; i < cloudSize; ++i) {
            /** clouds [i].center_index = median (clouds [i].dusts, i); */
            int index = median (clouds [i].dusts, table);
            #pragma omp critical
            clouds [i].center_index = index;
        }
    }
}
//@}

//@{
/**
 * @brief ...
 *
 * @param points_in ...
 * @param clusterSize_in ...
 * @param iter_in ...
 * @param method_in ...
 * @return SEXP
 **/
SEXP Kmeans (SEXP points_in, SEXP clusterSize_in,
             SEXP iter_in, SEXP method_in, SEXP epsilon_in)
{
    // inputs should be validated in R
    BEGIN_RCPP
    double epsilon = as<double> (epsilon_in);
    if (epsilon > 0 && epsilon < 1.0) {
        eps = epsilon;
    }
    List points (points_in);
    unsigned int pointsSize = points.size();
//     switch (as<int> (method_in)) {
//     case 1:
    dist = cosine;
//         break;
//     default:
//         dist = euclidean;
//     }
    fulltableType fulltable;
    int cloudSize = as<int> (clusterSize_in);
    std::tr1::hash<std::string> hashFun;
    /**
     * it's better not to use hash but exactly match
     */
    Map freqTable; // < index <<hash,count>,<hash,count>,.. > >
    std::vector<std::string>::iterator v_it;
    volatile std::size_t v_t = 0LU;
    timeval tv0;
    timeval tv1;
    //@{
    gettimeofday (&tv0, NULL);
    //@}
    /**
     * TODO parallelization requires careful handling
     * #pragma omp parallel for
     */
    {
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
    //@{
    gettimeofday (&tv1, NULL);
    std::cout << "freqTable: "
              << (tv1.tv_sec - tv0.tv_sec) * 1000000 + (tv1.tv_usec - tv0.tv_usec)
              << std::endl;
    //@}
    // initialize a cluster of size cluster
    // for std::random_shuffle
    srand ( (int) time (NULL));
    rand();
    // the value of cluster's element is the index of cloud that it belongs to
    std::vector<int> cluster (pointsSize), tmp;
    std::fill (cluster.begin(), cluster.end(), -1);
    std::vector<int> shuffle (pointsSize);
    for (int i = 0; i < pointsSize; ++i) {
        shuffle[i] = i;
    }
    std::random_shuffle (shuffle.begin(), shuffle.end());
    std::vector<int>::iterator shuffleIt = shuffle.begin();
    //@{
    gettimeofday (&tv0, NULL);
    //@}
    volatile bool skip = false;
    tmp.push_back (*shuffleIt);
    for (int i = 1; i < cloudSize; ++i) {
        // attempt to get a distinct point
        do {
            skip = false;
            shuffleIt++;
            if (shuffleIt == shuffle.end()) {
                std::cout << "Error on reach end" << std::endl;
                std::cout << "epsilon is too large" << std::endl;
                freqTable.clear();
                fulltable.clear();
                return 0;
            }
            for (int j = 0; j < i; ++j) {
                double d = 0.0;
                addToFullTable (&freqTable, fulltable, *shuffleIt, tmp[j]);
                getDistance (fulltable, *shuffleIt, tmp[j], d);
                if (d < eps / 10) {
                    skip = true;
                    break;
                }
            }
        } while (skip);
        if (shuffleIt != shuffle.end())
            tmp.push_back (*shuffleIt);
    }
    //@{
    gettimeofday (&tv1, NULL);
    std::cout << "TIME: "
              << (tv1.tv_sec - tv0.tv_sec) * 1000000 + (tv1.tv_usec - tv0.tv_usec)
              << std::endl;
    //@}
    shuffle.clear();
    //@{
    gettimeofday (&tv0, NULL);
    //@}
    // key is a placeholder
    std::map<int, cloud> clouds;
    for (int i = 0; i < cloudSize; ++i) {
        std::list<int> l;
        l.push_back (tmp[i]);
        clouds[i].dusts = l;
        clouds[i].center_index = tmp[i];
        cluster[tmp[i]] = i;
//         std::cout << "#" << tmp[i] << " ";
    }
    std::cout << std::endl;
    //@{
    gettimeofday (&tv1, NULL);
    std::cout << "TIME: "
              << (tv1.tv_sec - tv0.tv_sec) * 1000000
              + (tv1.tv_usec - tv0.tv_usec)
              << std::endl;
    //@}
    tmp.clear();
    volatile int which = 0;
    int iter = as<int> (iter_in); // assert positive in R level
    volatile int changed = 0;
    // BUG if points must be clustered no more than k groups
    // then a duplicated center would hold the lace as increase k
    // fortunately points is large and k is sufficient to avoid this case
    while (iter--) {
        changed = 0;
        //@{
        std::cout << "Loop Begin" << std::endl;
        gettimeofday (&tv0, NULL);
        //@}
        for (int i = 0; i < pointsSize; ++i) {
            if (cluster[i] == -1 /* when initial */ ||
                    clouds [cluster[i]].center_index != i) {
                which = whichClosest (&freqTable, fulltable, i, clouds, cloudSize);
                if (cluster[i] != which) {
                    /* change cluster[i] state */
                    /* remove */
                    if (removeFromCenters (clouds, cluster[i], i)) {
                        /* update */
                        cluster[i] = which;
                        addToCenters (clouds, &freqTable, fulltable, which, i);
                        // record which
                        changed ++;
                    }
                }
            }
        }
        //@{
        gettimeofday (&tv1, NULL);
        std::cout << (tv1.tv_sec - tv0.tv_sec) * 1000000
                  + (tv1.tv_usec - tv0.tv_usec)
                  << std::endl << "Loop End " << std::endl;
        //@}
        if (changed == 0) break;
//         for (int i = 0; i < cloudSize; ++i) {
//             std::cout << "# " << clouds[i].center_index << " ";
//         }
//         std::cout << std::endl;
        //@{
        std::cout << "updateCloudCenters Begin" << std::endl;
        gettimeofday (&tv0, NULL);
        //@}
        updateCloudCenters (clouds, cloudSize, fulltable);
        //@{
        gettimeofday (&tv1, NULL);
        std::cout << (tv1.tv_sec - tv0.tv_sec) * 1000000
                  + (tv1.tv_usec - tv0.tv_usec)
                  << std::endl << "updateCloudCenters End " << std::endl;
        //@}
//         for (int i = 0; i < cloudSize; ++i) {
//             std::cout << "# " << clouds[i].center_index << " ";
//         }
//         std::cout << std::endl;
    }
    std::vector<IntegerVector> x;
    for (int i = 0; i < cloudSize; ++i) {
        x.push_back (IntegerVector (clouds[i].dusts.begin(),
                                    clouds[i].dusts.end()));
    }
    List z =  List::create (_["clusters"] = x,
                            _["iterations"] = as<int> (iter_in) - iter);
    clouds.clear();
    freqTable.clear();
    fulltable.clear();
    return z;
    END_RCPP

}
//@}
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
