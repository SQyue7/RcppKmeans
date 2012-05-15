#include "rcpp_kmeans.h"
#include <functional>
using namespace Rcpp ;

//ISSUE

// TODO add more methods to calculate distances

typedef struct {
    int center_index; // one of dusts
    std::list<int> dusts; // we use remove
} cloud;

// typedef std::map<std::string, int> Point;
typedef std::map<size_t, int> Point;
typedef std::map<int, Point> Map;

// NOTE the table may be as big as i*j*8 bytes
// hmm, requires lots of RAM!!

//  <i ,<j, r> >
// static std::map<int, std::map<int, int> > *table;

static long double (*dist) (Point* xv, Point* yv);

static bool operator== (Point& l, Point& r)
{
    Point::iterator l_it;
    Point::iterator r_it;
    for (l_it = l.begin(); l_it != l.end(); ++l_it) {
        r_it = r.find (l_it->first);
        if (r_it == r.end()) {
            return false;
        }
    }
    for (r_it = r.begin(); r_it != r.end(); ++r_it) {
        l_it = l.find (r_it->first);
        if (l_it == l.end()) {
            return false;
        }
    }
    return true;
}

//@{
// Euclidean distance
long double euclidean (Point* xv, Point* yv)
{
    volatile long double sum = 0.0;
    Point::iterator xv_it;
    Point::iterator yv_it;
    // match xs against ys
    for (xv_it = xv->begin(); xv_it != xv->end(); ++xv_it) {
        yv_it = yv->find (xv_it->first);
        if (yv_it != yv->end()) {
            sum += (xv_it->second - yv_it->second) *
                   (xv_it->second - yv_it->second);
        } else {
            sum += (xv_it->second * xv_it->second);
        }
    }
    // match ys against xs, only take non-matched account
    for (yv_it = yv->begin(); yv_it != yv->end(); ++yv_it) {
        if (xv->find (yv_it->first) == xv->end()) {
            sum += (yv_it->second * yv_it->second);
        }
    }
    return sum;
}
//@}

//@{
// Cosine distance
// NOTE we return the multiplicative inverse of (cosine plus one)
long double cosine (Point* xv, Point* yv)
{
    //FIXME correct return type
    volatile long double sum = 0.0, xs = 0.0, ys = 0.0;
    Point::iterator xv_it;
    Point::iterator yv_it;
    // inner product, <x,y>, <x,x>
    for (xv_it = xv->begin(); xv_it != xv->end(); ++xv_it) {
        yv_it = yv->find (xv_it->first);
        if (yv_it != yv->end()) {
            sum += (xv_it->second * yv_it->second);
        }
        xs += (xv_it->second * xv_it->second);
    }
    // inner product, <y,y>
    for (yv_it = yv->begin(); yv_it != yv->end(); ++yv_it) {
        ys += (yv_it->second * yv_it->second);
    }
    // 1/(1+cos<x,y>)
    return 1.0 / (1.0 + sum * sum / xs * ys);
}
//@}



//@{
long double distance (Map* freqTable_in, int x, int y)
{
    int x_, y_;
    long double z;
    if (x == y) return 0.0;
    if (x < y) {
        x_ = x, y_ = y;
    } else {
        x_ = y, y_ = x;
    }
// NOTE it's better to use table for fast result access but it eat too much RAM
//     if (!table->empty() && (table->find (x_) != table->end())) {
//         if (table->at (x_).find (y_) != table->at (x_).end()) {
//             return table->at (x_).at (y_);
//         } else {
//             //@{
//             z = dist (&(*freqTable_in)[x_], &(*freqTable_in)[y_]);
//             //@}
//             table->at (x_).insert (std::pair<int, int> (y_, z));
//             return z;
//         }
//     } /*else {*/
    // TODO get impl according to method_in
    //@{
//     z = euclidean (& (*freqTable_in) [x_], & (*freqTable_in) [y_]);
    z = dist (& (*freqTable_in) [x_], & (*freqTable_in) [y_]);
    //@}
//     std::map<int, int> e;
//     e.insert (std::pair<int, int> (y_, z));
//     table->insert (std::pair<int, std::map<int, int> > (x_, e));
    return z;
}

//@}

//@{
// the median is one that gets closest to others
// \min_i\max_j S_{i,j}
// O(n^2)
int median (Map* freqTable_in, std::list<int> vec)
{
    int x = 0/*, y = 0*/;
    long double max = 0.0;
    long double min = std::numeric_limits<long double>::max();
    long double e = 0.0;
    std::list<int>::iterator x_it;
    std::list<int>::iterator y_it;
    // if euclidean return 0,
    // that is if happen to table is not calculated already
    for (x_it = vec.begin(); x_it != vec.end(); ++x_it) {
        for (y_it = x_it; y_it != vec.end(); ++y_it) {
            e = distance (freqTable_in, *x_it, *y_it);
            if (max < e) {
                max = e;
            }
        }
        if (min > max) {
            min = max;
            x = *x_it;
        }
    }
    return x;
}
//@}

//@{

// TODO accept method argument
/**
 * @brief ...
 *
 * @param freqTable_in ...
 * @param index ...
 * @param clouds ...
 * @param cloudSize ...
 * @return int the center of that cloud
 **/
int whichClosest (Map* freqTable_in, int index,
                  std::map<int, cloud> *clouds, int cloudSize)
{
    // avoid empty cloud
    // TODO keep or move one from another cloud that has a long distance from
    // its center into this (*clouds)[index]
    if ( (*clouds) [index].dusts.size() == 1) return index;
    long double dist_t = 0.0, dist = std::numeric_limits<long double>::max();
    int which = index;
    // Performance bottle kneck
    for (int i = 0 ; i < cloudSize; ++i) {
        dist_t = distance (freqTable_in, index, (*clouds) [i].center_index);
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
inline void removeFromCenters (std::map<int, cloud> * clouds,
                               const int key, const int value)
{

    if (clouds && key > -1) {
        (*clouds) [key].dusts.remove (value);
    }/** else {
        std::cout << "removeFromCenters" << std::endl;
    }**/
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
inline void addToCenters (std::map<int, cloud> *  clouds,
                          int key, const int value)
{
    if (clouds) {
        (*clouds) [key].dusts.push_back (value);
    }/** else {
        std::cout << __LINE__ << std::endl;
    }**/
}
//@}


//@{
void updateCloudCenters (Map* freqTable_in, std::map<int, cloud> * clouds,
                         int cloudSize)
{
    for (int i = 0; i < cloudSize; ++i) {
        (*clouds) [i].center_index = median (freqTable_in, (*clouds) [i].dusts);
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
    // the product table
    // if (table) delete table;
    // table = new std::map<int, std::map<int, int> >;

    List points (points_in);
    unsigned int pointsSize = points.size();

    // int method = as<int> (method_in);
    switch (as<int> (method_in)) {
    case 1:
        dist = cosine;
        break;
    default:
        dist = euclidean;
    }
    int cloudSize = as<int> (clusterSize_in);

    //@{
    std::hash< std::string > hashFun;
    Map freqTable; // < index <<hash,count>,<hash,count>,.. > >
    std::vector<std::string>::iterator v_it;
    std::size_t v_t = 0LU;
    for (int i = 0; i < pointsSize; ++i) {
        Point v;
        std::vector<std::string> s = as<std::vector<std::string> > (points[i]);
        for (v_it = s.begin(); v_it != s.end(); ++v_it) {
            v_t = hashFun (*v_it);
            if (v.find (v_t) != v.end()) {
                v[v_t] ++;
            } else {
                v[v_t] = 1;
            }
        }
        freqTable[i] = v;
    }
    //@}
    //@{
    // initialize a cluster of size cluster
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
        tmp.push_back (*sh_it);
    }
    shuffle.clear();
    std::map<int, cloud> clouds;
    for (int i = 0; i < cloudSize; ++i) {
        std::list<int> l;
        l.push_back (tmp[i]);
        clouds[i].dusts = l;
        clouds[i].center_index = tmp[i];
        cluster[tmp[i]] = i;
//         std::cout << tmp[i] << std::endl;
    }
    tmp.clear();
    //@}

    volatile int which = 0;
    int iter = as<int> (iter_in); // assert positive in R level
    volatile int changed = 0;
    int c = 0;
    while (iter--) {
        changed = 0;
        c++;
        for (int i = 0; i < pointsSize; ++i) {
            // one may change cluster alternatively cos' none is actually near
            // to it
            which = whichClosest (&freqTable, i, &clouds, cloudSize);
            if (cluster[i] != which) {
                /* change cluster[i] state */
                /* remove */
                removeFromCenters (&clouds, cluster[i], i);
                /* update */
                cluster[i] = which;
                addToCenters (&clouds, which, i);
                // record which
                changed ++;
            }
        }
        if (changed == 0) break;
        // TODO it's better that only update which and cluster[i]
        // use a vector to store (which and cluster[i])
        updateCloudCenters (&freqTable, &clouds, cloudSize);
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
    return z;
    END_RCPP

}
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
