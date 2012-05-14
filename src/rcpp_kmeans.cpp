#include "rcpp_kmeans.h"

using namespace Rcpp ;

//ISSUE
// which distance method is best to mesure texts

// TODO add more methods to calculate distances

typedef struct {
    int center_index; // one of dusts
    std::list<int> dusts; // we use remove
} cloud;

typedef std::map<std::string, int> Point;
typedef std::map<int, Point> Map;

// NOTE the table may be as big as i*j*8 bytes
// hmm, requires lots of RAM!!
//  <i ,<j, r> >
static std::map<int, std::map<int, int> > *table;

static double (*dist) (Point* xv, Point* yv);

//@{
// Euclidean distance
double euclidean (Point* xv, Point* yv)
{
    volatile double sum = 0;
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
// NOTE we return the multiplicative inverse of cosine
double cosine (Point* xv, Point* yv)
{
    //FIXME correct return type
    volatile double sum = 0.0, xs = 0.0, ys = 0.0;
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
    // NOTE we add 1 here to avoid sum == 1
    return (xs * ys / (sum * sum + 1.0));
}
//@}



//@{
double distance (Map* freqTable_in, int x, int y)
{
    int x_, y_, z;
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
    double max = 0;
    double min = 2147483647;
    int e = 0;
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
/**
 * @brief return the index of center that closest to @a points_in[@a index]
 *
 * @param points_in all points
 * @param index point's index
 * @param centers_in centers table
 * @return int the element position @a in centers_in
 **/
// TODO accept method argument
int whichClosest (Map* freqTable_in, int index,
                  std::map<int, cloud> *clouds, int cloudSize)
{
    int dist = 0x7FFFFFFF;
    int dist_t = 0;
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
                               int key, int value)
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
inline void addToCenters (std::map<int, cloud> *  clouds, int key, int value)
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
 * @brief cluster points into @a cluster using distance @a method
 *
 * @param points list of point, each point is a string vector
 * @param cluster the size of cluster
 * @param method the distance method to use
 * @return SEXP
 **/
SEXP Kmeans (SEXP points_in, SEXP clusterSize_in, SEXP iter_in, SEXP method_in)
{
    BEGIN_RCPP
    if (table) delete table;
    table = new std::map<int, std::map<int, int> >;

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
    Map freqTable;
    std::vector<std::string>::iterator v_it;
    for (int i = 0; i < pointsSize; ++i) {
        Point v;
        std::vector<std::string> s = as<std::vector<std::string> > (points[i]);
        for (v_it = s.begin(); v_it != s.end(); ++v_it) {
            // if found increase
            if (v.find (*v_it) != v.end()) {
                v[*v_it] += 1;
            } else {
                v[*v_it] = 1;
            }
        }
        freqTable[i] = v;
    }
    //@}

    //@{
    // for quickly find out which cluster a point belongs to
    // cluster to each points, -1 indicates not assigned yet
    std::vector<int> cluster (pointsSize);
    std::fill (cluster.begin(), cluster.end(), -1);
    // initialize a cluster of size cluster
    std::vector<int> tmp;
    for (int i = 0; i < pointsSize; ++i) {
        tmp.push_back (i);
    }

    std::random_shuffle (tmp.begin(), tmp.end());
    // TODO if two points close enough are assigned as centers, they might be
    // separated forever, so we'd better reject points that have center close to
    // it, regarding to a threshold.
    std::map<int, cloud> clouds;
    for (int i = 0; i < cloudSize; ++i) {
        std::list<int> l;
        l.push_back (tmp[i]);
        clouds[i].dusts = l;
        clouds[i].center_index = tmp[i];
        cluster[tmp[i]] = i;
    }
    //@}

    int which = 0;
    bool fixed = true;
    int iter = as<int> (iter_in); // assert positive in R level
    while (iter--) {
        fixed = true;
        for (int i = 0; i < pointsSize; ++i) {
            which = whichClosest (&freqTable, i, &clouds, cloudSize);
            if (cluster[i] != which) {
                /* NOTE is it better to exit if there is only one point */
                /* change cluster[i] state */
                /* remove */
                removeFromCenters (&clouds, cluster[i], i);
                /* update */
                cluster[i] = which;
                addToCenters (&clouds, which, i);
                // record which
                fixed = false;
            }
        }
        if (fixed) break;
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
    delete table;
    return z;
    END_RCPP

}
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
