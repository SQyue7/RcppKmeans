#include "rcpp_kmeans.h"
using namespace Rcpp ;

//ISSUE
// 1. which distance method is best to mesure texts
// 2. is it better to build a auxiliary dictionary

// TODO some methods to calculate distances



typedef struct {
    int center_index; // one of dusts
    std::list<int> dusts; // we use remove
} cloud;

typedef cloud* cloudAccess;

//  <i ,<j, r> >
static std::map<int, std::map<int, int> > *table;


//@{
// Euclidean distance


int euclidean (List points, int x, int y)
{
    int x_, y_, z;
    if (x < y) {
        x_ = x, y_ = y;
    } else {
        x_ = y, y_ = x;
    }

    if (!table->empty() && (table->find (x_) != table->end())) {
        if (table->at (x_).find (y_) != table->at (x_).end()) {
            return table->at (x_).at (y_);
        } else {
            // FIXME the following is .WRONG.
            std::vector<std::string> xs = as<std::vector<std::string> > (points[x_]);
            std::vector<std::string> ys = as<std::vector<std::string> > (points[y_]);
            std::vector<std::string>::iterator x_it, y_it;
            volatile int sum = 0;
            volatile int counter = 0;
            for (x_it = xs.begin(); x_it != xs.end(); ++x_it) {
                counter = 0;
                for (y_it = ys.begin(); y_it != ys.end(); ++y_it) {
                    if (not x_it->compare (*y_it)) {
                        counter++;
                    }
                }
                sum += counter * counter;
            }
            z = sum;
            table->at (x_).insert (std::pair<int, int> (y_, z));
            return z;
        }
    } else {
        std::vector<std::string> xs = as<std::vector<std::string> > (points[x_]);
        std::vector<std::string> ys = as<std::vector<std::string> > (points[y_]);
        std::vector<std::string>::iterator x_it, y_it;
        volatile int sum = 0;
        volatile int counter = 0;
        for (x_it = xs.begin(); x_it != xs.end(); ++x_it) {
            counter = 0;
            for (y_it = ys.begin(); y_it != ys.end(); ++y_it) {
                if (x_it->compare (*y_it)) {
                    counter++;
                }
            }
            sum += counter * counter;
        }
        z = sum;
        std::map<int, int> e;
        e.insert (std::pair<int, int> (y_, z));
        table->insert (std::pair<int, std::map<int, int> > (x_, e));
        return z;
    }
}

int euclidean (int x, int y)
{
    int x_, y_, z;
    if (x < y) {
        x_ = x, y_ = y;
    } else {
        x_ = y, y_ = x;
    }
    if (table->find (x_) != table->end()) {
        if (table->at (x_).find (y_) != table->at (x_).end()) {
            return table->at (x_).at (y_);
        }
    }
    return 0;
}
//@}

//@{
// the median is one that gets closest to others
// \min_i\max_j S_{i,j}
// O(n^2)
int median (std::list<int> vec)
{
    int x = 0, y = 0;
    int max = 0;
    int min = 2147483647;
    int e = 0;
    std::list<int>::iterator x_it;
    std::list<int>::iterator y_it;
    int count = 0;
    for (x_it = vec.begin(); x_it != vec.end(); ++x_it) {
        for (y_it = x_it; y_it != vec.end(); ++y_it) {
            e = euclidean (*x_it, *y_it);
            if (max < e) {
                max = e;
                y = *y_it;
            }
            count++;
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
int whichClosest (List points_in, int index, cloudAccess clouds, int cloudSize)
{
    int dist = 2147483647;
    int dist_t = 0;
    int which = index;
    for (int i = 0 ; i < cloudSize; ++i) {
        dist_t = euclidean (points_in, index, clouds[i].center_index);
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
 * @brief ...
 *
 * @param clouds ...
 * @param key the center index
 * @param value the dust that was included
 * @return void
 **/
inline void removeFromCenters (cloudAccess clouds, int key, int value)
{

    if (clouds && key > -1) {
        clouds[key].dusts.remove (value);
    } else {
//         std::cout << "removeFromCenters" << std::endl;
    }
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
inline void addToCenters (cloudAccess clouds, int key, int value)
{
    if (clouds) {
        clouds[key].dusts.push_back (value);
    } else {
        std::cout << "XXXXXX" << std::endl;
    }
}
//@}


//@{
void updateCloudCenters (cloudAccess clouds, int cloudSize)
{
    for (int i = 0; i < cloudSize; ++i) {
        (clouds[i]).center_index = median ( (clouds[i]).dusts);
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
    // TODO pass to distance dispatcher, not use right now, euclidean hard-coded
    // int method = as<int> (method_in);
    (void) method_in;

    int cloudSize = as<int> (clusterSize_in);
    cloudAccess clouds = new cloud[cloudSize];

    // for fast determine which cluster a point belongs to
    //@{
    // cluster to each points, -1 indicates not assigned yet
    std::vector<int> cluster (pointsSize);
    std::fill (cluster.begin(), cluster.end(), -1);
    // initialize a cluster of size cluster
    std::vector<int> tmp;
    for (int i = 0; i < pointsSize; ++i) {
        tmp.push_back (i);
    }
    std::random_shuffle (tmp.begin(), tmp.end());
    for (int i = 0; i < cloudSize; ++i) {
        (clouds[i]).center_index = tmp[i];
        (clouds[i]).dusts.push_back (tmp[i]);
        cluster[tmp[i]] = i;
    }
    //@}

    int which = 0;
    bool fixed = true;
    int iter = as<int> (iter_in); // assert positive in R level
    while (iter--) {
        for (int i = 0; i < pointsSize; ++i) {
            which = whichClosest (points, i, clouds, cloudSize);
//             std::cout << "XXXXXX " << __LINE__  << std::endl;
            if (cluster[i] != which) {
                /* change cluster[i] state */
                /* remove old state */
//                 std::cout << "XXXXXX " << __LINE__  << std::endl;
                removeFromCenters (clouds, cluster[i], i);
//                 std::cout << "XXXXXX  " << __LINE__  << std::endl;
                /* update state */
                cluster[i] = which;
                addToCenters (clouds, which, i);
//                 std::cout << "XXXXXX   " << __LINE__  << std::endl;
                fixed = false;
            }
        }
        if (fixed) break;
        fixed = true;
        updateCloudCenters (clouds, cloudSize);
//         std::cout << "XXXXXX    " << __LINE__  << std::endl;
    }
    List x = List::create (IntegerVector(), cloudSize);
    for (int i = 0; i < cloudSize; ++i) {
        IntegerVector nv (clouds[i].dusts.begin(), clouds[i].dusts.end());
        x[i] = nv;
    }
    List z =  List::create (_["clusters"] = x,
                            _["iterations"] = as<int> (iter_in) - iter);
    delete table;
    return z;
    END_RCPP

}
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
