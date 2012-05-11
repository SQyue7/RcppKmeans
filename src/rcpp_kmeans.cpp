#include "rcpp_kmeans.h"
using namespace Rcpp ;

// TODO some methods to calculate distances

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
    return 0;// BUG if you have built up the table
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
//     std::cout << "TICK " << count << std::endl;
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
 * @return int
 **/
// TODO accept method argument
int whichClosest (List points_in, int index,
                  std::map<int, std::list<int> > *centers_in)
{
    int dist = 2147483647;
    int dist_t = 0;
    int which = index;
    std::map<int, std::list<int> >::iterator it;
    for (it = centers_in->begin(); it != centers_in->end(); ++it) {
        dist_t = euclidean (points_in, index, it->first);
        if (dist_t < dist) {
            dist = dist_t;
            which = it->first;
        }
    }
    return which;
}
//@}


//@{
inline void removeFromCenters (const std::map<int, std::list<int> > *centers,
                               int key, int value)
{
    std::list<int> &  p = (std::list<int> &) centers->at (key);
    if (!p.empty()) {
        p.remove (value);
    } else {
        std::cout << "removeFromCenters" << std::endl;
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
inline void addToCenters (std::map<int, std::list<int> > *centers,
                          int key, int value)
{
    std::list<int> * p = & (centers->at (key));
    p->push_back (value);
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
    // make a dictionary from points
    List points (points_in);
    unsigned int pointsSize = points.size();
    // TODO pass to distance dispatcher, not use right now, euclidean hard-coded
    // int method = as<int> (method_in);
    (void) method_in;


    std::map<int, std::list<int> > centers;
    std::map<int, std::list<int> >::iterator centers_it;
    std::map<int, std::list<int> >::iterator centers_it_ref;
    int k = as<int> (clusterSize_in);
    // initialize a cluster of size cluster
    //@{

    std::vector<int> tmp;
    for (int i = 0; i < pointsSize; ++i) {
        tmp.push_back (i);
    }
    std::random_shuffle (tmp.begin(), tmp.end());
    for (int i = 0; i < k; ++i) {
        centers[tmp[i]] = std::list<int>();
    }
    //@}
    // TODO fill centers

    // cluster to each points, -1 indicates not assigned yet
    std::vector<int> cluster (pointsSize);
    std::fill (cluster.begin(), cluster.end(), -1);

    int first = 0;
    int which = 0;
    bool fixed = true;
    int iter = as<int> (iter_in); // assert positive in R level
    while (iter--) {
        for (int i = 0; i < pointsSize; ++i) {
            which = whichClosest (points, i, &centers);
            if (cluster[i] != which) {
                /* change cluster[i] state */
                /* remove old state */
                if (cluster[i] != -1)
                    removeFromCenters (&centers, cluster[i], i);
                /* update state */
                cluster[i] = which;
                addToCenters (&centers, which, i);
                fixed = false;
            }
        }
        std::cout << "X" << std::endl;
        // exits if condition meets
        if (fixed) break;
        fixed = true;
        // update each center
        // BUG
        // FIXME erase may change iterator
        // FIXME allocation issue
        centers_it =  centers.begin();
        while (centers_it != centers.end()) {
            centers_it_ref = centers_it;
            centers_it ++;
            if (!centers_it_ref->second.empty()) {
                first = median (centers_it_ref->second);
                if (centers_it_ref->first != first) {
                    std::swap<std::list<int> > (centers[first], centers_it_ref->second);
                    centers.erase (centers_it_ref);
                }
            }
        }
    }
    std::cout << "XX" << std::endl;

    //TODO output cluster centers  (as<int> iter_in - iter )
    // clean up!!
    List centers_vec;
    for (centers_it =  centers.begin();
            centers_it != centers.end();
            ++centers_it) {
        std::vector<int> vec (centers_it->second.begin(), centers_it->second.end());
        centers_vec[centers_it->first] = vec;
    }


    List z =  List::create (_["clusters"] = cluster,
                            _["centers"] = centers_vec,
                            _["iter"] = as<int> (iter_in) - iter);
    delete table;
    return z;
    END_RCPP

}



// kate: indent-mode cstyle; indent-width 4; replace-tabs on;

