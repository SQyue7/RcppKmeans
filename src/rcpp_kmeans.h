#ifndef _RCPP_KMEANS_H
#define _RCPP_KMEANS__H

#include <Rcpp.h>

RcppExport SEXP Kmeans(SEXP /*points*/, SEXP /*cluster*/,
                       SEXP /*iter*/, SEXP /*method*/, SEXP epsilon);

#endif
