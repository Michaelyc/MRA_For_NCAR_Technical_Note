#pragma once

#include <armadillo>

using namespace arma;

/*
Get cross-covariance matrix between two location sets

Input:
    locationsX1: X coordinates of location set 1
    locationsY1: Y coordinates of location set 1
    locationsX2: X coordinates of location set 2
    locationsY2: Y coordinates of location set 2
    nLocations1: the number of locations in lcs1
    nLocations2: the number of locations in lcs2
    sill: sill parameter
    range: range parameter
    nugget: nugget parameter

Output:
    covarianceMatrix: obtained covariances
*/
void evaluate_cross_covariance(double *covarianceMatrix, double *locationsX1, double *locationsY1, double *locationsX2, double *locationsY2, const unsigned long &nLocations1, const unsigned long &nLocations2, const double &sill, const double &range, const double &nugget);

void evaluate_variance_covariance(double *covarianceMatrix, double *locationsX, double *locationsY, const unsigned long &nLocations, const double &sill, const double &range, const double &nugget);