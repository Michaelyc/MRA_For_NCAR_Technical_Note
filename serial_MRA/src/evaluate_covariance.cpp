#include <armadillo>

#include "evaluate_covariance.hpp"

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
    covarianceMatrix: obtained covariances, of dimension nLocations2 * nLocations1
*/
void evaluate_cross_covariance(double *covarianceMatrix, double *locationsX1, double *locationsY1, double *locationsX2, double *locationsY2, const unsigned long &nLocations1, const unsigned long &nLocations2, const double &sill, const double &range, const double &nugget)
{
    for(unsigned long index1 = 0; index1 < nLocations1; index1++)
        for(unsigned long index2 = 0; index2 < nLocations2; index2++)
        {
            double xDifference=locationsX1[index1]-locationsX2[index2];
            double yDifference=locationsY1[index1]-locationsY2[index2];

            double distance = sqrt( xDifference*xDifference + yDifference*yDifference );
            if (distance < 1e-14)
                covarianceMatrix[index1*nLocations2+index2] = sill + nugget;
            else
                covarianceMatrix[index1*nLocations2+index2] = sill*exp(-distance/range);
        }
}

/*
Get variance-covariance matrix for one location set

Input:
    locationsX: X coordinates
    locationsY: Y coordinates
    nLocations: the number of locations
    sill: sill parameter
    range: range parameter
    nugget: nugget parameter

Output:
    covarianceMatrix: obtained covariances, of dimension nLocations * nLocations
*/
void evaluate_variance_covariance(double *covarianceMatrix, double *locationsX, double *locationsY, const unsigned long &nLocations, const double &sill, const double &range, const double &nugget)
{
    double offDiagonal;
    double diagonal = sill + nugget;
    for(unsigned long index1 = 0; index1 < nLocations; index1++)
    {
        for(unsigned long index2 = 0; index2 < index1; index2++)
        {
            double xDifference=locationsX[index1]-locationsX[index2];
            double yDifference=locationsY[index1]-locationsY[index2];

            double distance = sqrt( xDifference*xDifference + yDifference*yDifference );
            if (distance < 1e-14)
                offDiagonal = sill + nugget;
            else 
                offDiagonal = sill*exp(-distance/range);

            covarianceMatrix[index1*nLocations+index2] = offDiagonal;
            covarianceMatrix[index2*nLocations+index1] = offDiagonal;
        }
        covarianceMatrix[index1*nLocations+index1] = diagonal;
    }
}

