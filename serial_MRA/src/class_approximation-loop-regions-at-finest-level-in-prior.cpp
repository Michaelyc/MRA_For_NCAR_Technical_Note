#include <iostream>
#include <armadillo>
#include <mkl.h>

#include "class_data.hpp"
#include "class_approximation.hpp"
#include "class_partition.hpp"
#include "evaluate_covariance.hpp"
#include "constants.hpp"
#include <string.h>

void Approximation::loop_regions_at_finest_level_in_prior(double ***KCholTimesw)
{
    int currentLevel = NUM_LEVELS_M - 1;
    double loglikelihoodThisLevel = 0;
    
    unsigned long maxNumKnotsAtFinestLevel = 0;
    //Find the maximum number of knots at the finest level
    for(unsigned long iRegion = partition->nRegionsInTotal-partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion < partition->nRegionsInTotal; iRegion++)
        maxNumKnotsAtFinestLevel =  std::max(maxNumKnotsAtFinestLevel, partition->nKnotsAtFinestLevel[iRegion-partition->nRegionsInTotal+partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]]);

    /////// Define and allocate variables

    //Allocate memory for ATilde
    if(!SAVE_TO_DISK_FLAG)
        for(unsigned long iRegion = partition->nRegionsInTotal-partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion < partition->nRegionsInTotal; iRegion++)
            ATilde[iRegion].set_size(partition->nKnots, partition->nKnots, (NUM_LEVELS_M-1)*NUM_LEVELS_M/2);
    
    //Allocate memory for wTilde
    for(unsigned long iRegion = partition->nRegionsInTotal-partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion < partition->nRegionsInTotal; iRegion++)
        wTilde[iRegion].set_size(partition->nKnots,NUM_LEVELS_M-1);


    //Define indices of all the ancestors for the current region
    unsigned long *ancestorArray = new unsigned long [currentLevel];

    //Auxilliary variables for matrix computations
    int one=1;
    char L='L', N='N';
    int nKnots = partition->nKnots;
    int status;

    //Loop the regions 
    for(unsigned long iRegion = partition->nRegionsInTotal-partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion < partition->nRegionsInTotal; iRegion++)
    {
        mat **SicB;
        vec *Sicy;
        double *RCholThisLevel;

        unsigned long indexRegionAtThisLevel = iRegion-partition->nRegionsInTotal+partition->nRegionsAtEachLevel[NUM_LEVELS_M-1];
        
        //Assign the number of knots in the current region accordingly
        int nKnotsInCurrentRegion = partition->nKnotsAtFinestLevel[indexRegionAtThisLevel];

        if(nKnotsInCurrentRegion == 0) 
        {
            ATilde[iRegion].zeros();
            wTilde[iRegion].zeros();
            if(CALCULATION_MODE=="prediction") get_all_ancestors(ancestorArray, iRegion, currentLevel);
        }else
        {
            //Find the indices of all the ancestors
            get_all_ancestors(ancestorArray, iRegion, currentLevel);

            //Allocate memory for variables in the current region
            mat **covarianceMatrix = new mat* [currentLevel+1];

            double **covarianceMatrixMemory = new double* [currentLevel+1];
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
                covarianceMatrixMemory[jLevel] = new double [nKnotsInCurrentRegion*partition->nKnots];
            covarianceMatrixMemory[currentLevel] =  new double [nKnotsInCurrentRegion*nKnotsInCurrentRegion];

            SicB = new mat* [currentLevel];
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
                SicB[jLevel] = new mat(partition->nKnots,nKnotsInCurrentRegion);

            if(CALCULATION_MODE=="prediction") KCholTimesw[iRegion] = new double* [currentLevel];

            //Define and Allocate RCholThisLevel for the current region
            RCholThisLevel = new double [nKnotsInCurrentRegion*nKnotsInCurrentRegion];

            //Loop for all ancestors to get the cross-covariance matrix
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
            {
                //Calculate the covariances between current region and the ancestors
                covarianceMatrix[jLevel] =  new mat(covarianceMatrixMemory[jLevel],nKnotsInCurrentRegion,nKnots,false,true);
                evaluate_cross_covariance(covarianceMatrixMemory[jLevel], partition->knotsX[iRegion], partition->knotsY[iRegion], partition->knotsX[ancestorArray[jLevel]], partition->knotsY[ancestorArray[jLevel]], nKnotsInCurrentRegion, nKnots, sill, range, nugget);
            }

            //Get the current variance-covariance matrix
            covarianceMatrix[currentLevel] =  new mat(covarianceMatrixMemory[currentLevel],nKnotsInCurrentRegion,nKnotsInCurrentRegion,false,true);
            evaluate_variance_covariance(covarianceMatrixMemory[currentLevel], partition->knotsX[iRegion], partition->knotsY[iRegion], nKnotsInCurrentRegion, sill, range, nugget);
            
            //Loop for all ancestors to get the conditional cross-covariance matrix
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
            {
                get_conditional_covariance_matrix(jLevel, covarianceMatrixMemory, KCholTimesw[ancestorArray[jLevel]], covarianceMatrixMemory[jLevel], nKnotsInCurrentRegion,nKnots,nKnots);
                
                //Copy covarianceMatrix[jLevel] to SicB[jLevel] for later use
                memcpy(SicB[jLevel]->memptr(),covarianceMatrixMemory[jLevel],nKnotsInCurrentRegion*nKnots*sizeof(double));
                
                //Reassign covarianceMatrix[jLevel] = RChol[ancestorArray[jLevel]]^(-1)*covarianceMatrix[jLevel]
                dtrtrs_(&L,&N,&N,&nKnots,&nKnotsInCurrentRegion,RChol[ancestorArray[jLevel]].memptr(),&nKnots,covarianceMatrixMemory[jLevel],&nKnots,&status);

                if(CALCULATION_MODE=="prediction") KCholTimesw[iRegion][jLevel]=covarianceMatrixMemory[jLevel];
            }
            //Get the conditional variance-covariance matrix of this region
            get_conditional_covariance_matrix(currentLevel, covarianceMatrixMemory, covarianceMatrixMemory, covarianceMatrixMemory[currentLevel], nKnotsInCurrentRegion,nKnotsInCurrentRegion,nKnots);      
            
            //Cholesky factorization, covarianceMatrix[currentLevel] = RCholThisLevel*RCholThisLevel^T
            memcpy(RCholThisLevel,covarianceMatrixMemory[currentLevel],nKnotsInCurrentRegion*nKnotsInCurrentRegion*sizeof(double));
            dpotrf_(&L,&nKnotsInCurrentRegion,RCholThisLevel,&nKnotsInCurrentRegion,&status);   

            //Compute Sicy = RCholThisLevel^(-1)*partition->knotsResidual[indexRegionAtThisLevel]
            Sicy = new vec(partition->knotsResidual[indexRegionAtThisLevel],nKnotsInCurrentRegion,true,true);
            dtrtrs_(&L,&N,&N,&nKnotsInCurrentRegion,&one,RCholThisLevel,&nKnotsInCurrentRegion,Sicy->memptr(),&nKnotsInCurrentRegion,&status);
            
            //Compute SicB[jLevel] = RCholThisLevel^(-1)*covarianceMatrix[jLevel]
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
            {      
                inplace_trans(*SicB[jLevel]);
                dtrtrs_(&L,&N,&N,&nKnotsInCurrentRegion,&nKnots,RCholThisLevel,&nKnotsInCurrentRegion,SicB[jLevel]->memptr(),&nKnotsInCurrentRegion,&status);
            }

            get_wTilde_and_ATilde_in_prior(nKnots, currentLevel, iRegion, Sicy, SicB);

            if(CALCULATION_MODE!="prediction")
            {
                double loglikelihoodToAdd = 0; int pos = 0;
                for(int i = 0; i < nKnotsInCurrentRegion; i++)
                {
                    loglikelihoodToAdd += (*Sicy)[i]*(*Sicy)[i] + 2*log(RCholThisLevel[pos]);
                    pos += nKnotsInCurrentRegion + 1;
                }

                loglikelihoodThisLevel += loglikelihoodToAdd;

                for(int jLevel = 0; jLevel < currentLevel; jLevel++)
                {
                    delete SicB[jLevel];
                    delete[] covarianceMatrixMemory[jLevel];
                    delete covarianceMatrix[jLevel];
                }

                delete[] RCholThisLevel;
            }

            delete[] covarianceMatrixMemory[currentLevel];
            delete covarianceMatrix[currentLevel];
        }

        if(CALCULATION_MODE=="prediction")
        {
            int nPredictionsInCurrentRegion = partition->nPredictionsAtFinestLevel[indexRegionAtThisLevel];
            if(nPredictionsInCurrentRegion == 0) continue;

            //Allocate memory for variables in the current region
            mat **covarianceMatrix = new mat* [currentLevel+1];

            double **covarianceMatrixMemory = new double* [currentLevel];
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
                covarianceMatrixMemory[jLevel] = new double [nPredictionsInCurrentRegion*nKnots];

            double **posteriorPredictionKCholTimesw = new double* [currentLevel+1];
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
                posteriorPredictionKCholTimesw[jLevel]= new double [nPredictionsInCurrentRegion*nKnots];

            if(nKnotsInCurrentRegion > 0)
                posteriorPredictionKCholTimesw[currentLevel]= new double [nPredictionsInCurrentRegion*nKnotsInCurrentRegion];            

            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
            {
                //Calculate the covariances between prediction locations in the current region and knots in the ancestor regions
                covarianceMatrix[jLevel] =  new mat(covarianceMatrixMemory[jLevel],nKnots,nPredictionsInCurrentRegion,false,true);

                evaluate_cross_covariance(covarianceMatrixMemory[jLevel],  partition->predictionX[indexRegionAtThisLevel], partition->predictionY[indexRegionAtThisLevel], partition->knotsX[ancestorArray[jLevel]], partition->knotsY[ancestorArray[jLevel]], nPredictionsInCurrentRegion, nKnots, sill, range, nugget);
                
                get_conditional_covariance_matrix(jLevel, posteriorPredictionKCholTimesw, KCholTimesw[ancestorArray[jLevel]], covarianceMatrixMemory[jLevel],nPredictionsInCurrentRegion,nKnots,nKnots);

                //posteriorPredictionKCholTimesw[jLevel] = RChol[ancestorArray[jLevel]]^(-1)*covarianceMatrix[jLevel]
                memcpy(posteriorPredictionKCholTimesw[jLevel],covarianceMatrixMemory[jLevel],nKnots*nPredictionsInCurrentRegion*sizeof(double));
                dtrtrs_(&L,&N,&N,&nKnots,&nPredictionsInCurrentRegion,RChol[ancestorArray[jLevel]].memptr(),&nKnots,posteriorPredictionKCholTimesw[jLevel],&nKnots,&status);

            }

            if(nKnotsInCurrentRegion > 0)
            {
                //Calculate the covariances between prediction locations and the knots in the current region
                evaluate_cross_covariance(posteriorPredictionKCholTimesw[currentLevel],  partition->predictionX[indexRegionAtThisLevel], partition->predictionY[indexRegionAtThisLevel], partition->knotsX[iRegion], partition->knotsY[iRegion], nPredictionsInCurrentRegion, nKnotsInCurrentRegion, sill, range, nugget);

                get_conditional_covariance_matrix(currentLevel, posteriorPredictionKCholTimesw, KCholTimesw[iRegion], posteriorPredictionKCholTimesw[currentLevel],nPredictionsInCurrentRegion,nKnotsInCurrentRegion,nKnots);

                //posteriorPredictionKCholTimesw[currentLevel] := RCholThisLevel^(-1)*posteriorPredictionKCholTimesw[currentLevel]
                dtrtrs_(&L,&N,&N,&nKnotsInCurrentRegion,&nPredictionsInCurrentRegion,RCholThisLevel,&nKnotsInCurrentRegion,posteriorPredictionKCholTimesw[currentLevel],&nKnotsInCurrentRegion,&status);
            }

            //Get the variance-covariance matrix of the prediction locations in the current region
            mat *varianceCovarianceMatrix = new mat(nPredictionsInCurrentRegion,nPredictionsInCurrentRegion);

            evaluate_variance_covariance((*varianceCovarianceMatrix).memptr(),  partition->predictionX[indexRegionAtThisLevel], partition->predictionY[indexRegionAtThisLevel], nPredictionsInCurrentRegion, sill, range, nugget);

            get_conditional_covariance_matrix(currentLevel, posteriorPredictionKCholTimesw, posteriorPredictionKCholTimesw, (*varianceCovarianceMatrix).memptr(), nPredictionsInCurrentRegion, nPredictionsInCurrentRegion, nKnots);

            posteriorPredictionMean[indexRegionAtThisLevel].set_size(nPredictionsInCurrentRegion);

            if(nKnotsInCurrentRegion > 0)
            {
                cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,nPredictionsInCurrentRegion,nPredictionsInCurrentRegion,nKnotsInCurrentRegion,-1.0,posteriorPredictionKCholTimesw[currentLevel],nKnotsInCurrentRegion,posteriorPredictionKCholTimesw[currentLevel],nKnotsInCurrentRegion,1.0,(*varianceCovarianceMatrix).memptr(),nPredictionsInCurrentRegion);

                cblas_dgemv(CblasColMajor, CblasTrans, nKnotsInCurrentRegion, nPredictionsInCurrentRegion, 1.0, posteriorPredictionKCholTimesw[currentLevel], nKnotsInCurrentRegion, Sicy->memptr(), 1, 0.0, posteriorPredictionMean[indexRegionAtThisLevel].memptr(), 1);
            }
            else
                memset(posteriorPredictionMean[indexRegionAtThisLevel].memptr(),0,sizeof(double)*nPredictionsInCurrentRegion);

            posteriorPredictionVariance[indexRegionAtThisLevel] = (*varianceCovarianceMatrix).diag();

            //BTilde[indexRegionAtThisLevel][jLevel] := covarianceMatrix[jLevel]-posteriorPredictionKCholTimesw[currentLevel].t()*SicB[jLevel]
            BTilde[indexRegionAtThisLevel] = covarianceMatrixMemory;

            if(nKnotsInCurrentRegion > 0)
                for(int jLevel = 0; jLevel < currentLevel; jLevel++)
                    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,nKnots,nPredictionsInCurrentRegion,nKnotsInCurrentRegion,-1.0,SicB[jLevel]->memptr(),nKnotsInCurrentRegion,posteriorPredictionKCholTimesw[currentLevel],nKnotsInCurrentRegion,1.0,BTilde[indexRegionAtThisLevel][jLevel],nKnots);

            delete varianceCovarianceMatrix;
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
                delete[] posteriorPredictionKCholTimesw[jLevel];

            if(nKnotsInCurrentRegion > 0)
            {
                delete[] RCholThisLevel;
                delete[] posteriorPredictionKCholTimesw[currentLevel];
            }
        }

        if(nKnotsInCurrentRegion > 0) delete Sicy;
    }

    if(CALCULATION_MODE!="prediction") loglikelihood += loglikelihoodThisLevel;

    //Free memory for ancestorArray
    delete[] ancestorArray;
}