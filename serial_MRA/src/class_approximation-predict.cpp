#include <iostream>
#include <armadillo>
#include <mkl.h>

#include "class_data.hpp"
#include "class_approximation.hpp"
#include "class_partition.hpp"
#include "evaluate_covariance.hpp"
#include "constants.hpp"
#include <string.h>


void Approximation::predict()
{
    gettimeofday(&timeNow, NULL);
    cout<<"===========================> Predicting starts. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n";
			

    int currentLevel = NUM_LEVELS_M-1;
    unsigned long indexRegionAtThisLevel = 0;
    int nKnots = partition->nKnots;
    char L='L', N='N';

    //Define indices of all the ancestors for the current region
    unsigned long *ancestorArray = new unsigned long [currentLevel];

    for(unsigned long iRegion = REGION_START[NUM_LEVELS_M-1]; iRegion < REGION_END[NUM_LEVELS_M-1] + 1; iRegion++)
    {
        int status;
        unsigned long indexRegionAtThisLevel = iRegion-partition->nRegionsInTotal+partition->nRegionsAtEachLevel[NUM_LEVELS_M-1];
        int nPredictionsInCurrentRegion = partition->nPredictionsAtFinestLevel[indexRegionAtThisLevel];
        if(nPredictionsInCurrentRegion == 0) continue;

        //Find the indices of all the ancestors
        get_all_ancestors(ancestorArray, iRegion, currentLevel);
        
        int jLevel = currentLevel-1;
        mat* KcBTilde;
        mat** tmpBTilde = new mat* [jLevel];
        
        //KcBTilde := RChol[ancestorArray[jLevel]])^(-1)*BTilde[indexRegionAtThisLevel][jLevel]
        KcBTilde = new mat(BTilde[indexRegionAtThisLevel][jLevel],nKnots,nPredictionsInCurrentRegion,true,true);
        dtrtrs_(&L,&N,&N,&nKnots,&nPredictionsInCurrentRegion,RChol[ancestorArray[jLevel]].memptr(),&nKnots,KcBTilde->memptr(),&nKnots,&status);
        
        for(int kLevelBeforejLevel = 0 ; kLevelBeforejLevel < jLevel; kLevelBeforejLevel++)
        {
            //tmpBTilde[kLevelBeforejLevel] := BTilde[indexRegionAtThisLevel][kLevelBeforejLevel].t()-KCholTimesCurrentA[ancestorArray[currentLevel-1]].slice(kLevelBeforejLevel).t()*KcBTilde
            tmpBTilde[kLevelBeforejLevel] = new mat(BTilde[indexRegionAtThisLevel][kLevelBeforejLevel],nKnots,nPredictionsInCurrentRegion,false,true);
            cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,nKnots,nPredictionsInCurrentRegion,nKnots,-1.0,KCholTimesCurrentA[ancestorArray[currentLevel-1]].slice(kLevelBeforejLevel).memptr(),nKnots,KcBTilde->memptr(),nKnots,1.0,tmpBTilde[kLevelBeforejLevel]->memptr(),nKnots);
        }
        delete KcBTilde;
        
        for(int jLevel = currentLevel-2; jLevel > 0; jLevel--)
        {
            //KcBTilde := RChol[ancestorArray[jLevel]])^(-1)*tmpBTilde[jLevel])
            KcBTilde = new mat(tmpBTilde[jLevel]->memptr(),nKnots,nPredictionsInCurrentRegion,true,true);
            dtrtrs_(&L,&N,&N,&nKnots,&nPredictionsInCurrentRegion,RChol[ancestorArray[jLevel]].memptr(),&nKnots,KcBTilde->memptr(),&nKnots,&status);
            
            for(int kLevelBeforejLevel = 0 ; kLevelBeforejLevel < jLevel; kLevelBeforejLevel++)
            {
                //tmpBTilde[kLevelBeforejLevel] := tmpBTilde[kLevelBeforejLevel]-KCholTimesCurrentA[ancestorArray[jLevel]].slice(kLevelBeforejLevel).t()*KcBTilde
                cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,nKnots,nPredictionsInCurrentRegion,nKnots,-1.0,KCholTimesCurrentA[ancestorArray[jLevel]].slice(kLevelBeforejLevel).memptr(),nKnots,KcBTilde->memptr(),nKnots,1.0,tmpBTilde[kLevelBeforejLevel]->memptr(),nKnots);
            }

            delete KcBTilde;
        }
        
        for(jLevel = 0; jLevel < currentLevel-1; jLevel++)
        {
            //KcBTilde := RChol[ancestorArray[jLevel]])^(-1)*tmpBTilde[jLevel])
            KcBTilde = new mat(tmpBTilde[jLevel]->memptr(),nKnots,nPredictionsInCurrentRegion,false,true);
            
            dtrtrs_(&L,&N,&N,&nKnots,&nPredictionsInCurrentRegion,RChol[ancestorArray[jLevel]].memptr(),&nKnots,KcBTilde->memptr(),&nKnots,&status);
            
            //posteriorPredictionMean[indexRegionAtThisLevel] += KcBTilde.t()*KCholTimesCurrentw[ancestorArray[jLevel]]
            cblas_dgemv(CblasColMajor,CblasTrans,nKnots,nPredictionsInCurrentRegion,1.0,KcBTilde->memptr(),nKnots,KCholTimesCurrentw[ancestorArray[jLevel]].memptr(),1,1.0,posteriorPredictionMean[indexRegionAtThisLevel].memptr(),1);

            //posteriorPredictionVariance[indexRegionAtThisLevel] += columnSum(KcBTilde.^2);
            posteriorPredictionVariance[indexRegionAtThisLevel] += sum(square(*KcBTilde),0).t();
        }

        jLevel=currentLevel-1;
        
        //KcBTilde := RChol[ancestorArray[jLevel]])^(-1)*BTilde[indexRegionAtThisLevel][jLevel])
        KcBTilde = new mat(BTilde[indexRegionAtThisLevel][jLevel],nKnots,nPredictionsInCurrentRegion,false,true);
        dtrtrs_(&L,&N,&N,&nKnots,&nPredictionsInCurrentRegion,RChol[ancestorArray[jLevel]].memptr(),&nKnots,KcBTilde->memptr(),&nKnots,&status);

        //posteriorPredictionMean[indexRegionAtThisLevel] += KcBTilde.t()*KCholTimesCurrentw[ancestorArray[jLevel]]
        cblas_dgemv(CblasColMajor,CblasTrans,nKnots,nPredictionsInCurrentRegion,1.0,KcBTilde->memptr(),nKnots,KCholTimesCurrentw[ancestorArray[jLevel]].memptr(),1,1.0,posteriorPredictionMean[indexRegionAtThisLevel].memptr(),1);
        
        //posteriorPredictionVariance[indexRegionAtThisLevel] += columnSum(KcBTilde.^2);
        posteriorPredictionVariance[indexRegionAtThisLevel] += sum(square(*KcBTilde),0).t();
        
        delete KcBTilde;
        for(int jLevel = 0 ; jLevel < currentLevel-1; jLevel++)
        {
            delete[] BTilde[indexRegionAtThisLevel][jLevel];
            delete tmpBTilde[jLevel];
        }
        delete[] BTilde[indexRegionAtThisLevel][currentLevel-1];
        delete[] tmpBTilde;
    }

    //Free memory for temporary variables
    delete[] ancestorArray;

    gettimeofday(&timeNow, NULL);
    cout<<"===========================> Predicting is complete. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n\n";
}
