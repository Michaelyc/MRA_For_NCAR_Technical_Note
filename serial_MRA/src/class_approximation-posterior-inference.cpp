#include <iostream>
#include <armadillo>
#include <mkl.h>

#include "class_data.hpp"
#include "class_approximation.hpp"
#include "class_partition.hpp"
#include "evaluate_covariance.hpp"
#include "constants.hpp"
#include <string.h>

void Approximation::posterior_inference()
{
    gettimeofday(&timeNow, NULL);
    cout<<"===========================> Calculating posterior quantities starts. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n";

    //Allocate memory for KCholTimesCurrentw (vector size: partition->nKnots)
    KCholTimesCurrentw = new vec [partition->nRegionsInTotal];

    //Allocate memory for KCholTimesCurrentA (cube size: first dimension: partition->nKnots, second dimesnion: partition->nKnots, third dimension or the number of slices: (partition->nRegionsInTotal-1)*partition->nRegionsInTotal/2)
    KCholTimesCurrentA = new cube [partition->nRegionsInTotal-partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]];

    //Allocate memory for temporary variables
    vec w;
    mat A;
    double *tempMemory = new double [partition->nKnots*partition->nKnots];
    
    w.set_size(partition->nKnots);
    A.set_size(partition->nKnots,partition->nKnots);

    //Temporary variables if load ATilde from disk 
    int size=(NUM_LEVELS_M-1)*NUM_LEVELS_M/2;
    unsigned long offset = partition->nRegionsInTotal-partition->nRegionsAtEachLevel[NUM_LEVELS_M-1];

    //From the second finest to the coarsest level
    for(int currentLevel = NUM_LEVELS_M-2; currentLevel > -1 ; currentLevel--)
    {
        for(unsigned long iRegion = REGION_START[currentLevel]; iRegion < REGION_END[currentLevel] + 1; iRegion++)
        {   
            unsigned long indexReginAtThisLevel = iRegion - REGION_START[currentLevel];

            //Deal with this region
            if(SAVE_TO_DISK_FLAG) load_ATilde_from_disk(iRegion*NUM_PARTITIONS_J+1, iRegion*NUM_PARTITIONS_J+NUM_PARTITIONS_J, currentLevel == NUM_LEVELS_M-2, partition->nKnotsAtFinestLevel, offset, partition->nKnots, size);

            aggregate_w_and_A(w, A, tempMemory, currentLevel, (currentLevel+1)*(currentLevel+2)/2-1, iRegion*NUM_PARTITIONS_J+1, iRegion*NUM_PARTITIONS_J+NUM_PARTITIONS_J, iRegion);

            if(CALCULATION_MODE != "prediction")
                loglikelihood += -2*sum(log(RChol[iRegion].diag()));

            //Compute the posterior conditional variance-covariance matrix
            RChol[iRegion] = chol( RChol[iRegion]*RChol[iRegion].t()+A,"lower");

            //Compute KCholTimesCurrentw that is the Cholesky factor of the inverse of the conditional posterior variance-covariance matrix times w
            KCholTimesCurrentw[iRegion].set_size(partition->nKnots);
            KCholTimesCurrentw[iRegion] = solve(trimatl(RChol[iRegion]),w);

            //Compute KCholTimesCurrentA that is the Cholesky factor of the inverse of the conditional posterior variance-covariance matrix times A
            KCholTimesCurrentA[iRegion].set_size(partition->nKnots,partition->nKnots,currentLevel);
        

            //Get wTilde and ATilde
            if(currentLevel != 0)
                get_wTilde_and_ATilde_in_posterior(w,A, tempMemory, partition->nKnots, KCholTimesCurrentA[iRegion], KCholTimesCurrentw[iRegion], iRegion, currentLevel);

            if(CALCULATION_MODE != "prediction")
                KCholTimesCurrentA[iRegion].reset();

            //Free wTilde and ATilde of the children to save memory
            
            free_wTilde_and_ATilde(iRegion*NUM_PARTITIONS_J+1, iRegion*NUM_PARTITIONS_J+NUM_PARTITIONS_J);
            
            if(CALCULATION_MODE != "prediction")
            {
                loglikelihood += 2*sum(log(RChol[iRegion].diag()))-dot(KCholTimesCurrentw[iRegion],KCholTimesCurrentw[iRegion]);
                KCholTimesCurrentw[iRegion].reset();
            }
            
            if(SAVE_TO_DISK_FLAG)
            {
                //save ATilde of the current region to disk
                std::string fileName=TMP_DIRECTORY+"/"+std::to_string(iRegion)+".bin";
                ATilde[iRegion].save(fileName,arma_binary);
                ATilde[iRegion].reset();
            }
            
        }
        
        timeval timeNow;
        gettimeofday(&timeNow, NULL);
        cout<<"Posterior: Level "<<currentLevel+1<<" is complete. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n";

    }

    if(CALCULATION_MODE!="prediction")
        cout<<"The obtained loglikelihood is: "<<-0.5*(loglikelihood+NUM_OBSERVATIONS*log(2*3.14159265359))<<endl;
    
    //Deallocate temporary memory
    w.reset();
    A.reset();
    delete[] tempMemory;

    gettimeofday(&timeNow, NULL);
    cout<<"===========================>  Calculating posterior quantities is complete. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n\n";
}