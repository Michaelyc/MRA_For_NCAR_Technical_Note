#include <iostream>
#include <armadillo>
#include <mkl.h>

#include "class_data.hpp"
#include "class_approximation.hpp"
#include "class_partition.hpp"
#include "evaluate_covariance.hpp"
#include "constants.hpp"
#include <string.h>

/////////////////// Below are private member functions /////////////////

void Approximation::prior_allocate_memory()
{
    //Allocate memory for RChol
    RChol = new mat [partition->nRegionsInTotal];

    //Allocate memory for wTilde (first dimension: partition.nKnots, second dimesnion: partition.nLvs-1)
    wTilde = new mat [partition->nRegionsInTotal];

    //Allocate memory for ATilde (first dimension: partition.nKnots, second dimesnion: partition.nKnots, third dimension or the number of slices: (partition.nLvs-1)*partition.nLvs/2)
    ATilde = new cube [partition->nRegionsInTotal];
    
    //Assign paramter values if predicting
    if(CALCULATION_MODE == "prediction")
    {
        posteriorPredictionMean = new vec [partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]];
        posteriorPredictionVariance = new vec [partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]];
        BTilde = new double** [partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]];
    }
}

void Approximation::get_all_ancestors(unsigned long* ancestorArray, unsigned long currentRegion, const unsigned long &nAncestors)
{
    for(unsigned long iAncestor = 0; iAncestor < nAncestors; iAncestor++)
    {
        currentRegion = (currentRegion-1)/NUM_PARTITIONS_J;
        ancestorArray[nAncestors-iAncestor-1] = currentRegion;
    }
}

//Update: Changed the column-row for KcW and covariance matrix. Need to modify prediction routines accordingly
void Approximation::get_conditional_covariance_matrix(const int &ancestorLevel, double **ancestorNow, double **ancestorBefore, double* nowBefore, const int &numKnotsNow, const int &numKnotsBefore, const int &numKnotsAncestor)
{
    for(int iLevel = 0; iLevel < ancestorLevel; iLevel++)
        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,numKnotsBefore,numKnotsNow,numKnotsAncestor,-1.0,ancestorBefore[iLevel],numKnotsAncestor,ancestorNow[iLevel],numKnotsAncestor,1.0,nowBefore,numKnotsBefore);
}

void Approximation::get_wTilde_and_ATilde_in_prior(const int &nKnots, const int &currentLevel, const unsigned long &iRegion, vec *Sicy, mat **SicB)
{
    cube* ATildeCurrent;

    if(SAVE_TO_DISK_FLAG)
    {
        ATildeCurrent = new cube;
        (*ATildeCurrent).set_size(nKnots,nKnots,(NUM_LEVELS_M-1)*NUM_LEVELS_M/2);
    }

    unsigned long count = 0;
    for(int jLevel = 0; jLevel < currentLevel; jLevel++)
    {
        mat tmp = (*SicB[jLevel]).t();
        wTilde[iRegion].unsafe_col(jLevel) = tmp*(*Sicy);
        for(int kLevelAfterjLevel = jLevel; kLevelAfterjLevel < currentLevel; kLevelAfterjLevel++)
        {
            if(SAVE_TO_DISK_FLAG)
                (*ATildeCurrent).slice(count++) = tmp*(*SicB[kLevelAfterjLevel]);
            else
                ATilde[iRegion].slice(count++) = tmp*(*SicB[kLevelAfterjLevel]);
        }
    }

    if(SAVE_TO_DISK_FLAG)
    {
        std::string fileName=TMP_DIRECTORY+"/"+std::to_string(iRegion)+".bin";
        ATildeCurrent->save(fileName,arma_binary);
        delete ATildeCurrent;
    }
}

void Approximation::loop_regions_before_finest_level_in_prior(double ***KCholTimesw)
{
    unsigned long *ancestorArray = new unsigned long [NUM_LEVELS_M-2];
    
    //Allocate memory for the current region of KCholTimesw
    for(int currentLevel = 0; currentLevel < NUM_LEVELS_M-1; currentLevel++)
    {
        for(unsigned long iRegion = REGION_START[currentLevel]; iRegion < REGION_END[currentLevel] + 1; iRegion++)
        {
            KCholTimesw[iRegion] = new double* [currentLevel];
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
                KCholTimesw[iRegion][jLevel]= new double [partition->nKnots*partition->nKnots];

            RChol[iRegion].set_size(partition->nKnots,partition->nKnots);
        }
    }

    //Auxilliary variables for matrix computations
    char L='L', N='N';
    int nKnots = partition->nKnots;

    for(int currentLevel = 0; currentLevel < NUM_LEVELS_M-1; currentLevel++)
    {
        for(unsigned long iRegion = REGION_START[currentLevel]; iRegion < REGION_END[currentLevel] + 1; iRegion++)
        {
            int status;
            //Find the indices of all the ancestors

            get_all_ancestors(ancestorArray, iRegion, currentLevel);

            //Loop for all ancestors to get the conditional cross-covariance matrix
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
            {
                evaluate_cross_covariance(KCholTimesw[iRegion][jLevel], partition->knotsX[iRegion], partition->knotsY[iRegion], partition->knotsX[ancestorArray[jLevel]], partition->knotsY[ancestorArray[jLevel]], partition->nKnots, partition->nKnots, sill, range, nugget);

                get_conditional_covariance_matrix(jLevel, KCholTimesw[iRegion], KCholTimesw[ancestorArray[jLevel]],KCholTimesw[iRegion][jLevel],nKnots,nKnots,nKnots);
                
                dtrtrs_(&L,&N,&N,&nKnots,&nKnots,(RChol[ancestorArray[jLevel]]).memptr(),&nKnots,KCholTimesw[iRegion][jLevel],&nKnots,&status);
            }

            //Get the conditional variance-covariance matrix of this region
            evaluate_variance_covariance(RChol[iRegion].memptr(), partition->knotsX[iRegion], partition->knotsY[iRegion], partition->nKnots, sill, range, nugget);

            get_conditional_covariance_matrix(currentLevel, KCholTimesw[iRegion], KCholTimesw[iRegion], RChol[iRegion].memptr(),nKnots,nKnots,nKnots);

            //dpotrf_(&L,&nKnots,(RChol[iRegion]).memptr(),&nKnots,&status);  
            RChol[iRegion] = chol(RChol[iRegion],"lower");

        }

        timeval timeNow;
        gettimeofday(&timeNow, NULL);

        cout<<"Prior: Level "<<currentLevel+1<<" is complete. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n";
    }

    delete []ancestorArray; ancestorArray = NULL;
}


void Approximation::aggregate_A(mat &A, double *tempMemory, const unsigned long &indexA, const unsigned long &indexStart, const unsigned long &indexEnd, unsigned long iRegion)
{
    bool firstChild = true;

    for(unsigned long jChild = indexStart; jChild < indexEnd + 1; jChild++)
    {
        if(firstChild)
        {
            A = ATilde[jChild].slice(indexA);
            firstChild = false;
        }
        else
            A += ATilde[jChild].slice(indexA);
    }

}
void Approximation::aggregate_w_and_A(vec &w, mat &A, double *tempMemory, const int &jLevel, const unsigned long &indexA, const unsigned long &indexStart, const unsigned long &indexEnd, unsigned long iRegion)
{
    bool firstChild = true;
    for(unsigned long jChild = indexStart; jChild < indexEnd+1; jChild++)
    {
        if(firstChild)
        {
            w = wTilde[jChild].unsafe_col(jLevel);
            A = ATilde[jChild].slice(indexA);
            firstChild = false;
        }
        else
        {
            w += wTilde[jChild].unsafe_col(jLevel);
            A += ATilde[jChild].slice(indexA);
        }
    }
}

void Approximation::get_wTilde_and_ATilde_in_posterior(vec &w, mat &A, double *tempMemory, const unsigned long &nKnots, cube &KCholTimesCurrentA, const vec &KCholTimesCurrentw, const unsigned long &iRegion, const int &currentLevel)
{
    if(currentLevel == 0) return;

    wTilde[iRegion].set_size(nKnots,currentLevel);
    ATilde[iRegion].set_size(nKnots,nKnots,(currentLevel+1)*currentLevel/2);

    unsigned long indexAChild = (currentLevel+1)*(currentLevel+2)/2-1;
    unsigned long indexACurrent = currentLevel*(currentLevel+1)/2;
    for(long jLevel = currentLevel-1; jLevel > -1; jLevel--)
    {
        aggregate_w_and_A(w, A, tempMemory, jLevel, --indexAChild, iRegion*NUM_PARTITIONS_J+1, iRegion*NUM_PARTITIONS_J+NUM_PARTITIONS_J, iRegion);

        KCholTimesCurrentA.slice(jLevel) = solve(trimatl(RChol[iRegion]),A.t());
        wTilde[iRegion].unsafe_col(jLevel) = w-KCholTimesCurrentA.slice(jLevel).t()*KCholTimesCurrentw;

        for(long kLevelAfterjLevel = currentLevel-1; kLevelAfterjLevel > jLevel-1; kLevelAfterjLevel--)
        {
            aggregate_A(A, tempMemory, --indexAChild, iRegion*NUM_PARTITIONS_J+1, iRegion*NUM_PARTITIONS_J+NUM_PARTITIONS_J, iRegion);

            ATilde[iRegion].slice(--indexACurrent) = A-KCholTimesCurrentA.slice(jLevel).t()*KCholTimesCurrentA.slice(kLevelAfterjLevel);
        }
    }
}

void Approximation::free_wTilde_and_ATilde(const unsigned long &indexStartFreeingATilde, const unsigned long &indexEndFreeingATilde)
{
    for(unsigned long i = indexStartFreeingATilde; i < indexEndFreeingATilde + 1; i++)
    {
        ATilde[i].reset();
        wTilde[i].reset();
    }
}

void Approximation::load_ATilde_from_disk(const unsigned long &indexChildrenStart, const unsigned long &indexChildrenEnd, bool secondLastLevel, unsigned long* nKnotsAtFinestLevel, unsigned long offset, unsigned long size1, unsigned long size2)
{
    for(unsigned long i = indexChildrenStart; i < indexChildrenEnd + 1; i++)
    {
        if( secondLastLevel && nKnotsAtFinestLevel[i-offset]==0) 
        {
            ATilde[i].set_size(size1,size1,size2);
            ATilde[i].zeros();
            continue;
        }
        std::string fileName=TMP_DIRECTORY+"/"+std::to_string(i)+".bin";
        ATilde[i].load(fileName,arma_binary);
    }
}

/////////////////// Above are private member functions /////////////////

/////////////////// Below are public member functions /////////////////

Approximation::Approximation(const double &sill, const double &range, const double &nugget)
{
    (*this).sill=sill;
    (*this).range=range;
    (*this).nugget=nugget;
    loglikelihood = 0;
}

void Approximation::create_prior()
{
    gettimeofday(&timeNow, NULL);
    cout<<"===========================> Calculating prior quantities starts. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n";

    prior_allocate_memory();

    //Allocate memory for temporary variables storing K^{1/2}*W
    double ***KCholTimesw = new double** [partition->nRegionsInTotal];

    //Loop each region before the finest level
    loop_regions_before_finest_level_in_prior(KCholTimesw);
    
    //Loop each region in the finest level
    loop_regions_at_finest_level_in_prior(KCholTimesw);
    
    timeval timeNow;
    gettimeofday(&timeNow, NULL);

    cout<<"Prior: Level "<<NUM_LEVELS_M<<" is complete. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n";

    for(int currentLevel = 0; currentLevel < NUM_LEVELS_M-1; currentLevel++)
    {
        for(unsigned long iRegion = REGION_START[currentLevel]; iRegion < REGION_END[currentLevel] + 1; iRegion++)
        {
            for(int jLevel = 0; jLevel < currentLevel; jLevel++)
                delete[] KCholTimesw[iRegion][jLevel];
            delete[] KCholTimesw[iRegion];
        }
    }

    if(CALCULATION_MODE=="prediction")
    {
        for(unsigned long iRegion = REGION_START[NUM_LEVELS_M-1]; iRegion < REGION_END[NUM_LEVELS_M-1] + 1; iRegion++)
        {
            if(partition->nKnotsAtFinestLevel[iRegion-partition->nRegionsInTotal+partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]] == 0) continue;

            for(int jLevel = 0; jLevel < NUM_LEVELS_M - 1; jLevel++)
                delete[] KCholTimesw[iRegion][jLevel];
            delete[] KCholTimesw[iRegion];
        }
    }

    delete [] KCholTimesw; KCholTimesw=NULL;

    gettimeofday(&timeNow, NULL);
    cout<<"===========================> Calculating prior quantities is complete. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n\n";
}


void Approximation::likelihood()
{
	create_prior();
	posterior_inference();
}


void Approximation::dump_prediction_result()
{
    gettimeofday(&timeNow, NULL);
    cout<<"===========================> Dumping prediction results starts. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n";
				
    std::ofstream file;
    file.open(PREDICTION_RESULTS_FILE_NAME.c_str(), ios::binary);

    unsigned long offset = partition->nRegionsAtEachLevel[NUM_LEVELS_M-1]-partition->nRegionsInTotal;
    if(file.is_open())
    {
        //Write the number of prediction locations
        unsigned long n = 0;
        for(unsigned long indexRegionAtThisLevel = REGION_START[NUM_LEVELS_M-1]+offset; indexRegionAtThisLevel < REGION_END[NUM_LEVELS_M-1]+1+offset; indexRegionAtThisLevel++)
            n += partition->nPredictionsAtFinestLevel[indexRegionAtThisLevel];
        file.write((char*)&n,sizeof(unsigned long));
        
        //Write the prediction longtitude
        for(unsigned long indexRegionAtThisLevel = REGION_START[NUM_LEVELS_M-1]+offset; indexRegionAtThisLevel < REGION_END[NUM_LEVELS_M-1]+1+offset; indexRegionAtThisLevel++)
        {
            if(partition->nPredictionsAtFinestLevel[indexRegionAtThisLevel]==0) continue;
            file.write((char*)partition->predictionX[indexRegionAtThisLevel],sizeof(double)*partition->nPredictionsAtFinestLevel[indexRegionAtThisLevel]);
        }

        //Write the prediction lattitude
        for(unsigned long indexRegionAtThisLevel = REGION_START[NUM_LEVELS_M-1]+offset; indexRegionAtThisLevel < REGION_END[NUM_LEVELS_M-1]+1+offset; indexRegionAtThisLevel++)
        {
            if(partition->nPredictionsAtFinestLevel[indexRegionAtThisLevel]==0) continue;
            file.write((char*)partition->predictionY[indexRegionAtThisLevel],sizeof(double)*partition->nPredictionsAtFinestLevel[indexRegionAtThisLevel]);
        }

        //Write the prediction mean
        for(unsigned long indexRegionAtThisLevel = REGION_START[NUM_LEVELS_M-1]+offset; indexRegionAtThisLevel < REGION_END[NUM_LEVELS_M-1]+1+offset; indexRegionAtThisLevel++)
            for(unsigned long jPrediction = 0; jPrediction < partition->nPredictionsAtFinestLevel[indexRegionAtThisLevel]; jPrediction++)
            {
                vec designVec(3);
                designVec[1] = partition->predictionX[indexRegionAtThisLevel][jPrediction];
                designVec[2] = partition->predictionY[indexRegionAtThisLevel][jPrediction];
                designVec[0] = 1.0;
                double mean = posteriorPredictionMean[indexRegionAtThisLevel][jPrediction]+dot(designVec,data->coefficients);
                file.write((char*)&mean,sizeof(double));
            }
        
        //Write the prediction variance
        for(unsigned long indexRegionAtThisLevel = REGION_START[NUM_LEVELS_M-1]+offset; indexRegionAtThisLevel < REGION_END[NUM_LEVELS_M-1]+1+offset; indexRegionAtThisLevel++)
            for(unsigned long jPrediction = 0; jPrediction < partition->nPredictionsAtFinestLevel[indexRegionAtThisLevel]; jPrediction++)
                file.write((char*)&posteriorPredictionVariance[indexRegionAtThisLevel][jPrediction],sizeof(double));
        
        //Close the data file
        file.close();
    }else 
    {
        cout<<"Program exits with an error: fail to open the output file "<<PREDICTION_RESULTS_FILE_NAME<<endl;
        exit(EXIT_FAILURE);
    }

    gettimeofday(&timeNow, NULL);
    cout<<"===========================> Dumping prediction results is complete. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n\n";
}

Approximation::~Approximation()
{
    delete[] RChol;
        
    delete[] wTilde;
    delete[] ATilde;

    delete[] KCholTimesCurrentw;
    delete[] KCholTimesCurrentA;

    if(CALCULATION_MODE == "prediction")
    {
        delete[] BTilde;

        delete[] posteriorPredictionMean;
        delete[] posteriorPredictionVariance;
    }
    //cout<<"Approximation is deleted\n";
}



