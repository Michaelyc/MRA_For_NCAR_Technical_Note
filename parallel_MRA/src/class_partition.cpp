#include <iostream>
#include <fstream>
#include <armadillo>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <unordered_set>

#include "class_partition.hpp"
#include "constants.hpp"

/////////////////// Below are private member functions /////////////////

/*   
Creats knots in the current region
Input: 
    xMin, xMax, yMin, yMax: the boudary of the current region
    nKnotsX, nKnotsY: the number of knots in the current region in each direction
Output: 
    currentKnotsX: the array of the coordinates in longitude for each knot
    currentKnotsY: the array of the coordinates in latitude for each knot
*/
void Partition::create_knots(double *&currentKnotsX, double *&currentKnotsY, const double &xMin, const double &xMax, const int &nKnotsX, const double &yMin, const double &yMax, const int &nKnotsY)
{
    double offsetX = (xMax-xMin)*OFFSET;
    double offsetY = (yMax-yMin)*OFFSET;

    double xStart = xMin+offsetX;
    double xEnd = xMax-offsetX;
    double yStart = yMin+offsetY;
    double yEnd = yMax-offsetY;

    double xIncrement = ( nKnotsX==1 ? 0 : (xEnd-xStart)/(nKnotsX-1));
    double yIncrement = ( nKnotsY==1 ? 0 : (yEnd-yStart)/(nKnotsY-1));

    #pragma omp parallel for
    for(int iKnotX = 0; iKnotX < nKnotsX; iKnotX++)
        for(int jKnotY = 0; jKnotY < nKnotsY; jKnotY++)
        {
            currentKnotsX[iKnotX*nKnotsY+jKnotY] = xEnd-(nKnotsX-iKnotX-1)*xIncrement;//xStart+iKnotX*xIncrement;
            currentKnotsY[iKnotX*nKnotsY+jKnotY] = yEnd-(nKnotsY-jKnotY-1)*yIncrement;//yStart+jKnotY*yIncrement;
        }
}

/*   
Creats partitions in the current region
*/
void Partition::create_partitions(unsigned long region, double *partitionXMin, double *partitionXMax, double *partitionYMin, double *partitionYMax)
{
    double xMin=partitionXMin[region];
    double xMax=partitionXMax[region];
    double yMin=partitionYMin[region];
    double yMax=partitionYMax[region];

    if(NUM_PARTITIONS_J==2)
    {
        if((xMax-xMin)>=(yMax-yMin))
        {
            double xMid = (xMax+xMin)/2;
            partitionXMin[region*2+1]=xMin; partitionXMin[region*2+2]=xMid;
            partitionXMax[region*2+1]=xMid; partitionXMax[region*2+2]=xMax;

            partitionYMin[region*2+1]=yMin; partitionYMin[region*2+2]=yMin;
            partitionYMax[region*2+1]=yMax; partitionYMax[region*2+2]=yMax;
        }else
        {
            double yMid = (yMax+yMin)/2;
            partitionXMin[region*2+1]=xMin; partitionXMin[region*2+2]=xMin;
            partitionXMax[region*2+1]=xMax; partitionXMax[region*2+2]=xMax;

            partitionYMin[region*2+1]=yMin; partitionYMin[region*2+2]=yMid;
            partitionYMax[region*2+1]=yMid; partitionYMax[region*2+2]=yMax;
        }
    }else
    {
        double xMid = (xMax+xMin)/2;
        double yMid = (yMax+yMin)/2;

        partitionXMin[region*4+1]=xMin; partitionXMin[region*4+2]=xMid; partitionXMin[region*4+3]=xMin; partitionXMin[region*4+4]=xMid;

        partitionXMax[region*4+1]=xMid; partitionXMax[region*4+2]=xMax; partitionXMax[region*4+3]=xMid; partitionXMax[region*4+4]=xMax;

        partitionYMin[region*4+1]=yMin; partitionYMin[region*4+2]=yMin; partitionYMin[region*4+3]=yMid; partitionYMin[region*4+4]=yMid;

        partitionYMax[region*4+1]=yMid; partitionYMax[region*4+2]=yMid; partitionYMax[region*4+3]=yMax; partitionYMax[region*4+4]=yMax;
    }
}

/////////////////// Above are private member functions /////////////////

/////////////////// Below are public member functions /////////////////

/*Implement the constructor of Partition. Assign values to nRegionsAtEachLevel for the array of the number of regions in each level, nRegionsInTotal for the total number of regions, and MPI quantities.
*/
Partition::Partition()
{
    //If user did not predetermine the number of levels, find the number of levels such that the average number of observations per region similar to the number of knots
    if(NUM_LEVELS_M == -99)//-99 stands for the default method to automatically determine NUM_LEVELS_M 
    {
        //Find the number of regions required at finest level, nRegionsAtFinestLevel 
        unsigned long nRegionsAtFinestLevel = ceil((double)NUM_OBSERVATIONS/(double)NUM_KNOTS_r);

        //Find the number of levels, NUM_LEVELS_M, which satisfies NUM_PARTITIONS_J^(NUM_LEVELS_M-1)>=nRegionsAtFinestLevel
        NUM_LEVELS_M = ceil(log(nRegionsAtFinestLevel)/log(NUM_PARTITIONS_J))+1;

        if(NUM_LEVELS_M < 2)
	    {
            if(WORKER == 0) cout<<"Program exits with an error: the NUM_LEVELS_M calculated by default is "<<NUM_LEVELS_M<<", which is required to be larger than 1. Please decrease NUM_KNOTS_r to get a larger NUM_LEVELS_M.\n";
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(EXIT_FAILURE);
	    }
    }

    //Find the array of the number of regions in each level, nRgnsVec, and the total number of regions, nRegionsInTotal
    nRegionsAtEachLevel = new unsigned long [NUM_LEVELS_M];
    nRegionsAtEachLevel[0] = 1;
    nRegionsInTotal = 1;
    for(int iLevel = 1; iLevel < NUM_LEVELS_M; iLevel++)
    {
        nRegionsAtEachLevel[iLevel] = nRegionsAtEachLevel[iLevel-1]*NUM_PARTITIONS_J;
        nRegionsInTotal += nRegionsAtEachLevel[iLevel];
    }

    //Calculate MPI information
        WORKING_REGION_FLAG = new bool [nRegionsInTotal]();

        unsigned long indexRegionAtFinestLevelStartThisWorker, indexRegionAtFinestLevelEndThisWorker;
        WORKERS_FOR_EACH_REGION = new std::set<unsigned long>[nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1]];

        //Assign regions in the finest level to the worker
        unsigned long nRegionsPerWorker = nRegionsAtEachLevel[NUM_LEVELS_M-1] / MPI_SIZE;
        unsigned long nWorkersWithAdditionalRegion = nRegionsAtEachLevel[NUM_LEVELS_M-1] - nRegionsPerWorker * MPI_SIZE;


        if(WORKER < nWorkersWithAdditionalRegion)
        {
            indexRegionAtFinestLevelStartThisWorker =  WORKER * (nRegionsPerWorker+1);
            indexRegionAtFinestLevelEndThisWorker =  indexRegionAtFinestLevelStartThisWorker + nRegionsPerWorker;
        }
        else
        {
            indexRegionAtFinestLevelStartThisWorker =  nWorkersWithAdditionalRegion * (nRegionsPerWorker+1) + (WORKER -  nWorkersWithAdditionalRegion)* nRegionsPerWorker;
            indexRegionAtFinestLevelEndThisWorker =  indexRegionAtFinestLevelStartThisWorker + nRegionsPerWorker - 1;
        }

        indexRegionAtFinestLevelStartThisWorker += nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1];
        indexRegionAtFinestLevelEndThisWorker += nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1];

        for(unsigned long i = indexRegionAtFinestLevelStartThisWorker; i < indexRegionAtFinestLevelEndThisWorker+1; i++)
            WORKING_REGION_FLAG[i] = true;

        //Find the regions dealt with by this worker at the second finest level
        for(unsigned long index = indexRegionAtFinestLevelStartThisWorker; index <= indexRegionAtFinestLevelEndThisWorker; index++)
        {
            unsigned long ancestor = index;
            while(ancestor > 0)
            {
                ancestor = (ancestor-1)/NUM_PARTITIONS_J;
                //INDICES_REGIONS.insert(ancestor);
                WORKING_REGION_FLAG[ancestor] = true;
            }
            INDICES_REGIONS_AT_CURRENT_LEVEL.insert((index-1)/NUM_PARTITIONS_J);
        }

        //Assign WORKERS_FOR_EACH_REGION at the second finest level
        unsigned long indexStartFinestLevel = nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1];
        unsigned long tmpIndexStart = indexStartFinestLevel, tmpIndexEnd;
        if(DYNAMIC_SCHEDULE_FLAG) origninalWorker = new int [nRegionsInTotal];

        for(int iWorker = 0; iWorker < nWorkersWithAdditionalRegion; iWorker++)
        {
            tmpIndexEnd =  tmpIndexStart + nRegionsPerWorker + 1;
            for(unsigned long index = tmpIndexStart; index < tmpIndexEnd; index++)
            {
                unsigned long ancestor = (index-1)/NUM_PARTITIONS_J;
                WORKERS_FOR_EACH_REGION[ancestor].insert(iWorker);
                if(DYNAMIC_SCHEDULE_FLAG) 
                {
                    origninalWorker[index] = iWorker;
                    origninalWorker[ancestor] = iWorker;
                    while(ancestor > 0)
                    {
                        ancestor = (ancestor-1)/NUM_PARTITIONS_J;
                        origninalWorker[ancestor] = iWorker;
                    }
                }
            }
            tmpIndexStart += nRegionsPerWorker+1;
            tmpIndexEnd += nRegionsPerWorker+1;
        }
        for(int iWorker = nWorkersWithAdditionalRegion; iWorker < MPI_SIZE; iWorker++)
        {
            tmpIndexEnd =  tmpIndexStart + nRegionsPerWorker;
            for(unsigned long index = tmpIndexStart; index < tmpIndexEnd; index++)
            {
                unsigned long ancestor = (index-1)/NUM_PARTITIONS_J;
                WORKERS_FOR_EACH_REGION[ancestor].insert(iWorker);
                if(DYNAMIC_SCHEDULE_FLAG) 
                {
                    origninalWorker[index] = iWorker;
                    origninalWorker[ancestor] = iWorker;
                    while(ancestor > 0)
                    {
                        ancestor = (ancestor-1)/NUM_PARTITIONS_J;
                        origninalWorker[ancestor] = iWorker;
                    }
                }
            }
            tmpIndexStart += nRegionsPerWorker;
            tmpIndexEnd += nRegionsPerWorker;
        }

        //Assign regionStart and regionEnd for each level
        REGION_START = new unsigned long [NUM_LEVELS_M];
        REGION_END = new unsigned long [NUM_LEVELS_M];
        REGION_START[NUM_LEVELS_M-1]=indexRegionAtFinestLevelStartThisWorker;
        REGION_END[NUM_LEVELS_M-1]=indexRegionAtFinestLevelEndThisWorker;
        for(int iLevel = NUM_LEVELS_M-2; iLevel > -1; iLevel--)
        {
            REGION_START[iLevel] = (REGION_START[iLevel+1]-1)/NUM_PARTITIONS_J;
            REGION_END[iLevel] = (REGION_END[iLevel+1]-1)/NUM_PARTITIONS_J;
        }
}

/*
Implement the member function of Partition, build_partition, which builds the hierarchical partition
*/

void Partition::build_partition()
{
    if(WORKER == 0)
	{
		gettimeofday(&timeNow, NULL);
		cout<<"===========================> Processor 1: building hierarchical grid starts. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n";
	}

    //The number of knots in each direction
    int nKnotsX = ceil(sqrt(NUM_KNOTS_r));
    int nKnotsY = (int)(NUM_KNOTS_r/nKnotsX);
    nKnots = nKnotsX*nKnotsY;

    //Allocate memory for the coordinates of knots
    knotsX = new double* [nRegionsInTotal];
    knotsY = new double* [nRegionsInTotal];

    //Array of coordinates for the partition in each region with dimension [nRegionsInTotal]
    double *partitionXMin, *partitionYMin, *partitionXMax, *partitionYMax;
    partitionXMin = new double [nRegionsInTotal];
    partitionXMax = new double [nRegionsInTotal];
    partitionYMin = new double [nRegionsInTotal];
    partitionYMax = new double [nRegionsInTotal];

    partitionXMin[0] = data->domainBoundaries[0];
    partitionXMax[0] = data->domainBoundaries[1];
    partitionYMin[0] = data->domainBoundaries[2];
    partitionYMax[0] = data->domainBoundaries[3];

    for(int iLevel = 0; iLevel < NUM_LEVELS_M-1; iLevel++)
    {
        #pragma omp parallel for
        for(unsigned long jRegion = REGION_START[iLevel]; jRegion < REGION_END[iLevel] + 1; jRegion++)
            create_partitions(jRegion, partitionXMin, partitionXMax, partitionYMin, partitionYMax);
    }

    #pragma omp parallel for schedule(dynamic,1)
    for(unsigned long iRegion = 0; iRegion < nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion++)
    {
        if(WORKING_REGION_FLAG[iRegion])
        {
            knotsX[iRegion] = new double [nKnots];
            knotsY[iRegion] = new double [nKnots];
            create_knots(knotsX[iRegion], knotsY[iRegion], partitionXMin[iRegion], partitionXMax[iRegion], nKnotsX, partitionYMin[iRegion], partitionYMax[iRegion], nKnotsY);
        }
    }

    //Allocate memory for nKnotsAtFinestLevel, nPredictionsAtFinestLevel
    nKnotsAtFinestLevel = new unsigned long [nRegionsAtEachLevel[NUM_LEVELS_M-1]];
    nPredictionsAtFinestLevel = new unsigned long [nRegionsAtEachLevel[NUM_LEVELS_M-1]];

    //Allocate memory for nKnotsAtFinestLevel, nPredictionsAtFinestLevel
    knotsResidual = new double* [nRegionsAtEachLevel[NUM_LEVELS_M-1]];
    
    //Allocate memory for the coordinates of predictions
    predictionX = new double* [nRegionsAtEachLevel[NUM_LEVELS_M-1]];
    predictionY = new double* [nRegionsAtEachLevel[NUM_LEVELS_M-1]];
    
    //Get coordinates of the knots and predictions at the regions in the finest level
    int maxOpenMPThreads=omp_get_max_threads();

    #pragma omp parallel for num_threads(maxOpenMPThreads) schedule(dynamic,1)
    for(unsigned long iRegion = nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion < nRegionsInTotal; iRegion++)
    {
        std::set<unsigned long> *indexObservationsInThisRegion = new std::set<unsigned long>;
        std::set<unsigned long> *indexPredictionsInThisRegion = new std::set<unsigned long>;

        unsigned long indexStartThisLevel = iRegion+nRegionsAtEachLevel[NUM_LEVELS_M-1]-nRegionsInTotal;
        if(WORKING_REGION_FLAG[iRegion])
        {
            double xMin=partitionXMin[iRegion];
            double xMax=partitionXMax[iRegion];
            double yMin=partitionYMin[iRegion];
            double yMax=partitionYMax[iRegion];

            indexObservationsInThisRegion->clear();

            for(unsigned long jObservation = 0; jObservation < NUM_OBSERVATIONS; jObservation++)
            {
                if(data->observationLon[jObservation] >= xMin && data->observationLon[jObservation] < xMax && data->observationLat[jObservation] >= yMin && data->observationLat[jObservation] < yMax)
                {
                    
                    bool skip = false;
                    //Skip if the observation location coincides with knots in levels above
                    unsigned long ancestor = iRegion;
                    for(int kLevel = 0; kLevel < NUM_LEVELS_M-1; kLevel++)
                    {
                        ancestor = (ancestor-1)/NUM_PARTITIONS_J;
                        for(int pKnot = 0; pKnot < nKnots; pKnot++)
                        {
                            if(knotsX[ancestor][pKnot] == data->observationLon[jObservation] && knotsY[ancestor][pKnot] == data->observationLat[jObservation])
                            {
                                skip = true;
                                break;
                            }
                        }
                        if(skip) break;
                    }
                    
                    //Add this observation
                    if(!skip) indexObservationsInThisRegion->insert(jObservation);
                }
            }

            nKnotsAtFinestLevel[indexStartThisLevel] = indexObservationsInThisRegion->size();
            knotsX[iRegion] = new double [nKnotsAtFinestLevel[indexStartThisLevel]];            
            knotsY[iRegion] = new double [nKnotsAtFinestLevel[indexStartThisLevel]];
            knotsResidual[indexStartThisLevel] = new double [nKnotsAtFinestLevel[indexStartThisLevel]];
       
            unsigned long tmpIndex=0;
            for(std::set<unsigned long>::iterator jObservation = indexObservationsInThisRegion->begin(); jObservation != indexObservationsInThisRegion->end(); jObservation++)
            {
                knotsX[iRegion][tmpIndex] = data->observationLon[*jObservation];
                knotsY[iRegion][tmpIndex] = data->observationLat[*jObservation];
                knotsResidual[indexStartThisLevel][tmpIndex++] = data->observationResiduals[*jObservation];
            }

            if(CALCULATION_MODE=="prediction")
            {
                indexPredictionsInThisRegion->clear();
                for(unsigned long jPrediction = 0; jPrediction < NUM_PREDICTIONS; jPrediction++)
                {
                    if(data->predictionLon[jPrediction] >= xMin && data->predictionLon[jPrediction] < xMax && data->predictionLat[jPrediction] >= yMin && data->predictionLat[jPrediction] < yMax)
                        indexPredictionsInThisRegion->insert(jPrediction);
                }

                nPredictionsAtFinestLevel[indexStartThisLevel] = indexPredictionsInThisRegion->size();
                predictionX[indexStartThisLevel] = new double [nPredictionsAtFinestLevel[indexStartThisLevel]];            
                predictionY[indexStartThisLevel] = new double [nPredictionsAtFinestLevel[indexStartThisLevel]];

                unsigned long tmpIndex=0;
                for(std::set<unsigned long>::iterator jPrediction = indexPredictionsInThisRegion->begin(); jPrediction != indexPredictionsInThisRegion->end(); jPrediction++)
                {
                    predictionX[indexStartThisLevel][tmpIndex] = data->predictionLon[*jPrediction];
                    predictionY[indexStartThisLevel][tmpIndex++] = data->predictionLat[*jPrediction];
                }
            }
        }    
        indexObservationsInThisRegion->clear();
        indexPredictionsInThisRegion->clear();
    }
    

    //Release the memory for the temporary arrays
    delete[] partitionXMin;
    delete[] partitionYMin;
    delete[] partitionXMax;
    delete[] partitionYMax;

    //Delete variables that are not used anymore
    delete[] data->observationLon;
    delete[] data->observationLat;
    delete[] data->observationResiduals;

    if(DYNAMIC_SCHEDULE_FLAG) 
		dynamic_schedule();
	else
		WORLD = MPI_COMM_WORLD;

    if(WORKER == 0)
	{
		gettimeofday(&timeNow, NULL);
		cout<<"===========================> Processor 1: building hierarchical grid is complete. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n\n";
		if(PRINT_DETAIL_FLAG) print_partition_summary();
	}
}

//Implement the member function of Partition, dynamic_schedule, which assign the work load for each worker by dynamic scheduling
void Partition::dynamic_schedule()
{
    //Synchronize the knots of regions before the finest level
    for(unsigned long iRegion = 0; iRegion < nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion++)
    {
        if(!WORKING_REGION_FLAG[iRegion])
        {
            knotsX[iRegion] = new double [nKnots];
            knotsY[iRegion] = new double [nKnots];
        }
    }

    for(unsigned long iRegion = 0; iRegion < nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion++)
    {
        MPI_Bcast(knotsX[iRegion],nKnots,MPI_DOUBLE,origninalWorker[iRegion],MPI_COMM_WORLD);
        MPI_Bcast(knotsY[iRegion],nKnots,MPI_DOUBLE,origninalWorker[iRegion],MPI_COMM_WORLD);
    }

    //Synchronize the knots of regions at the finest level and predictions if predicting
    unsigned long indexStartFinestLevel = nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1];
    for(unsigned long iRegion = indexStartFinestLevel; iRegion < nRegionsInTotal; iRegion++)
    {

        MPI_Bcast(&nKnotsAtFinestLevel[iRegion-indexStartFinestLevel],1,MPI_UNSIGNED_LONG,origninalWorker[iRegion],MPI_COMM_WORLD);
        if(!WORKING_REGION_FLAG[iRegion])
        {
            knotsX[iRegion] = new double [nKnotsAtFinestLevel[iRegion-indexStartFinestLevel]];
            knotsY[iRegion] = new double [nKnotsAtFinestLevel[iRegion-indexStartFinestLevel]];
            knotsResidual[iRegion-indexStartFinestLevel] = new double [nKnotsAtFinestLevel[iRegion-indexStartFinestLevel]];
        }

        if(CALCULATION_MODE=="prediction")
        {
            MPI_Bcast(&nPredictionsAtFinestLevel[iRegion-indexStartFinestLevel],1,MPI_UNSIGNED_LONG,origninalWorker[iRegion],MPI_COMM_WORLD);
            if(!WORKING_REGION_FLAG[iRegion])
            {
                predictionX[iRegion-indexStartFinestLevel] = new double [nPredictionsAtFinestLevel[iRegion-indexStartFinestLevel]];
                predictionY[iRegion-indexStartFinestLevel] = new double [nPredictionsAtFinestLevel[iRegion-indexStartFinestLevel]];
            }
        }

        MPI_Bcast(knotsX[iRegion],nKnotsAtFinestLevel[iRegion-indexStartFinestLevel],MPI_DOUBLE,origninalWorker[iRegion],MPI_COMM_WORLD);
        MPI_Bcast(knotsY[iRegion],nKnotsAtFinestLevel[iRegion-indexStartFinestLevel],MPI_DOUBLE,origninalWorker[iRegion],MPI_COMM_WORLD);
        MPI_Bcast(knotsResidual[iRegion-indexStartFinestLevel],nKnotsAtFinestLevel[iRegion-indexStartFinestLevel],MPI_DOUBLE,origninalWorker[iRegion],MPI_COMM_WORLD);

        if(CALCULATION_MODE=="prediction")
        {
            MPI_Bcast(predictionX[iRegion-indexStartFinestLevel],nPredictionsAtFinestLevel[iRegion-indexStartFinestLevel],MPI_DOUBLE,origninalWorker[iRegion],MPI_COMM_WORLD);
            MPI_Bcast(predictionY[iRegion-indexStartFinestLevel],nPredictionsAtFinestLevel[iRegion-indexStartFinestLevel],MPI_DOUBLE,origninalWorker[iRegion],MPI_COMM_WORLD);
        }
    }

    //Dynamic schedule work load for each worker
    
    ///Clear WORKERS_FOR_EACH_REGION at the second finest level
    for(unsigned long iRegion = nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1] - nRegionsAtEachLevel[NUM_LEVELS_M-2]; iRegion < nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion++)
        WORKERS_FOR_EACH_REGION[iRegion].clear();

    ///Clear WORKING_REGION_FLAG for all regions
    memset(WORKING_REGION_FLAG,false,nRegionsInTotal);
    
    ///Clear INDICES_REGIONS_AT_CURRENT_LEVEL
    INDICES_REGIONS_AT_CURRENT_LEVEL.clear();

    ///Calculate work load
    double totalComplexity = 0;
    for(unsigned long index = 0; index < nRegionsAtEachLevel[NUM_LEVELS_M-1]; index++)
    {
        double tmp = nKnotsAtFinestLevel[index];
        totalComplexity += tmp*tmp;//*tmp;
    }
    double complexityEachWorker = totalComplexity/MPI_SIZE;

    ///Assign working regions for each worker
    REGION_START[NUM_LEVELS_M-1]=nRegionsInTotal;
    REGION_END[NUM_LEVELS_M-1]=0;
    int thisWorker=-1;
    double complexityThisWorker=0;  
    unsigned long indexRegionAtFinestLevelStartThisWorker = nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1], indexRegionAtFinestLevelEndThisWorker;
    for(unsigned long index = 0; index < nRegionsAtEachLevel[NUM_LEVELS_M-1]; index++)
    {
        double tmp = nKnotsAtFinestLevel[index];
        double complexityThisRegion = tmp*tmp;//*tmp;
        complexityThisWorker += complexityThisRegion;
        //if(WORKER==0) cout<<complexityThisRegion<<" "<<complexityEachWorker<<endl;
        if(complexityThisWorker >= complexityEachWorker || index == nRegionsAtEachLevel[NUM_LEVELS_M-1]-1)
        {
            thisWorker++;
            if(thisWorker == MPI_SIZE) thisWorker = MPI_SIZE-1;

            indexRegionAtFinestLevelEndThisWorker = indexStartFinestLevel+index;
            
            //Check whether this region should be handled by this worker or next worker
            if(complexityEachWorker - complexityThisWorker + complexityThisRegion < complexityThisWorker - complexityEachWorker 
                && index != nRegionsAtEachLevel[NUM_LEVELS_M-1]-1
                && indexRegionAtFinestLevelEndThisWorker > indexRegionAtFinestLevelStartThisWorker)
            {
                //Next worker handle this region
                indexRegionAtFinestLevelEndThisWorker--;
                complexityThisWorker = complexityThisRegion;
            }
            else
            {
                //This worker handle this region
                complexityThisWorker = 0;
            }

            //if(WORKER==0) cout<<thisWorker<<" "<<indexRegionAtFinestLevelStartThisWorker<<" "<<indexRegionAtFinestLevelEndThisWorker<<endl;

            for(unsigned long jRegion = indexRegionAtFinestLevelStartThisWorker; jRegion < indexRegionAtFinestLevelEndThisWorker+1; jRegion++)
            {
                unsigned long ancestor = (jRegion-1)/NUM_PARTITIONS_J;
                WORKERS_FOR_EACH_REGION[ancestor].insert(thisWorker);
            }

            if(thisWorker == WORKER)
            {
                for(unsigned long jRegion = indexRegionAtFinestLevelStartThisWorker; jRegion < indexRegionAtFinestLevelEndThisWorker+1; jRegion++)
                {
                    unsigned long ancestor = jRegion;
                    WORKING_REGION_FLAG[ancestor] = true;
                    while(ancestor > 0)
                    {
                        ancestor = (ancestor-1)/NUM_PARTITIONS_J;
                        WORKING_REGION_FLAG[ancestor] = true;
                    }
                    INDICES_REGIONS_AT_CURRENT_LEVEL.insert((jRegion-1)/NUM_PARTITIONS_J);
                }
                if(indexRegionAtFinestLevelStartThisWorker<REGION_START[NUM_LEVELS_M-1]) REGION_START[NUM_LEVELS_M-1]=indexRegionAtFinestLevelStartThisWorker;
                if(indexRegionAtFinestLevelEndThisWorker>REGION_END[NUM_LEVELS_M-1]) REGION_END[NUM_LEVELS_M-1]=indexRegionAtFinestLevelEndThisWorker;
            }

            indexRegionAtFinestLevelStartThisWorker = indexRegionAtFinestLevelEndThisWorker+1;
        }
    }

    //Deactivate MPI processes that are not assigned any regions, which rarely happens
    if(thisWorker<MPI_SIZE-1) 
    {
        MPI_Group worldGroup;
        MPI_Comm_group(MPI_COMM_WORLD,&worldGroup);

        MPI_Group newWorldGroup;
        int ranges[1][3]={{thisWorker+1,MPI_SIZE-1,1}};
        MPI_Group_range_excl(worldGroup,1,ranges,&newWorldGroup);

        MPI_Comm_create(MPI_COMM_WORLD, newWorldGroup, &WORLD);

        MPI_SIZE=thisWorker+1;

        if(WORKER>thisWorker)
        {
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }
    else
        WORLD = MPI_COMM_WORLD;
    
    for(int iLevel = NUM_LEVELS_M-2; iLevel > -1; iLevel--)
    {
        REGION_START[iLevel] = (REGION_START[iLevel+1]-1)/NUM_PARTITIONS_J;
        REGION_END[iLevel] = (REGION_END[iLevel+1]-1)/NUM_PARTITIONS_J;
    }

    //Delete regions that are not handled by this WORKER
    for(unsigned long iRegion = 0; iRegion < nRegionsInTotal; iRegion++)
        if(!WORKING_REGION_FLAG[iRegion])
        {
            delete[] knotsX[iRegion]; knotsX[iRegion]=NULL;
            delete[] knotsY[iRegion]; knotsY[iRegion]=NULL;
        }

    for(unsigned long iRegion = nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion < nRegionsInTotal; iRegion++)
        if(!WORKING_REGION_FLAG[iRegion])
        {
            delete[] knotsResidual[iRegion-indexStartFinestLevel]; knotsResidual[iRegion-indexStartFinestLevel]=NULL;

            if(CALCULATION_MODE=="prediction")
            {
                delete[] predictionX[iRegion-indexStartFinestLevel]; predictionX[iRegion-indexStartFinestLevel]=NULL;
                delete[] predictionY[iRegion-indexStartFinestLevel]; predictionY[iRegion-indexStartFinestLevel]=NULL;
            }
        }

}

//Implement the member function of Partition, print_partition_summary, which shows a brief summary of the partition
void Partition::print_partition_summary()
{
    cout<<">>The number of levels: "<<NUM_LEVELS_M<<endl<<endl;
    cout<<">>The total number of regions: "<<nRegionsInTotal<<endl<<endl;
    cout<<">>The number of regions in each level: \n    ";
    for(int iLevel = 0; iLevel < NUM_LEVELS_M; iLevel++) cout<<nRegionsAtEachLevel[iLevel]<<" ";
    cout<<endl<<endl;

    cout<<">>Knot coordinates in the coarsest region, up to ten knots are shown: \n    longitude: ";
    for(int iKnot = 0; iKnot < std::min(nKnots,10); iKnot++) 
        printf("%8.2lf ",knotsX[0][iKnot]);
    cout<<"\n    latitude:  ";
    for(int iKnot = 0; iKnot < std::min(nKnots,10); iKnot++) 
        printf("%8.2lf ",knotsY[0][iKnot]);
    cout<<endl<<endl;

    if(CALCULATION_MODE=="predict")
        cout<<">>The number of predictions after eliminating duplicates: "<<NUM_PREDICTIONS<<endl<<endl;

}

void Partition::dump_structure_information()
{
    //Open the file "structure_information.txt" to store structure information
    std::ofstream file;
    file.open("structure_information.txt");

    //Dump the number of levels 
    file<<"The number of levels: "<<NUM_LEVELS_M<<endl<<endl;

    //Dump the total number of regions
    file<<"The total number of regions: "<<nRegionsInTotal<<endl<<endl;


    //Dump the number of knots in regions before the finest level
    file<<"The actual number of knots in regions before the finest level: "<<nKnots<<endl<<endl;

    //Dump some statistics of the number of knots in each region at the finest level
    unsigned long nKnotsMin = nKnotsAtFinestLevel[0];
    unsigned long nKnotsMax = nKnotsAtFinestLevel[0];
    unsigned long nRegionsHaveZeroKnots = 0;

    for(unsigned long iRegion = 0; iRegion < nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion++)
    {
        if(nKnotsMin>nKnotsAtFinestLevel[iRegion]) nKnotsMin=nKnotsAtFinestLevel[iRegion];
        if(nKnotsMax<nKnotsAtFinestLevel[iRegion]) nKnotsMax=nKnotsAtFinestLevel[iRegion];
        if(nKnotsAtFinestLevel[iRegion] == 0) nRegionsHaveZeroKnots++;
    }

    file<<"The minimal number of knots in regions at the finest level: "<<nKnotsMin<<endl<<endl;
    file<<"The maximal number of knots in regions at the finest level: "<<nKnotsMax<<endl<<endl;
    file<<"The number of in regions at the finest level that have zero knots: "<<nRegionsHaveZeroKnots<<endl<<endl;

    unsigned long thresholdLeft, thresholdRight;
    unsigned long nRegionsBetweenThresholds;
    double interval = (nKnotsMax-nKnotsMin)/10;
    
    thresholdLeft=nKnotsMin;
    file<<"The number of regions at the finest level with number of knots in the following intervals:"<<endl;
    
    for(int i = 0; i < 10; i++)
    {
        nRegionsBetweenThresholds = 0;
        thresholdRight=thresholdLeft+interval;
        if(i == 9) thresholdRight = nKnotsMax;
        for(unsigned long iRegion = 0; iRegion < nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion++)
           if(thresholdLeft<=nKnotsAtFinestLevel[iRegion] && nKnotsAtFinestLevel[iRegion]<=thresholdRight) nRegionsBetweenThresholds++;
        file<<"[ "<<thresholdLeft<<" , "<<thresholdRight<<" ]: "<<nRegionsBetweenThresholds<<endl;
        thresholdLeft = thresholdRight+1;
    }    
    file<<endl;

    //Dump the number of regions in each level
    file<<"The number of regions in each level: \n";
    for(int iLevel = 0; iLevel < NUM_LEVELS_M; iLevel++) 
        file<<"Level "<<iLevel+1<<": "<<nRegionsAtEachLevel[iLevel]<<endl;
    file<<endl;

    //Dump the number of knots in each region at the finest level
    file<<"The number of knots in each region at the finest level: \n";
    for(unsigned long iRegion = 0; iRegion < nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion++)
        file<<"Region "<<nRegionsInTotal-nRegionsAtEachLevel[NUM_LEVELS_M-1]+iRegion<<": "<<nKnotsAtFinestLevel[iRegion]<<endl;


    //Close the file
    file.close();
}

Partition::~Partition()
{
    delete[] nRegionsAtEachLevel;

    delete[] knotsX;
    delete[] knotsY;
    delete[] knotsResidual;

    delete[] predictionX, 
    delete[] predictionY;

    delete[] nKnotsAtFinestLevel, 
    delete[] nPredictionsAtFinestLevel;
    
    if(DYNAMIC_SCHEDULE_FLAG) delete[] origninalWorker;
    //cout<<"Partition is deleted\n";
}