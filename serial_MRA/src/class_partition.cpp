#include <iostream>
#include <fstream>
#include <armadillo>
#include <math.h>
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
            cout<<"Program exits with an error: the NUM_LEVELS_M calculated by default is "<<NUM_LEVELS_M<<", which is required to be larger than 1. Please decrease NUM_KNOTS_r to get a larger NUM_LEVELS_M.\n";
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

    //Assign regionStart and regionEnd for each level
    REGION_START = new unsigned long [NUM_LEVELS_M];
    REGION_END = new unsigned long [NUM_LEVELS_M];
    REGION_START[NUM_LEVELS_M-1] = nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1];
    REGION_END[NUM_LEVELS_M-1] = nRegionsInTotal - 1;
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
    gettimeofday(&timeNow, NULL);
    cout<<"===========================> Building hierarchical grid starts. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n";

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
        for(unsigned long jRegion = REGION_START[iLevel]; jRegion < REGION_END[iLevel] + 1; jRegion++)
            create_partitions(jRegion, partitionXMin, partitionXMax, partitionYMin, partitionYMax);
    }

    for(unsigned long iRegion = 0; iRegion < nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion++)
    {
        knotsX[iRegion] = new double [nKnots];
        knotsY[iRegion] = new double [nKnots];
        create_knots(knotsX[iRegion], knotsY[iRegion], partitionXMin[iRegion], partitionXMax[iRegion], nKnotsX, partitionYMin[iRegion], partitionYMax[iRegion], nKnotsY);
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
    for(unsigned long iRegion = nRegionsInTotal - nRegionsAtEachLevel[NUM_LEVELS_M-1]; iRegion < nRegionsInTotal; iRegion++)
    {
        std::set<unsigned long> *indexObservationsInThisRegion = new std::set<unsigned long>;
        std::set<unsigned long> *indexPredictionsInThisRegion = new std::set<unsigned long>;

        unsigned long indexStartThisLevel = iRegion+nRegionsAtEachLevel[NUM_LEVELS_M-1]-nRegionsInTotal;

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

    gettimeofday(&timeNow, NULL);
    cout<<"===========================> Building hierarchical grid is complete. Elapsed time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n\n";
    if(PRINT_DETAIL_FLAG) print_partition_summary();
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