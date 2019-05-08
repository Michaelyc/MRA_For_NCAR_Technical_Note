#include <iostream>
#include <mpi.h>

#include "dlib/optimization.h"
#include "dlib/global_optimization.h"

#include "class_data.hpp"
#include "class_partition.hpp"
#include "class_approximation.hpp"
#include "constants.hpp"
#include "optimization.hpp"

double loglikelihood=0;

double objective(const dlib::matrix<double,0,1>& parameter)
{
	double alpha=parameter(0);
	double beta=parameter(1);
	double tau=parameter(2);

	OPTIMIZATION_ITERATION++;
	if(WORKER==0)
	{
		cout<<"##########################################################################\n";
		cout<<"############## Iteration "<<OPTIMIZATION_ITERATION<<". alpha: "<<alpha<<", beta: "<<beta<<", tau: "<<tau<<endl;
		cout<<"##########################################################################\n\n";
	}

	INDICES_REGIONS_AT_CURRENT_LEVEL.clear();
	for(unsigned long iRegion = REGION_START[NUM_LEVELS_M-1]; iRegion < REGION_END[NUM_LEVELS_M-1]+1; iRegion++)
	{
		unsigned long ancestor = iRegion;
		WORKING_REGION_FLAG[ancestor] = true;
		while(ancestor > 0)
		{
			ancestor = (ancestor-1)/NUM_PARTITIONS_J;
			WORKING_REGION_FLAG[ancestor] = true;
		}
		INDICES_REGIONS_AT_CURRENT_LEVEL.insert((iRegion-1)/NUM_PARTITIONS_J);
	}

	MPI_Barrier(WORLD);
	Approximation approximation(alpha,beta,tau);
	approximation.likelihood();

	MPI_Bcast(&approximation.loglikelihood,1,MPI_DOUBLE,0,WORLD);
	MPI_Barrier(WORLD);

	//printf("\nIter %d, WORKER %d: %16.16lf %16.16lf %16.16lf %16.16lf\n",OPTIMIZATION_ITERATION+1, WORKER,alpha,beta,tau, approximation.loglikelihood);

	loglikelihood=-0.5*(approximation.loglikelihood+NUM_OBSERVATIONS*log(2*3.14159265359));

	return approximation.loglikelihood;
}


void get_optimal_parameters()
{
	OPTIMIZATION_ITERATION = 0;
	
	dlib::matrix<double,0,1> initialGuess={ALPHA_INITIAL_GUESS,BETA_INITIAL_GUESS,TAU_INITIAL_GUESS};
	dlib::matrix<double,0,1> lowerBound = {ALPHA_LOWER_BOUND,BETA_LOWER_BOUND,TAU_LOWER_BOUND};
	dlib::matrix<double,0,1> upperBound={ALPHA_UPPER_BOUND,BETA_UPPER_BOUND,TAU_UPPER_BOUND};

	double rhoBegin = 0.49*min(upperBound-lowerBound);


	bool optimizationSuccess=false;
	try
	{
		dlib::find_min_bobyqa(objective,initialGuess,7,lowerBound,upperBound,0.0001,1e-8,MAX_ITERATIONS);
		optimizationSuccess = true;

	}
	catch (dlib::bobyqa_failure error)
	{
		if(WORKER==0)
		{
			cout<<"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n";
			cout<<"\\\\\\\\\\\\\\The optimal parameter values have not been found. Reason: "<<error.what()<<endl;
		}
	}

	if(optimizationSuccess) 
	{
		if(WORKER==0)
		{
			cout<<"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n";
			cout<<"\\\\\\\\\\\\\\The optimal parameter values are, ALPHA: "<<initialGuess(0)<<", BETA: "<<initialGuess(1)<<", TAU: "<<initialGuess(2)<<"; The maximum log likelihood is "<<loglikelihood<<".\n";
		}
	}

	if(WORKER==0) cout<<"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n\n";
}
