#include <iostream>
#include <math.h>
#include <sys/time.h>

#include "class_data.hpp"
#include "class_partition.hpp"
#include "class_approximation.hpp"
#include "read_user_parameters.hpp"
#include "constants.hpp"
#include "optimization.hpp"


int main(int argc, char* argv[])
{
	time_t currentTime;
	time(&currentTime);

	gettimeofday(&timeBegin, NULL);
	cout<<"\n===========================> Program starts at: "<<ctime(&currentTime);

	//Read user parameters
	read_user_parameters();

	//Load data
	data = new Data;
	data->load_data();

	//Build partition
	partition = new Partition;
	partition->build_partition();

	//Do the calculation mode
	if(CALCULATION_MODE == "build_structure_only") partition->dump_structure_information();
	
	if(CALCULATION_MODE == "likelihood")
	{
		Approximation approximation(ALPHA,BETA,TAU);
		approximation.likelihood();
	}
	
	if(CALCULATION_MODE=="prediction")
	{
		Approximation approximation(ALPHA,BETA,TAU);
		approximation.likelihood();
		approximation.predict();
		if(DUMP_PREDICTION_RESULTS_FLAG) approximation.dump_prediction_result();
	}

	if(CALCULATION_MODE == "optimization")
		get_optimal_parameters();

	time(&currentTime);
	gettimeofday(&timeNow, NULL);
	cout<<"===========================> Program ends at: "<<ctime(&currentTime);
	cout<<"===========================> Running time: "<<(double)timeNow.tv_sec-(double)timeBegin.tv_sec+((double)timeNow.tv_usec-(double)timeBegin.tv_usec)/1000000.0<<" seconds.\n\n";
	
	return 0;
}
