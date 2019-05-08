#include <iostream>

#include "constants.hpp"
#include "class_data.hpp"
#include "class_partition.hpp"

//User parameters
    //Name for the data file.
    std::string DATA_FILE_NAME;

    //The flag for whether to eliminate duplicates in raw data. If it is guaranteed that there are no duplicates in the raw data, set the flag to "false" to save computation time.
    bool ELIMINATION_DUPLICATES_FLAG;

    //The ratio of the smallest distance between any knot and the region boundary to the length of the region in each dimension. Requires a numerical value or "default" to use the default value exp/100.
    double OFFSET;

    //The number of partitions in each level. Requires to be either 2 or 4.
    int NUM_PARTITIONS_J;

    //The number of knots in each region before the finest resolution level. The number of knots at the finest resolution level is automatically determined by the data and the built structure, and the number can be different in different regions at the finest resolution level.
    int NUM_KNOTS_r;

    //The number of levels. Requires a positive integer or "default" to use the number that is automatically determined from NUM_PARTITIONS_J and NUM_KNOTS_r by making the average number of observations per region at the finest resolution level similar to the number of knots in regions before the finest resolution level.
    int NUM_LEVELS_M;

    //The type of calculations. Must be one of "prediction", "optimization", "likelihood", "build_structure_only".
    std::string CALCULATION_MODE;

    //The flag for whether to print the detailed information when running the code.
    bool PRINT_DETAIL_FLAG;

    //The way of specifying the prediction locations. Must be one of 'D' for all the locations in DATA_FILE_NAME no matter whether the associated observation is a valid value or NaN, 'N' for locations in DATA_FILE_NAME that have NaN values, 'A' for locations specified in PREDICTION_LOCATION_FILE.
    char PREDICTION_LOCATION_MODE;

    //Name for the prediction location file. Only used when PREDICTION_LOCATION_MODE='A'.
    std::string PREDICTION_LOCATION_FILE;

    //The flag for whether to dump prediction results to output file.
    bool DUMP_PREDICTION_RESULTS_FLAG;

    //The output file for storing prediction results. Only used when DUMP_PREDICTION_RESULTS_FLAG="true" and CALCULATION_MODE="predict".
    std::string PREDICTION_RESULTS_FILE_NAME;

    //The flag for whether to save temporary files to disk.
	bool SAVE_TO_DISK_FLAG;

    // The flag for the dynamic scheduling of workload for each worker
    bool DYNAMIC_SCHEDULE_FLAG;

    //The directory for saving the temporary files on disk. Only used when "SAVE_TO_DISK"="true".
	std::string TMP_DIRECTORY;

    //Sill parameter in the covariance function. Requires to be a positive numerical value. Will be ignored if CALCULATION_MODE="optimization" where the optimal values that maximizes the likelihood will be found.
	double ALPHA;

    //Range parameter in the covariance function. Requires to be a positive numerical value. Will be ignored if CALCULATION_MODE="optimization" where the optimal values that maximizes the likelihood will be found.
	double BETA;

    //Nugget parameter in the covariance function. Requires to be a positive numerical value. Will be ignored if CALCULATION_MODE="optimization" where the optimal values that maximizes the likelihood will be found.
	double TAU;

    //Maximum number of iterations allowed in the optimization. Only used when CALCULATION_MODE="optimization".
	int MAX_ITERATIONS;

    //Initial guess of sill parameter in the covariance function in optimization. Only used when CALCULATION_MODE="optimization".
	double ALPHA_INITIAL_GUESS;

    //Lower bound for sill parameter in the covariance function in optimization. Only used when CALCULATION_MODE="optimization".
	double ALPHA_LOWER_BOUND;

    //Upper bound for sill parameter in the covariance function in optimization. Only used when CALCULATION_MODE="optimization".
	double ALPHA_UPPER_BOUND;

    //Initial guess of range parameter in the covariance function in optimization. Only used when CALCULATION_MODE="optimization".
	double BETA_INITIAL_GUESS;

    //Lower bound for range parameter in the covariance function in optimization. Only used when CALCULATION_MODE="optimization".
	double BETA_LOWER_BOUND;

    //Upper bound for range parameter in the covariance function in optimization. Only used when CALCULATION_MODE="optimization".
	double BETA_UPPER_BOUND;

    //Initial guess of nugget parameter in the covariance function in optimization. Only used when CALCULATION_MODE="optimization".
	double TAU_INITIAL_GUESS;

    //Lower bound for nugget parameter in the covariance function in optimization. Only used when CALCULATION_MODE="optimization".
	double TAU_LOWER_BOUND;

    //Upper bound for nugget parameter in the covariance function in optimization. Only used when CALCULATION_MODE="optimization".
	double TAU_UPPER_BOUND;

//Other constants
    //Data instance
    Data *data;

    //Partition instance;
    Partition *partition;

    //The total number of unique observation locations that have valid measurements (i.e, not NaN)
    unsigned long NUM_OBSERVATIONS;

    //The total number of unique prediction locations
    unsigned long NUM_PREDICTIONS;

    //Iteration for optimization
    int OPTIMIZATION_ITERATION;

    //Variables for showing execution time
    timeval timeBegin, timeNow;    

    //The starting index of region at each level
    unsigned long *REGION_START;

    //The ending index of region at each level
    unsigned long *REGION_END;