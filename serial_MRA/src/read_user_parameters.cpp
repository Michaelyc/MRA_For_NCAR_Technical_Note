#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "read_user_parameters.hpp"
#include "constants.hpp"
#include "sys/stat.h"

//Function for reading user_parameters
void read_user_parameters()
{
	//Open the configuration file "option"
    std::ifstream file;
    file.open("user_parameters");

    //Read data
    if(file.is_open())
	{
		std::string line;
		while(getline(file,line))
		{
			if(line[0] == '#' || line[0] == 0) continue;

			std::string name = "";
			std::string value = "";
			bool nameOrValue = true; //true for name and false for value
			for(unsigned int iCharacter = 0; iCharacter < line.length(); iCharacter++)
			{
				//Skip irrelevant symbols, only relevant if the character is one of '-', '.', ,'/', '0'-'9', ':'; 'A'-'Z'; 'a'-'z'; "\", "=",
				//if (line[iCharacter]>44&&line[iCharacter]<59) || (line[iCharacter]>64&&line[iCharacter]<91) || (line[iCharacter]>96&&line[iCharacter]<123) || (line[iCharacter]==47) || (line[iCharacter]==92) || (line[iCharacter]==61) || (line[iCharacter]==47)

				//Skip irrelevant symbols including space, tab, ', "
				if(line[iCharacter] == ' ' || line[iCharacter] == '\t' || line[iCharacter] == 39 || line[iCharacter] == 34) continue;

				//Check if reading name or value
				if(line[iCharacter] == '=')
				{
					nameOrValue = false;
					continue;
				}

				//Get the string for name or value
				if(nameOrValue)
					name += line[iCharacter];
				else
					value += line[iCharacter];
			}

			if(name=="") continue;

			//Assign value to name
			if(name == "DATA_FILE_NAME")
			{
				DATA_FILE_NAME = value;
				continue;
			}
			if(name == "PREDICTION_RESULTS_FILE_NAME")
			{
				PREDICTION_RESULTS_FILE_NAME = value;
				continue;
			}
			if(name == "CALCULATION_MODE")
			{
				CALCULATION_MODE = value;
				continue;
			}
			if(name == "PREDICTION_LOCATION_MODE")
			{
				PREDICTION_LOCATION_MODE = value[0];
				continue;
			}
			if(name == "PREDICTION_LOCATION_FILE")
			{
				PREDICTION_LOCATION_FILE = value;
				continue;
			}
			if(name == "PRINT_DETAIL_FLAG")
			{
				std::istringstream tmp(value);
				tmp>>std::boolalpha>>PRINT_DETAIL_FLAG;
				continue;
			}
			if(name == "DUMP_PREDICTION_RESULTS_FLAG")
			{
				std::istringstream tmp(value);
				tmp>>std::boolalpha>>DUMP_PREDICTION_RESULTS_FLAG;
				continue;
			}
			if(name == "SAVE_TO_DISK_FLAG")
			{
				std::istringstream tmp(value);
				tmp>>std::boolalpha>>SAVE_TO_DISK_FLAG;
				continue;
			}
			if(name == "ELIMINATION_DUPLICATES_FLAG")
			{
				std::istringstream tmp(value);
				tmp>>std::boolalpha>>ELIMINATION_DUPLICATES_FLAG;
				continue;
			}
			if(name == "TMP_DIRECTORY")
			{
				TMP_DIRECTORY = value;
				continue;
			}
			if(name == "OFFSET")
			{
				if(value == "default")
				{
					OFFSET = exp(1)/100;
					continue;
				}
				
				std::istringstream tmp(value);
				tmp>>OFFSET;
				
				if(OFFSET<=0 || OFFSET >=1)
				{
					cout<<"Program exits with an error: the value of OFFSET can only be in (0,1), while the currently assigned value is "<<OFFSET<<".\n";
					exit(EXIT_FAILURE);
				}
				continue;
			}
			if(name == "NUM_PARTITIONS_J")
			{
				std::istringstream tmp(value);
				tmp>>NUM_PARTITIONS_J;
				continue;
			}
			if(name == "NUM_KNOTS_r")
			{
				std::istringstream tmp(value);
				tmp>>NUM_KNOTS_r;
				continue;
			}
			if(name == "NUM_LEVELS_M")
			{
				if(value=="default") 
					NUM_LEVELS_M = -99;//-99 stands for the default method to automatically determine NUM_LEVELS_M  
				else
				{
					std::istringstream tmp(value);
					tmp>>NUM_LEVELS_M;
				}
				continue;
			}
			if(name == "ALPHA")
			{
				std::istringstream tmp(value);
				tmp>>ALPHA;
				continue;
			}
			if(name == "BETA")
			{
				std::istringstream tmp(value);
				tmp>>BETA;
				continue;
			}
			if(name == "TAU")
			{
				std::istringstream tmp(value);
				tmp>>TAU;
				continue;
			}
			if(name == "MAX_ITERATIONS")
			{
				std::istringstream tmp(value);
				tmp>>MAX_ITERATIONS;
				continue;
			}
			if(name == "ALPHA_LOWER_BOUND")
			{
				std::istringstream tmp(value);
				tmp>>ALPHA_LOWER_BOUND;
				continue;
			}
			if(name == "BETA_LOWER_BOUND")
			{
				std::istringstream tmp(value);
				tmp>>BETA_LOWER_BOUND;
				continue;
			}
			if(name == "TAU_LOWER_BOUND")
			{
				std::istringstream tmp(value);
				tmp>>TAU_LOWER_BOUND;
				continue;
			}
			if(name == "ALPHA_UPPER_BOUND")
			{
				std::istringstream tmp(value);
				tmp>>ALPHA_UPPER_BOUND;
				continue;
			}
			if(name == "BETA_UPPER_BOUND")
			{
				std::istringstream tmp(value);
				tmp>>BETA_UPPER_BOUND;
				continue;
			}
			if(name == "TAU_UPPER_BOUND")
			{
				std::istringstream tmp(value);
				tmp>>TAU_UPPER_BOUND;
				continue;
			}			
			if(name == "ALPHA_INITIAL_GUESS")
			{
				std::istringstream tmp(value);
				tmp>>ALPHA_INITIAL_GUESS;
				continue;
			}
			if(name == "BETA_INITIAL_GUESS")
			{
				std::istringstream tmp(value);
				tmp>>BETA_INITIAL_GUESS;
				continue;
			}
			if(name == "TAU_INITIAL_GUESS")
			{
				std::istringstream tmp(value);
				tmp>>TAU_INITIAL_GUESS;
				continue;
			}			
			cout<<"Program exits with an error: The variable with name \""+name+"\" in the configuration file \"user_parameters\" is not defined. Check again with the variable name by looking at the comments.\n";
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		cout<<"Program exits with an error: fail to open the configuration file \"user_parameters\"\n";
		exit(EXIT_FAILURE);
	}

	if(NUM_PARTITIONS_J!=2&&NUM_PARTITIONS_J!=4)
	{
		cout<<"Program exits with an error: the specified NUM_PARTITIONS_J is "<<NUM_PARTITIONS_J<<". Only implemented for NUM_PARTITIONS_J equal to 2 or 4.\n";
		exit(EXIT_FAILURE);
	}

	if(NUM_LEVELS_M < 2 && NUM_LEVELS_M != -99)//-99 stands for the default method to automatically determine NUM_LEVELS_M 
	{
		cout<<"Program exits with an error: the specified NUM_LEVELS_M is "<<NUM_LEVELS_M<<", which is required to be larger than 1.\n";
		exit(EXIT_FAILURE);
	}

	if(CALCULATION_MODE!="prediction"&&CALCULATION_MODE!="optimization"&&CALCULATION_MODE!="likelihood"&&CALCULATION_MODE!="build_structure_only")
	{
		cout<<"Program exits with an error: the specified CALCULATION_MODE is \""<<CALCULATION_MODE<<"\", which must be one of \"predict\", \"optimization\", \"likelihood\", and \"build_structure_only\".\n";
		exit(EXIT_FAILURE);
	}

	if(PREDICTION_LOCATION_MODE!='D'&&PREDICTION_LOCATION_MODE!='N'&&PREDICTION_LOCATION_MODE!='A')
	{
		cout<<"Program exits with an error: the specified PREDICTION_LOCATION_MODE is '"<<PREDICTION_LOCATION_MODE<<"', which must be one of 'D', 'N', and 'A'\n";
		exit(EXIT_FAILURE);
	}

	if(SAVE_TO_DISK_FLAG)
	{
		//Check whether The temporary directory that stores variables in runtime exists
		struct stat directoryStatus;
		if(stat(TMP_DIRECTORY.c_str(),&directoryStatus) == -1)
		{
			cout<<"Program exits with an error: the directory "<<TMP_DIRECTORY<<" does not exist. Please create it before running the program.\n";
			exit(EXIT_FAILURE);
		}
		
		if(!S_ISDIR(directoryStatus.st_mode))
		{
			cout<<"Program exits with an error: "<<TMP_DIRECTORY<<" is not a directory. Please create a directory with name "<<TMP_DIRECTORY<<" before running the program.\n";
			exit(EXIT_FAILURE);
		}
	}
	
	if(CALCULATION_MODE == "optimization")
	{
		if(ALPHA_UPPER_BOUND<ALPHA_LOWER_BOUND)
		{
			cout<<"Program exits with an error: the specified ALPHA_UPPER_BOUND "<<ALPHA_UPPER_BOUND<<" is smaller than the specified ALPHA_LOWER_BOUND "<<ALPHA_LOWER_BOUND<<".\n";
			exit(EXIT_FAILURE);
		}
		if(BETA_UPPER_BOUND<BETA_LOWER_BOUND)
		{
			cout<<"Program exits with an error: the specified BETA_UPPER_BOUND "<<BETA_UPPER_BOUND<<" is smaller than the specified BETA_LOWER_BOUND "<<BETA_LOWER_BOUND<<".\n";
			exit(EXIT_FAILURE);
		}
		if(TAU_UPPER_BOUND<TAU_LOWER_BOUND)
		{
			cout<<"Program exits with an error: the specified TAU_UPPER_BOUND "<<TAU_UPPER_BOUND<<" is smaller than the specified TAU_LOWER_BOUND "<<TAU_LOWER_BOUND<<".\n";
			exit(EXIT_FAILURE);
		}

		if(ALPHA_INITIAL_GUESS<ALPHA_LOWER_BOUND)
		{
			cout<<"Program exits with an error: the specified ALPHA_INITIAL_GUESS "<<ALPHA_INITIAL_GUESS<<" is smaller than the specified ALPHA_LOWER_BOUND "<<ALPHA_LOWER_BOUND<<".\n";
			exit(EXIT_FAILURE);
		}
		if(BETA_INITIAL_GUESS<BETA_LOWER_BOUND)
		{
			cout<<"Program exits with an error: the specified BETA_UPPER_BOUND "<<BETA_INITIAL_GUESS<<" is smaller than the specified BETA_LOWER_BOUND "<<BETA_LOWER_BOUND<<".\n";
			exit(EXIT_FAILURE);
		}
		if(TAU_INITIAL_GUESS<TAU_LOWER_BOUND)
		{
			cout<<"Program exits with an error: the specified TAU_INITIAL_GUESS "<<TAU_INITIAL_GUESS<<" is smaller than the specified TAU_LOWER_BOUND "<<TAU_LOWER_BOUND<<".\n";
			exit(EXIT_FAILURE);
		}

		if(ALPHA_INITIAL_GUESS>ALPHA_UPPER_BOUND)
		{
			cout<<"Program exits with an error: the specified ALPHA_INITIAL_GUESS "<<ALPHA_INITIAL_GUESS<<" is smaller than the specified ALPHA_UPPER_BOUND "<<ALPHA_UPPER_BOUND<<".\n";
			exit(EXIT_FAILURE);
		}
		if(BETA_INITIAL_GUESS>BETA_UPPER_BOUND)
		{
			cout<<"Program exits with an error: the specified BETA_UPPER_BOUND "<<BETA_INITIAL_GUESS<<" is smaller than the specified BETA_UPPER_BOUND "<<BETA_UPPER_BOUND<<".\n";
			exit(EXIT_FAILURE);
		}
		if(TAU_INITIAL_GUESS>TAU_UPPER_BOUND)
		{
			cout<<"Program exits with an error: the specified TAU_INITIAL_GUESS "<<TAU_INITIAL_GUESS<<" is smaller than the specified TAU_UPPER_BOUND "<<TAU_UPPER_BOUND<<".\n";
			exit(EXIT_FAILURE);
		}
	}
}