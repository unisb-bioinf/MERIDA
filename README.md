


## Installation
To install MERIDA, you need to have the following programs/libraries installed:
- Eigen3
- Cplex
- Boost (including the libraries)

If everything is installed, create a build directory, change into the created build directory
and issue cmake:

mkdir build
cd build
cmake -DEIGEN_PATH="/usr/include/eigen3" -DCPLEX_PATH="/share/data/opt/opt/cplex1262" ..


## Usage
MERIDA can be used in two ways: 
- Without using the implemented CV functionality: here, a config file and a 'no' as second parameter is required as input (see Example_Config1.txt)
- Using the implemented CV functionality: then a different parameter setting file is required, the second argument is a 'yes' and the third argument is the numer of folds for cross validation (see Example_Config2.txt)
