
## MERIDA

Welcome to the GitHub page of MERIDA (MEthod  for  Rule Identification  with multi-omics  DAta), a tool for drug sensitivity prediction in cancer.

For issues and questions, please contact Kerstin Lenhof (klenhof[at]bioinf.uni-sb.de)
If you use MERIDA or the code in this repository, please cite out paper

## Installation
To install MERIDA, you need to have the following programs/libraries installed:
- Eigen3
- Cplex
- Boost (including the libraries)

If everything is installed, create a build directory, change into the created build directory
and issue cmake:

mkdir build
cd build
cmake -DEIGEN_PATH=/path/to/eigen3 -DCPLEX_PATH=/path/to/cplex ..


## Usage
MERIDA can be used in two ways: 
- Without using the implemented CV functionality: here, a config file and a 'no' as second parameter is required as input (see Example_Config1.txt)
- Using the implemented CV functionality: then a different parameter setting file is required, the second argument is a 'yes' and the third argument is the numer of folds for cross validation (see Example_Config2.txt)
