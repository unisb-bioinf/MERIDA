#include <iostream>
#include <sstream>
#include <DenseMatrixReader.h>
#include <DenseMatrix.h>
#include <fstream>
#include "MERIDA_ILP.h"
#include "LogicalModelCrossValidator.h"

//what is expected from parameter file
std::string data_file;
std::string result_directory;
int number_of_sens_elements;
int number_of_res_elements;
std::string weightFunction;
double threshold;
std::string ic50_file;
std::string pre_selected_features;

/*
 * Checks whether the given parameters are valid
 * Returns boolean
 * */
bool checkParameters(){

    return true;
}

/*
 * Parses the parameters from the parameter file
 * Returns whether the parameters could be successfully parsed
 * */
bool parseParameters(int argc, char** argv ){
    if(argc != 3){
        std::cout << "Please provide the name of a file with parameter settings, and whether you want to perform cross validation." << std::endl;
    }

    std::cout << "Entered parseParameters" << std::endl;
    std::ifstream data_strm(argv[1]);
    std::string line;

    //read from parameter file
    if(data_strm.is_open()){



        if(!data_strm || data_strm.eof()){
            std::cout << "First field in the parameter setting file should be the name of a data file." << std::endl;
            return false;

        }

        //First parameter
        getline(data_strm, line);
        std::stringstream lstrm(line);
        std::string parameter_name;
        std::getline(lstrm, parameter_name, '\t'); //not used at the moment
        std::getline(lstrm, data_file, '\n'); // \n will not be found but this is not neccessary anyway



        if(!data_strm || data_strm.eof()){
            std::cout << "Second field in the parameter setting file should be the directory for the results." << std::endl;
            return false;

        }

        //Second parameter
        getline(data_strm, line);
        std::stringstream lstrm2(line);
        std::getline(lstrm2, parameter_name, '\t'); //not used at the moment
        std::getline(lstrm2, result_directory, '\n'); // will not be found but this is not neccessary anyway



        if(!data_strm || data_strm.eof()){
            std::cout << "Third field in the parameter setting file should be a positive integer." << std::endl;
            return false;

        }
        //Third parameter
        getline(data_strm, line);
        std::stringstream lstrm3(line);
        std::getline(lstrm3, parameter_name, '\t'); //not used at the moment
        lstrm3 >> number_of_sens_elements;


        if(!data_strm || data_strm.eof()){
            std::cout << "Fourth field in the parameter setting file should be a positive integer." << std::endl;
            return false;

        }

        //Fourth parameter
        getline(data_strm, line);
        std::stringstream lstrm4(line);
        std::getline(lstrm4, parameter_name, '\t'); //not used at the moment
        lstrm4 >> number_of_res_elements;

        if(!data_strm || data_strm.eof()){
            std::cout << "Fifth field in the parameter setting file should be either the the weight function to be used or no." << std::endl;
            return false;

        }

        //Fifth parameter
        getline(data_strm, line);
        std::stringstream lstrm7(line);
        std::getline(lstrm7, parameter_name, '\t'); //not used at the moment
        lstrm7 >> weightFunction;


        if(!data_strm || data_strm.eof()){
            std::cout << "Sixth field in the parameter setting file should be either the path to a file with ic50 values for the cell lines in the data matrix or no in case you provide the weights by yourself." << std::endl;
            return false;

        }
        //Sixth parameter
        getline(data_strm, line);
        std::stringstream lstrm8(line);
        std::getline(lstrm8, parameter_name, '\t'); //not used at the moment
        std::getline(lstrm8, ic50_file, '\n'); // will not be found but this is not neccessary anyway



        if(!data_strm || data_strm.eof()){
            std::cout << "Seventh field in the parameter setting file should either be the threshold between sensitive and resistant cell lines or 0.0 in case you provide the weights by yourself." << std::endl;
            return false;

        }
        //Seventh parameter
        getline(data_strm, line);
        std::stringstream lstrm9(line);
        std::getline(lstrm9, parameter_name, '\t'); //not used at the moment
        lstrm9 >> threshold;




        if(!data_strm || data_strm.eof()){
            std::cout << "Eighth field in the parameter setting file should be the file name of a file for preselected features." << std::endl;
            return false;

        }

        //Eighth parameter
        getline(data_strm, line);
        std::stringstream lstrm11(line);
        std::getline(lstrm11, parameter_name, '\t'); //not used at the moment
        std::getline(lstrm11, pre_selected_features, '\n'); // \n will not be found but this is not neccessary anyway

        if(getline(data_strm, line)){
            std::cout << "Warning: not all fields in your file are used." << std::endl;
        }

        //Sanity check
        return checkParameters();

    }

    else{
        std::cout << "Could not open the file " << argv[1] << std::endl;
        return false;
    }


}

int main(int argc, char** argv) {


    for(int i = 0; i < argc; i++){
        std::cout << argv[i] << "------";
    }
    std::cout << std::endl;
    if(argc < 3){
        std::cout << "This program can be used in two ways: \n"
                "- either for the training of a logical model, then a parameter setting file is required and the second argument should be a 'no' \n"
                "- or for cross validating, then a different parameter setting file is required, the second argument is a 'yes' and the third "
                "argument is the numer of folds for cross validation.\n";
        return 12;
    }

    std::string cv_yes_no(argv[2]);
    if(cv_yes_no == "yes"){
        if(argc != 4){
            std::cout << argc << std::endl;
            std::cout << "This program can be used in two ways: \n"
                    "- either for the training of a logical model, then a parameter setting file is required and the second argument should be a 'no'\n"
                    "- or for cross validating, then a different parameter setting file is required, the second argument is a 'yes' and the third "
                    "argument is the numer of folds for cross validation.\n";
            return 12;
        }
        else{
            LogicalModelCrossValidator cv{argv[1], atoi(argv[3])};
            cv.crossvalidate<MERIDA_ILP>();
            return 0;
        }
    }

    //If no cross validation is performed
    //Check and read parameters
    if(!parseParameters(argc, argv)){
        std::cout << "Please modify the parameter setting file such that calculations can be performed." << std::endl;
        return 1;
    }


    //Start LOBICO computations/model building process
    std::cout << "Your parameters were" << std::endl;
    std::cout << "Max number of sensitivity determining elements: " << number_of_sens_elements << std::endl;
    std::cout << "Max number of resistance determining elements: " << number_of_res_elements << std::endl;
    std::cout << "Pre-selected features file: " << pre_selected_features << std::endl;
    std::cout << "Output directory " << result_directory << std::endl;

    std::cout << "Starting MERIDA_ILP computations" << std::endl;


    MERIDA_ILP mdl(data_file, result_directory, pre_selected_features, number_of_sens_elements);//number of sensitive elements and number of resistant elements should be equal, only one is needed
    if(mdl.fileValid()){
        std::cout << "Next: build model" << std::endl;

        if(weightFunction != "no"){
            std::cout << "Should be here" << std::endl;
            GeneTrail::DenseMatrix ic50_values(1,1);

            std::ifstream ic50_strm(ic50_file);
            GeneTrail::DenseMatrixReader reader2;

            if(!ic50_strm){
                std::cout << "IC50 filename was invalid. " << std::endl;
                std::cout << "filename was " << ic50_file << std::endl;
                return 0;
            }
            else{
                ic50_values = reader2.read(ic50_strm, GeneTrail::DenseMatrixReader::READ_COL_NAMES | GeneTrail::DenseMatrixReader::READ_ROW_NAMES);
                if(ic50_values.cols()!=1){
                    std::cout << "IC50 value matrix is invalid." << std::endl;
                    std::cout << "Matrix had " + std::to_string(ic50_values.cols()) + " many columns but needs 1." << std::endl;
                    return 0;
                }

            }
            ic50_strm.close();
            mdl.calculateWeights(weightFunction, ic50_values, threshold);
        }

        mdl.buildModel();
        std::cout << "Next: get model" << std::endl;
        mdl.getModel();
        mdl.setDebugstrm(false);
        std::cout << "Next: solve model" << std::endl;
        mdl.solveModel();
        std::cout << "End of main!" << std::endl;
    }
    else{
        std::cout << "Filename inside parameter file was invalid" << std::endl;
        std::cout << "Filename was " <<  data_file << std::endl;
    }
    return 0;
}
