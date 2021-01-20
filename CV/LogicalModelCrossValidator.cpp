//
// Created by klenhof on 13.03.18.
//

#include "LogicalModelCrossValidator.h"

auto LogicalModelCrossValidator::readParamFile() -> void {

    std::ifstream data_strm(paramfile);
    std::string line;

    std::cout << "name of parameter file : " << paramfile << std::endl;
    std::cout << "?????????????????????????????????????" << std::endl;
    //read from parameter file
    if(data_strm.is_open()){

        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        //First parameter
        getline(data_strm, line);
        std::stringstream lstrm(line);
        std::string parameter_name;
        std::getline(lstrm, parameter_name, '\t'); //not used at the moment
        std::getline(lstrm, matrixFile, '\n'); // will not be found but this is not neccessary anyway
        std::cout  << matrixFile;

        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Second parameter (result directory) was not found in the file, file is too small\n";
            return;

        }

        //Second parameter
        getline(data_strm, line);
        std::stringstream lstrm2(line);
        std::getline(lstrm2, parameter_name, '\t'); //not used at the moment
        std::getline(lstrm2, resultDirectory, '\n'); // will not be found but this is not neccessary anyway
        std::cout << resultDirectory;

        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Third parameter (first parameter for logical model) was not found in the file, file is too small.\n";
            return;

        }


        //Third parameter
        getline(data_strm, line);
        std::stringstream lstrm3(line);
        std::getline(lstrm3, parameter_name, '\t'); //not used at the moment
        lstrm3 >> first_p;

        std::cout << first_p;
        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Fourth parameter (second parameter for logical model) was not found in the file, file is too small\n";
            return;

        }

        //Fourth parameter
        getline(data_strm, line);
        std::stringstream lstrm4(line);
        std::getline(lstrm4, parameter_name, '\t'); //not used at the moment
        lstrm4 >> second_p;

        std::cout << second_p;
        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Fifth parameter (sensitivity) was not found in the file, file is too small.\n";
            return;

        }

        //Fifth parameter
        getline(data_strm, line);
        std::stringstream lstrm5(line);
        std::getline(lstrm5, parameter_name, '\t'); //not used at the moment
        lstrm5 >> sensitivity;

        std::cout << sensitivity;

        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Sixth parameter (specificity) was not found in the file, file is too small.\n";
            return;

        }

        //Sixth parameter
        getline(data_strm, line);
        std::stringstream lstrm6(line);
        std::getline(lstrm6, parameter_name, '\t'); //not used at the moment
        lstrm6 >> specificity;

        std::cout << specificity;

        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Seventh parameter (weight calculation) was not found in the file, file is too small.\n";
            return;

        }

        //Seventh parameter
        getline(data_strm, line);
        std::stringstream lstrm7(line);
        std::getline(lstrm7, parameter_name, '\t'); //not used at the moment
        lstrm7 >> weightMode;


        std::cout << weightMode;
        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Eigth parameter (error calculation) was not found in the file, file is too small.\n";
            return;

        }

        //Eigth parameter
        getline(data_strm, line);
        std::stringstream lstrm8(line);
        std::getline(lstrm8, parameter_name, '\t'); //not used at the moment
        lstrm8 >> errorMode;

        std::cout << errorMode;
        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Ninth parameter (binarization threshold) was not found in the file, file is too small\n";
            return;

        }

        //Ninth parameter
        getline(data_strm, line);
        std::stringstream lstrm9(line);
        std::getline(lstrm9, parameter_name, '\t'); //not used at the moment
        lstrm9 >> threshold;

        std::cout << threshold;
        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Tenth parameter (IC50 file) was not found in the file, file is too small.\n";
            return;

        }

        //Tenth parameter
        getline(data_strm, line);
        std::stringstream lstrm10(line);
        std::getline(lstrm10, parameter_name, '\t'); //not used at the moment
        std::getline(lstrm10, ic50File, '\n');

        std::cout << ic50File;

        if(!data_strm || data_strm.eof()){
            valid = false;
            errorMsg += "Eleventh parameter (CV mode) was not found in file, file is too small.\n";
            return;
        }
        //Eleventh parameter
        getline(data_strm, line);
        std::stringstream lstrm12(line);
        std::getline(lstrm12, parameter_name, '\t');
        std::getline(lstrm12, cvMode, '\n');

        if(getline(data_strm, line)){
            std::stringstream lstrm11(line);
            std::getline(lstrm11, parameter_name, '\t');
            lstrm11 >> alpha;
            std::cout << "First case " << alpha;
        }
        else{
            alpha = 0.0;
            std::cout << "Second case " << alpha;
        }

        if(getline(data_strm, line)){
            std::stringstream lstrm13(line);
            std::getline(lstrm13, parameter_name, '\t');
            std::getline(lstrm13, preselected_features, '\t');
        }
        else{
            preselected_features = "none";
        }

        if(getline(data_strm, line)){
            std::stringstream lstrm14(line);
            std::getline(lstrm14, parameter_name, '\t');
            std::getline(lstrm14, pathway_definition_filename, '\t');
        }
        else{
            pathway_definition_filename = "none";
        }

        if(getline(data_strm, line)){
            std::cout << "Warning: Not all parameters in your parameter file are used." << std::endl;
            return;
        }

    }

    else{
        valid = false;
        errorMsg + "Could not open the file " + paramfile + "\n";
        return;
    }
    return;
}

auto LogicalModelCrossValidator::readData(GeneTrail::DenseMatrix & mtx, GeneTrail::DenseMatrix & ic50vals) -> void {

    std::ifstream data_strm(matrixFile);
    GeneTrail::DenseMatrixReader reader;

    if(!data_strm){
        errorMsg +=  "Filename for matrix file was invalid.\n";
        errorMsg +=  "Filename was " +  matrixFile + "\n";
        valid = false;
    }
    else{
        mtx = reader.read(data_strm, GeneTrail::DenseMatrixReader::READ_COL_NAMES | GeneTrail::DenseMatrixReader::READ_ROW_NAMES);
        number_of_features = mtx.cols() -2; //weights and response must be subtracted
        number_of_samples = mtx.rows();

    }

    data_strm.close();

    std::ifstream ic50_strm(ic50File);
    GeneTrail::DenseMatrixReader reader2;

    if(!ic50_strm){
        errorMsg +=  "Filename for ic50 file was invalid.\n";
        errorMsg +=  "Filename was " +  ic50File + "\n";
        valid = false;
    }
    else{
        ic50vals = reader2.read(ic50_strm, GeneTrail::DenseMatrixReader::READ_COL_NAMES | GeneTrail::DenseMatrixReader::READ_ROW_NAMES);
        if(ic50vals.cols()!=1){
            errorMsg += "IC50 value matrix is invalid\n";
            errorMsg += "Matrix had " + std::to_string(ic50vals.cols()) + " many columns but needs 1. \n";
            valid = false;
        }

    }
    ic50_strm.close();
    return;
}


auto LogicalModelCrossValidator::printError(std::vector<std::vector<double>>& error) -> void{

    std::string delim("_");
    std::string filename{};
    if(first_p == 0 && second_p == 0){

        if(alpha == 0.0){ //SRS model
            filename = resultDirectory + "CVerror" + delim + delim + "Lambda1" + std::to_string(sensitivity) + delim + "Lambda2" + std::to_string(specificity) + delim + "Feat" + std::to_string(number_of_features) + delim + "Samp" + std::to_string(number_of_samples) + ".txt";
        }
        else{ //SRS2 model
            filename =  filename = resultDirectory + "CVerror" + delim + delim + "Lambda1" + std::to_string(sensitivity) + delim + "Lambda2" + std::to_string(specificity) + delim +  "Alph" + std::to_string(alpha) +  delim + "Feat" + std::to_string(number_of_features) + delim + "Samp" + std::to_string(number_of_samples) + ".txt";
        }

    }
    else{
        filename = resultDirectory + "CVerror" + delim + "FirstParameter" + std::to_string(first_p) +
                   delim + "SecondParameter" + std::to_string(second_p) +
                   delim + "Sens" + std::to_string(sensitivity) +
                   delim + "Spec" + std::to_string(specificity) +
                   delim + "Feat" + std::to_string(number_of_features) +
                   delim + "Samp" + std::to_string(number_of_samples) + ".txt";
    }

    std::cout << "CV output has been written to " << filename << std::endl;
    std::ofstream output(filename);

    output << "Error_function\tfold\terror" << std::endl;
    for(unsigned int i = 0; i<error.size(); ++i){
        output << "ObjF\t" << i  << "\t" << error[i][0] << "\n";
    }

    for(unsigned int i = 0; i<error.size(); ++i){
        output << "Misclassification\t" << i << "\t" << error[i][1] << "\n";
    }

    for(unsigned int i = 0; i<error.size(); ++i){
        output << "Sensitivity\t" << i << "\t" << error[i][2] << "\n";
    }

    for(unsigned int i = 0; i<error.size(); ++i){
        output << "Specificity\t" << i << "\t" << error[i][3] << "\n";
    }

    for(unsigned int i = 0; i<error.size(); ++i){
        output << "TN\t" << i << "\t" << error[i][4] << "\n";
    }

    for(unsigned int i = 0; i<error.size(); ++i){
        output << "TP\t" << i << "\t" << error[i][5] << "\n";
    }

    for(unsigned int i = 0; i<error.size(); ++i){
        output << "FN\t" << i << "\t" << error[i][6] << "\n";
    }

    for(unsigned int i = 0; i<error.size(); ++i){
        output << "FP\t" << i << "\t" << error[i][7] << "\n";
    }

    output.close();
}
