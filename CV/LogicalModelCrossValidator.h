//
// Created by klenhof on 13.03.18.
//


#ifndef ACT_LOMO_LOGICALMODELCROSSVALIDATOR_H
#define ACT_LOMO_LOGICALMODELCROSSVALIDATOR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <DenseMatrixReader.h>
#include <DenseMatrix.h>
#include <algorithm>
#include <math.h>
#include <matrixTools.h>
#include "CrossValidator.h"


class LogicalModelCrossValidator: public CrossValidator<LogicalModelCrossValidator>{
    std::string cvMode;//perform x-fold cross validation or cross validation leaving out x samples?
    std::string weightMode; //which function should be used to calculate weight?
    std::string errorMode; //which error function should be used for cv folds?
    std::string matrixFile; //file where full data matrix is stored
    std::string ic50File; //file where ic50 values of cell lines are stored
    std::string resultDirectory; //result directory
    std::string pathway_definition_filename; //filename of a file with pathway definitions in gmt format
    double sensitivity; //sensitivity and specificity are currently not used, but could be used to control logical model
    double specificity;
    double threshold; //binarization threshold for IC50 values
    double alpha; //parameter for SRS2
    int first_p; //first parameter for respective logical model
    int second_p; //second parameter for respective logical model
    int number_of_samples;
    int number_of_features;
    bool valid; //is set to false if one of the files could not be opened or read
    std::string errorMsg;
    std::string preselected_features;
    auto readParamFile() -> void;
    auto readData(GeneTrail::DenseMatrix & mtx, GeneTrail::DenseMatrix & ic50vals) -> void;
    auto printError(std::vector<std::vector<double>>& error) -> void;

public:
    LogicalModelCrossValidator(std::string filename, int foldsize): CrossValidator(filename, foldsize), valid(true){ std::cout << filename << "!!!" << std::endl;}
    template<typename Logmod>
    auto crossvalidate() -> void{

        readParamFile();


        std::cout << "Your parameters are " << std::endl;
        std::cout << "CV modus: " << cvMode << std::endl;
        std::cout << "Weight function: " << weightMode << std::endl;
        std::cout << "Error function: " << errorMode << std::endl;
        std::cout << "First parameter (specificity, lambda1): " << specificity << std::endl; //TODO: clean this up
        std::cout << "Second parameter (specificity, lambda2: " << sensitivity << std::endl;
        std::cout << "Matrix filename: " << matrixFile << std::endl;
        std::cout << "IC50 value filename: " << ic50File << std::endl;
        std::cout << "M1: " << first_p << std::endl;
        std::cout << "M2: " << second_p << std::endl;
        std::cout << "Threshold: " << threshold << std::endl;
        std::cout << "Alpha: " << alpha << std::endl;
        std::cout << "Preselected features: " << preselected_features << std::endl;
        std::cout << "PWDatabase: " << pathway_definition_filename << std::endl;


        //Read in data from data file
        GeneTrail::DenseMatrix mtx(1,1);
        GeneTrail::DenseMatrix ic50vals(1,1);

        readData(mtx, ic50vals);

        if(!valid){
            std::cout << "Cross validation not possible" << std::endl;
            std::cout << errorMsg << std::endl;
            return;
        }

        //Randomly shuffle data because it could be ordered
        std::vector <unsigned int> indices(number_of_samples);
        unsigned int count{0};
        for(auto& i : indices){
            i = count;
            ++count;
        }
        std::random_shuffle(indices.begin(), indices.end());

        //print out content
        /*std::cout << "myvector contains: ";
        for(std::vector<unsigned int>::iterator it=indices.begin(); it!=indices.end(); ++it){
            std::cout << ' ' << *it;
            std::cout <<'\n';
        }*/

        mtx.shuffleRows(indices);

        //Make CV folds according to foldsize and available samples
        //If CV mode is "fold" then foldsize must be calculated

        if(cvMode == "fold"){
            foldsize = int(floor(double(number_of_samples)/double(foldsize))); //foldsize is cv_folds in this case and must be changed to foldsize
            std::cout << "FOLDSIZE!!!!!!!!!!!!!!!!!!! " << foldsize << std::endl;
        }

        int cv_folds(floor(double(number_of_samples)/double(foldsize)));
        //Error vector for folds plus final averaged error
        std::vector<std::vector<double>> fold_error(cv_folds, std::vector<double>(8));

        for(int i = 0; i<cv_folds; ++i){

            //Training set
            std::vector<unsigned int> training_rows{};
            for(int j = 0; j< i*foldsize; ++j){
                training_rows.push_back(j);
            }

            for(int j = i*foldsize + foldsize; j< number_of_samples; j++){
                training_rows.push_back(j);
            }

            //Test set
            std::vector<unsigned int> test_rows{};

            for(int j = i*foldsize; j < i*foldsize + foldsize; j++) {
                test_rows.push_back(j);
            }

            std::cout << "Fold " << i << std::endl;
            std::cout << "Training rows " << std::endl;
            //Convert indices to names
            std::vector<std::string> tr_rows_names{};
            for(auto j: training_rows){
                tr_rows_names.push_back(mtx.rowName(j));
                std::cout << mtx.rowName(j) << std::endl;

            }

            std::cout << "Test rows" << std::endl;
            std::vector<std::string> test_rows_names{};
            for(auto j: test_rows){
                test_rows_names.push_back(mtx.rowName(j));
                std::cout << mtx.rowName(j) << std::endl;
            }

            std::cout << "Length of test row vector " << test_rows_names.size() << std::endl;
            //Create training and test set
            std::tuple<GeneTrail::DenseMatrixSubset, GeneTrail::DenseMatrixSubset> train_test{GeneTrail::splitMatrixRows(mtx, tr_rows_names, test_rows_names)};

            std::cout << "Going to build model" << std::endl;
            //Train model on training data set
            Logmod train_mdl(first_p, second_p, sensitivity, specificity, alpha, std::get<0>(train_test), preselected_features, pathway_definition_filename); //TODO: remove first_p and second_p when time is right :D

            std::cout << "Trained model" << std::endl;
            train_mdl.calculateWeights(weightMode, ic50vals,threshold);
            std::cout << "Calculated weights" << std::endl;
            //Need to use internal code information -> would be nice to change this one day
            std::string fold_i{resultDirectory + "Fold" + std::to_string(i)};
            train_mdl.setResultDirectory(fold_i);
            std::cout << "Next: build model" << std::endl;
            train_mdl.buildModel();
            std::cout << "Next: get model" << std::endl;
            train_mdl.getModel(); //delete if tested
            train_mdl.setDebugstrm(false);
            std::cout << "Next: solve model" << std::endl;
            train_mdl.solveModel();

            std::cout << "Came here!!!!!" << std::endl;
            if(!train_mdl.isFeasible()){
                (fold_error[i])[0] = nan("");
                (fold_error[i])[1] = nan("");
                (fold_error[i])[2] = nan("");
                (fold_error[i])[3] = nan("");
                (fold_error[i])[4] = nan("");
                (fold_error[i])[5] = nan("");
                (fold_error[i])[6] = nan("");
                (fold_error[i])[7] = nan("");
                std::cout << "Infeasible ILP for fold " << i << std::endl;
                break;
            }

            //Predict output for build model
            GeneTrail::DenseMatrixSubset& test_set{std::get<1>(train_test)};

            train_mdl.calculateWeights(test_set, weightMode, ic50vals, threshold);

            std::vector<int> prediction(train_mdl.predict(test_set));


            /*std::cout << "Prediction for fold " << i << std::endl;

            for(int j= 0; j< prediction.size(); j++){
                std::cout << "sample " <<  j << "\t";
                std::cout << prediction[j] << std::endl;
            }*/
            //Report error, sensitivity, and specificity
            (fold_error[i])[0] = train_mdl.calculateError(test_set, prediction, "objectiveFunction");
            (fold_error[i])[1] = train_mdl.calculateError(test_set, prediction, "misclassification");
            (fold_error[i])[2] = train_mdl.calculateSensitivity(test_set, prediction);
            (fold_error[i])[3] = train_mdl.calculateSpecificity(test_set, prediction);
            (fold_error[i])[4] = train_mdl.calculateTN(test_set, prediction);
            (fold_error[i])[5] = train_mdl.calculateTP(test_set, prediction);
            (fold_error[i])[6] = train_mdl.calculateFN(test_set, prediction);
            (fold_error[i])[7] = train_mdl.calculateFP(test_set, prediction);


            std::cout << "fold " <<  i << "\t" << "objectivef, miscl, sensitivity, specificity, TN, TP, FN, FP" << "\t";
            std::cout << fold_error[i][0] << "," <<  fold_error[i][1] << "," << fold_error[i][2] << "," << fold_error[i][3];
            std::cout << fold_error[i][4] << "," << fold_error[i][5] << "," << fold_error[i][6] << "," << fold_error[i][7] << std::endl;



        }

        printError(fold_error);

        //Train final model
        Logmod final_mdl(first_p, second_p, sensitivity, specificity, alpha,  mtx, preselected_features, pathway_definition_filename);//TODO: one day the first_p and second_p parameter should be removed --> define a class that holds an object with all the data in it



        final_mdl.calculateWeights(weightMode, ic50vals, threshold);
        final_mdl.setResultDirectory(resultDirectory);
        std::cout << "Next: build model" << std::endl;
        final_mdl.buildModel();
        std::cout << "Next: get model" << std::endl;
        final_mdl.getModel(); //delete if tested
        final_mdl.setDebugstrm(false);
        std::cout << "Next: solve model" << std::endl;
        final_mdl.solveModel();
    }

};


#endif //ACT_LOMO_LOGICALMODELCROSSVALIDATOR_H
