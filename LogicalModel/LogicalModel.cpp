#include "LogicalModel.h"

LogicalModel::LogicalModel(std::string filename): data(1,1),  data_mtx(){


    //Read in data from data file
    std::ifstream data_strm(filename);
    GeneTrail::DenseMatrixReader reader;

    if(!data_strm){
        valid_file = false;
    }
    else{
        valid_file = true;
        data = reader.read(data_strm, GeneTrail::DenseMatrixReader::READ_COL_NAMES | GeneTrail::DenseMatrixReader::READ_ROW_NAMES);
        number_of_features = data.cols() -2; //weights and response must be subtracted
        number_of_samples = data.rows();
        std::vector<unsigned int> s(number_of_samples);
        unsigned int count{0};
        for(auto& i : s){
            i = count;
            ++count;
        }

        std::vector<unsigned int> f(number_of_features + 2);
        unsigned int count2{0};
        for(auto& i: f){
            i = count2;
            ++count2;
        }

        data_mtx = std::move(GeneTrail::DenseMatrixSubset(&data, s, f));


    }


    data_strm.close();

}

LogicalModel::LogicalModel(GeneTrail::DenseMatrix & mtx):  valid_file(true), data(mtx), data_mtx(){
    number_of_features = data.cols() -2; //weights and response must be subtracted
    number_of_samples = data.rows();
    std::vector<unsigned int> s(number_of_samples);
    unsigned int count{0};
    for(auto& i : s){
        i = count;
        ++count;
    }
    std::vector<unsigned int> f(number_of_features + 2);
    unsigned int count2{0};
    for(auto& i: f){
        i = count2;
        ++count2;
    }
    data_mtx = std::move(GeneTrail::DenseMatrixSubset(&data, s, f));
}

LogicalModel::LogicalModel(GeneTrail::DenseMatrixSubset & mtx): valid_file(true), data(1,1), data_mtx(mtx) {

    number_of_features = data_mtx.cols()-2;
    number_of_samples = data_mtx.rows();
}



auto LogicalModel::calculateWeights(GeneTrail::DenseMatrixSubset& mtx, std::string weightMode, GeneTrail::DenseMatrix& ic50_values, double threshold) ->  void {

    if(weightMode == "none"){
        calculateNoWeights(mtx, ic50_values, threshold);
        return;
    }

    if(weightMode == "linear"){
        calculateLinearWeights(mtx, ic50_values, threshold);
        return;
    }
    if(weightMode == "quadratic"){
        calculateQuadraticWeights(mtx, ic50_values, threshold);
        return;
    }
    if(weightMode == "cubic"){
        calculateCubicWeights(mtx, ic50_values, threshold);
        return;
    }
    if(weightMode == "exponential"){
        calculateExponentialWeights(mtx, ic50_values, threshold);
        return;
    }
    else{
        std::cout << "..." << weightMode << "..." << std::endl;
        std::cout << "Only no, linear, quadratic, cubic and exponential weights supported, calculate linear weights now" << std::endl;
        calculateLinearWeights(mtx, ic50_values, threshold);
        return;
    }

}

auto LogicalModel::calculateNoWeights(GeneTrail::DenseMatrixSubset & mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void{

    for(unsigned int i = 0; i< mtx.rows(); ++i){
        mtx(i, mtx.cols()-2) = 1.0;
    }
}

auto LogicalModel::calculateLinearWeights(GeneTrail::DenseMatrixSubset& mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void{

    std::vector<double> intermediate_results(mtx.rows());

    double sensitive{0.0};
    double resistant{0.0};

    for(unsigned int i = 0; i< mtx.rows(); ++i){
        //std::cout << "???" << std::endl;
        //std::cout << i << std::endl;
        //std::cout << mtx.rowName(i) << std::endl;
        //std::cout << ic50_values.rowIndex(mtx.rowName(i)) << std::endl;
        intermediate_results[i] = std::fabs(ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold);
        //std::cout << i << std::endl;
        //std::cout << intermediate_results[i] << std::endl;

        if(doubleEpsilonEqual(mtx(i, mtx.cols()-1), 1.0)){
            sensitive += std::fabs(ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold);
        }
        else{
            resistant += std::fabs(ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold);
        }
    }


    for(unsigned int i = 0; i< mtx.rows(); ++i){
        if(doubleEpsilonEqual(mtx(i, mtx.cols()-1), 1.0)){
            mtx(i, mtx.cols()-2) = intermediate_results[i]/ (2*sensitive);
            //std::cout << "weight for sample " << i  << " " << mtx(i, mtx.cols()-2) << std::endl;
        }
        else{
            mtx(i, mtx.cols()-2) = intermediate_results[i]/ (2*resistant);
            //std::cout << "weight for sample " << i << " " << mtx(i, mtx.cols()-2) << std::endl;
        }
    }

}

auto LogicalModel::calculateQuadraticWeights(GeneTrail::DenseMatrixSubset& mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void{


    std::vector<double> intermediate_results(mtx.rows());

    double sensitive{0.0};
    double resistant{0.0};

    for(unsigned int i = 0; i< mtx.rows(); ++i){

        intermediate_results[i] = pow(ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold, 2);

        //std::cout << intermediate_results[i] << std::endl;
        if(doubleEpsilonEqual(mtx(i, mtx.cols()-1), 1.0)){
            sensitive += pow(ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold,2);
        }
        else{
            resistant += pow(ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold,2);
        }
    }

    for(unsigned int i = 0; i< mtx.rows(); ++i){
        if(doubleEpsilonEqual(mtx(i, mtx.cols()-1), 1.0)){
            mtx(i, mtx.cols()-2) = intermediate_results[i]/ (2*sensitive);
            //std::cout << "weight for sample " << i  << " " << mtx(i, mtx.cols()-2) << std::endl;
        }
        else{
            mtx(i, mtx.cols()-2) = intermediate_results[i]/ (2*resistant);
            //std::cout << "weight for sample " << i  << " " << mtx(i, mtx.cols()-2) << std::endl;
        }
    }


}

auto LogicalModel::calculateCubicWeights(GeneTrail::DenseMatrixSubset & mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void{


    std::vector<double> intermediate_results(mtx.rows());

    double sensitive{0.0};
    double resistant{0.0};

    for(unsigned int i = 0; i< mtx.rows(); ++i){

        intermediate_results[i] = pow(std::fabs(ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold), 3);

        std::cout << intermediate_results[i] << std::endl;
        if(doubleEpsilonEqual(mtx(i, mtx.cols()-1), 1.0)){
            sensitive += pow(std::fabs(ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold),3);
        }
        else{
            resistant += pow(std::fabs(ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold),3);
        }
    }

    for(unsigned int i = 0; i< mtx.rows(); ++i){
        if(doubleEpsilonEqual(mtx(i, mtx.cols()-1), 1.0)){
            mtx(i, mtx.cols()-2) = intermediate_results[i]/ (2*sensitive);
            std::cout << "weight for sample " << i  << " " << mtx(i, mtx.cols()-2) << std::endl;
        }
        else{
            mtx(i, mtx.cols()-2) = intermediate_results[i]/ (2*resistant);
            std::cout << "weight for sample " << i  << " " << mtx(i, mtx.cols()-2) << std::endl;
        }
    }
}

auto LogicalModel::calculateExponentialWeights(GeneTrail::DenseMatrixSubset & mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void{


    std::vector<double> intermediate_results(mtx.rows());

    double sensitive{0.0};
    double resistant{0.0};

    for(unsigned int i = 0; i< mtx.rows(); ++i){

        intermediate_results[i] = pow(10, ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold);

        std::cout << intermediate_results[i] << std::endl;
        if(doubleEpsilonEqual(mtx(i, mtx.cols()-1), 1.0)){
            sensitive += pow(10, ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold);
        }
        else{
            resistant += pow(10, ic50_values(ic50_values.rowIndex(mtx.rowName(i)), 0) - threshold);
        }
    }

    for(unsigned int i = 0; i< mtx.rows(); ++i){
        if(doubleEpsilonEqual(mtx(i, mtx.cols()-1), 1.0)){
            mtx(i, mtx.cols()-2) = intermediate_results[i]/ (2*sensitive);
            std::cout << "weight for sample " << i  << " " << mtx(i, mtx.cols()-2) << std::endl;
        }
        else{
            mtx(i, mtx.cols()-2) = intermediate_results[i]/ (2*resistant);
            std::cout << "weight for sample " << i  << " " << mtx(i, mtx.cols()-2) << std::endl;
        }
    }
}
auto LogicalModel::calculateWeights(std::string weightMode, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void {

    calculateWeights(data_mtx, weightMode, ic50_values, threshold);
    return;

}


auto LogicalModel::calculateMisclassifcationRate(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>& prediction) -> double{

    double miscl{0.0};
    std::cout <<  "Number of rows in prediction matrix " << mtx_response.rows() << std::endl;
    for(unsigned int i = 0; i< mtx_response.rows(); ++i){


        if(!((doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1),1.0)) && (prediction[i] == 1)) && !((doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1),0.0)) && (prediction[i] == 0))){
            std::cout << "Entered the misclassification part" << std::endl;
            miscl+=1.0;
        }
    }
    return miscl/mtx_response.rows();
}

auto LogicalModel::calculateObjectiveFunction(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>& prediction) -> double{

    double resistant{0.0};
    double sensitive{0.0};
    for(unsigned int i = 0; i< mtx_response.rows(); ++i){

        if(doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1),0.0)){
            resistant +=  mtx_response(i, mtx_response.cols()-2)*prediction[i];
        }
        else{
            sensitive += mtx_response(i, mtx_response.cols()-2)*prediction[i];
        }
    }
    return resistant-sensitive;
}


auto LogicalModel::calculateError(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>& prediction, std::string errorMode) -> double{

    if(errorMode == "misclassification"){
        return calculateMisclassifcationRate(mtx_response, prediction);
    }
    if(errorMode == "objectiveFunction"){
        return calculateObjectiveFunction(mtx_response, prediction);
    }
    else{
        std::cout << "Only misclassificaton and objectiveFunction are allowed, calculate misclassification rate now" << std::endl;
        return calculateMisclassifcationRate(mtx_response, prediction);
    }
}


auto LogicalModel::calculateError(std::vector<int>& prediction, std::string errorMode) -> double{


    if(errorMode == "misclassification"){
        return calculateMisclassifcationRate(data_mtx, prediction);
    }
    if(errorMode == "objectiveFunction"){
        return calculateObjectiveFunction(data_mtx, prediction);
    }
    else{
        std::cout << "Only misclassificaton and objectiveFunction are allowed, calculate misclassification rate now" << std::endl;
        return calculateMisclassifcationRate(data_mtx, prediction);
    }
}

auto LogicalModel::calculateSensitivity(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>&prediction) -> double{

    double nr_sensitives(0.0);
    double sensitivity(0.0);
    for(int i = 0; i< mtx_response.rows(); i++){

        if(doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1), 1.0)){
            nr_sensitives += 1.0;

            if(prediction[i] == 1){
                sensitivity += 1.0;
            }
        }

    }

    if(doubleEpsilonEqual(nr_sensitives, 0.0)){
        return 1;
    }
    return sensitivity/nr_sensitives;

}

auto LogicalModel::calculateSpecificity(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>&prediction) -> double{

    double nr_resistants(0.0);
    double specificity(0.0);

    for(int i = 0; i< mtx_response.rows(); i++){
        if(doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1),0.0)){
            nr_resistants += 1.0;

            if(prediction[i] == 0){
                specificity += 1.0;
            }
        }
    }

    if(doubleEpsilonEqual(nr_resistants, 0.0)){
        return 1;
    }

    return specificity/nr_resistants;
}


auto LogicalModel::calculateTP(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>&prediction) -> double{

    double TP(0.0);
    for(int i = 0; i< mtx_response.rows(); i++){

        if(doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1), 1.0)){

            if(prediction[i] == 1){
                TP += 1.0;
            }
        }

    }


    return TP;


}

auto LogicalModel::calculateFP(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>&prediction) -> double{

    double FP(0.0);
    for(int i = 0; i< mtx_response.rows(); i++){

        if(doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1), 0.0)){

            if(prediction[i] == 1){
                FP += 1.0;
            }
        }

    }


    return FP;


}

auto LogicalModel::calculateTN(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>&prediction) -> double{

    double TN(0.0);
    for(int i = 0; i< mtx_response.rows(); i++){

        if(doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1), 0.0)){

            if(prediction[i] == 0){
                TN += 1.0;
            }
        }

    }
    return TN;
}

auto LogicalModel::calculateFN(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>&prediction) -> double{

    double FN(0.0);
    for(int i = 0; i< mtx_response.rows(); i++){

        if(doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1), 1.0)){

            if(prediction[i] == 0){
                FN += 1.0;
            }
        }

    }
    return FN;
}




auto LogicalModel::doubleEpsilonEqual(double a, double b) -> bool{

    return std::fabs(a-b) < 0.1;
}
