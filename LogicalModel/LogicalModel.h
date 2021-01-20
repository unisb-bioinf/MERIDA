

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <DenseMatrixReader.h>
#include <DenseMatrix.h>
#include <DenseMatrixSubset.h>

#ifndef LOGICALMODEL_LOGICALMODEL_H
#define LOGICALMODEL_LOGICALMODEL_H

class LogicalModel{

    auto calculateNoWeights(GeneTrail::DenseMatrixSubset& mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void; //ic50_values and threshold can be invalid, because they are not used anyway
    auto calculateLinearWeights(GeneTrail::DenseMatrixSubset& mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void;
    auto calculateQuadraticWeights(GeneTrail::DenseMatrixSubset& mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void;
    auto calculateCubicWeights(GeneTrail::DenseMatrixSubset& mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void;
    auto calculateExponentialWeights(GeneTrail::DenseMatrixSubset& mtx, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void;
    auto calculateMisclassifcationRate(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>& prediction) -> double;
    virtual auto calculateObjectiveFunction(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>& prediction) -> double;

protected:
    bool valid_file;
    int number_of_samples;
    int number_of_features;
    GeneTrail::DenseMatrix data;
    GeneTrail::DenseMatrixSubset data_mtx;
    auto doubleEpsilonEqual(double a, double b) -> bool;

public:
    LogicalModel(std::string filename);
    LogicalModel(GeneTrail::DenseMatrix & mtx);
    LogicalModel(GeneTrail::DenseMatrixSubset & data_mtx);

    virtual auto buildModel() -> bool = 0;
    virtual auto solveModel() -> void = 0;
    virtual auto predict(GeneTrail::DenseMatrixSubset& mtx) -> std::vector<int> = 0;
    auto calculateWeights(GeneTrail::DenseMatrixSubset& mtx, std::string weightMode, GeneTrail::DenseMatrix& ic50_values, double threshold) ->  void; //actually changes the given matrix
    auto calculateWeights(std::string weightMode, GeneTrail::DenseMatrix& ic50_values, double threshold) -> void; //actually changes the matrix used in this class (data_mtx)
    auto calculateError(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>& prediction, std::string errorMode) -> double;
    auto calculateError(std::vector<int>& prediction, std::string errorMode) -> double;
    auto calculateSensitivity(GeneTrail::DenseMatrixSubset& mtx, std::vector<int> & prediction) -> double;
    auto calculateSpecificity(GeneTrail::DenseMatrixSubset& mtx, std::vector<int> & prediction) -> double;
    auto calculateTP(GeneTrail::DenseMatrixSubset& mtx, std::vector<int> &prediction) -> double;
    auto calculateFP(GeneTrail::DenseMatrixSubset& mtx, std::vector<int> &prediction) -> double;
    auto calculateTN(GeneTrail::DenseMatrixSubset& mtx, std::vector<int> &prediction) -> double;
    auto calculateFN(GeneTrail::DenseMatrixSubset& mtx, std::vector<int> &prediction) -> double;

};

#endif
