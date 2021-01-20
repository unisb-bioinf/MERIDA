//
// Created by klenhof on 05.03.19.
//



#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <DenseMatrixReader.h>
#include <DenseMatrix.h>
#include <DenseMatrixSubset.h>
#include <LogicalModel.h>
#define IL_STD
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#ifndef MERIDA_ILP_MERIDA_ILP_H
#define MERIDA_ILP_MERIDA_ILP_H


class MERIDA_ILP: public LogicalModel {
    int max_feat;
    bool testo;
    bool sol_requested; //needed to ensure that solutions are only generated once
    bool feasible;
    double elapsed_time;
    std::string result_directory;
    std::ofstream out;
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloNumVarArray x_i;
    IloNumVarArray y_i;
    IloNumVarArray s_k;
    IloNumVarArray r_k;
    IloNumVarArray z_k;
    IloNumArray sol_x_i;
    IloNumArray sol_y_i;
    IloNumArray sol_s_k;
    IloNumArray sol_r_k;
    std::vector<std::string> sensitive_features;
    std::vector<std::string> resistant_features;
    GeneTrail::DenseMatrix preselected_features;

    auto requestSol() -> void;
    auto printSol() -> void;
    auto printHeader(std::ofstream & output) -> void;
    auto printShortForm(std::ofstream & output) -> void;
    auto print_x_i(std::ofstream & output) -> void;
    auto print_y_i(std::ofstream & output) -> void;
    auto print_s_k(std::ofstream & output) -> void;
    auto print_r_k(std::ofstream & output) -> void;
    auto print_y_k(std::ofstream & output) -> void;

    auto calculateObjectiveFunction(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>& prediction) -> double override;

public:
    MERIDA_ILP(std::string filename, std::string result_directory,  std::string preselected_feature_file, int max_feat);
    MERIDA_ILP(int max_feat, GeneTrail::DenseMatrixSubset & mtx,  std::string preselected_feature_file); //now data is not any useful structure
    MERIDA_ILP(int max_feat,  GeneTrail::DenseMatrix & mtx,  std::string preselected_feature_file);
    MERIDA_ILP(int u1, int u2, double lambda_1, double lambda_2, double alpha, GeneTrail::DenseMatrixSubset & mtx, std::string pr, std::string pwd):
            MERIDA_ILP(u1, mtx, pr){} //now data is not any useful structure
    MERIDA_ILP(int u1, int u2, double lambda_1, double lambda_2, double alpha, GeneTrail::DenseMatrix & mtx, std::string pr, std::string pwd):
            MERIDA_ILP(u1, mtx, pr){}
    ~MERIDA_ILP();
    //Forbidden until decided otherwise
    MERIDA_ILP(MERIDA_ILP & srs) = delete; //copy ctor
    MERIDA_ILP(MERIDA_ILP && srs) = delete; //move ctor
    MERIDA_ILP & operator=(MERIDA_ILP && srs) = delete; //move operator

    auto fileValid() -> bool;
    auto getDataMatrix() -> GeneTrail::DenseMatrix&; //should be DenseMatrixSubset
    auto getResultDirectory() -> std::string;
    auto setServer(bool b) -> void;
    auto getModel() -> void;
    auto setResultDirectory(std::string dir) -> void;
    auto isFeasible() -> bool;
    /*Must be used after buildModel but before solveModel
     * redirects output into file with name 'name'
     * can append info to existing file if desired*/
    auto setDebugstrm(bool append) -> void;
    auto buildModel() -> bool;
    auto solveModel() -> void;
    auto mytime() -> double;
    auto totalMemory() -> int;
    auto usedMemory() -> int;
    auto predict(GeneTrail::DenseMatrixSubset& mtx) -> std::vector<int>;
};


#endif //MERIDA_ILP_MERIDA_ILP_H
