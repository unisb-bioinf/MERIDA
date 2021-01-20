//
// Created by klenhof on 05.03.19.
//


#include "MERIDA_ILP.h"





MERIDA_ILP::MERIDA_ILP(std::string filename, std::string result_directory,  std::string preselected_features_file, int max_feat): LogicalModel(filename),
                                                                                                                                   max_feat(max_feat), testo(false), sol_requested(false), feasible(false), result_directory(result_directory), model(env), cplex(env),
                                                                                                                                   x_i(env), y_i(env), s_k(env), r_k(env), z_k(env),
                                                                                                                                   sol_x_i(env), sol_y_i(env), sol_s_k(env), sol_r_k(env), sensitive_features(), resistant_features(), preselected_features(1,1){

   
    //Read in data from data file

    if(preselected_features_file.empty()){
        return;

    }

    else{
        std::ifstream data_strm(preselected_features_file);
        GeneTrail::DenseMatrixReader reader;

        if(!data_strm){
            std::cout << "Preselected feature file was invalid!" << std::endl;
        }
        else{
            preselected_features = reader.read(data_strm, GeneTrail::DenseMatrixReader::READ_COL_NAMES | GeneTrail::DenseMatrixReader::READ_ROW_NAMES);

        }


        data_strm.close();
    }


}

MERIDA_ILP::MERIDA_ILP(int max_feat, GeneTrail::DenseMatrix & mtx, std::string preselected_features_file): LogicalModel(mtx),
                                                                                                            max_feat(max_feat), testo(false), sol_requested(false), feasible(false),
                                                                                                            model(env), cplex(env), x_i(env), y_i(env), s_k(env), r_k(env), z_k(env), sol_x_i(env), sol_y_i(env), sol_s_k(env), sol_r_k(env), sensitive_features(), resistant_features(), preselected_features(1,1){


    //Read in data from data file

    if(preselected_features_file.empty()){//no features pre selected
        return;
    }

    else{
        std::ifstream data_strm(preselected_features_file);
        GeneTrail::DenseMatrixReader reader;

        if(!data_strm){
            std::cout << "Preselected feature file was invalid!" << std::endl;
        }
        else{
            preselected_features = reader.read(data_strm, GeneTrail::DenseMatrixReader::READ_COL_NAMES | GeneTrail::DenseMatrixReader::READ_ROW_NAMES);

        }


        data_strm.close();

    }


}

MERIDA_ILP::MERIDA_ILP(int max_feat, GeneTrail::DenseMatrixSubset & mtx, std::string preselected_features_file): LogicalModel(mtx),
                                                                                                                  max_feat(max_feat), testo(false), sol_requested(false), feasible(false),
                                                                                                                  model(env), cplex(env), x_i(env), y_i(env), s_k(env), r_k(env), z_k(env),
                                                                                                                  sol_x_i(env), sol_y_i(env), sol_s_k(env), sol_r_k(env), sensitive_features(), resistant_features(), preselected_features(1,1){

    //Read in data from data file

    if(preselected_features_file.empty()){ //no features pre selected
        return;
    }

    else{
        std::ifstream data_strm(preselected_features_file);
        GeneTrail::DenseMatrixReader reader;

        if(!data_strm){
            std::cout << "Preselected feature file was invalid!" << std::endl;
        }
        else{
            preselected_features = reader.read(data_strm, GeneTrail::DenseMatrixReader::READ_COL_NAMES | GeneTrail::DenseMatrixReader::READ_ROW_NAMES);

        }


        data_strm.close();

    }


}



MERIDA_ILP::~MERIDA_ILP() {
    if(out.is_open()){
        out.close();
    }
    env.end();
}

auto MERIDA_ILP::getDataMatrix() -> GeneTrail::DenseMatrix& {
    return data; //needs to be changed one day
}


auto MERIDA_ILP::getResultDirectory() -> std::string {
    return result_directory;
}

auto MERIDA_ILP::fileValid() -> bool {
    return valid_file;
}

auto MERIDA_ILP::getModel() -> void {

    std::string filename{result_directory +
                         "M_" + std::to_string(max_feat) +
                         "Feat_" + std::to_string(number_of_features) +
                         "Samp_" + std::to_string(number_of_samples) + "Model.lp"};
    cplex.extract(model);

    try{
        cplex.exportModel(filename.c_str());
    }
    catch(IloCplex::Exception & e){
        std::cout << e.getMessage();
    }

}

auto MERIDA_ILP::buildModel() -> bool {



    //Precompute all possible things
    std::vector<std::vector<std::string>> samples;

    //List of cell lines with their associated alteration in the matrix
    for(int i = 0; i< data_mtx.rows(); i++){

        std::vector<std::string> alterations;
        for(int j = 0; j< number_of_features; j++){
            if(doubleEpsilonEqual(data_mtx(i,j), 1.0)){
                alterations.push_back(data_mtx.colName(j));
            }
        }

        samples.push_back(std::move(alterations));
    }


    //Precalculate stuff
    std::vector<std::string> G;
    for(int i = 0; i< number_of_features; i++){
        G.push_back(data_mtx.colName(i));
    }

    //Build ILP

    //x_i and y_i
    x_i = IloNumVarArray(env, number_of_features);
    y_i = IloNumVarArray(env, number_of_features);

    for(int i = 0; i< number_of_features; i++){
        x_i[i] = IloBoolVar(env, G.at(i).c_str());
    }
    for(int i = 0; i< number_of_features; i++){
        y_i[i] = IloBoolVar(env, G.at(i).c_str());
    }

    //Constraints for x_i and y_i
    for(int i = 0; i< number_of_features; i++){
        IloExpr expr(env);
        expr += x_i[i] + y_i[i] -1;
        model.add(expr <= 0);
        expr.end();
    }

    //Only a fixed number of elements is allowed to be chosen for x_i and y_i
    IloExpr expr_sens_feat(env);
    IloExpr expr_res_feat(env);
    for(int i = 0; i< number_of_features; i++){
        expr_sens_feat += x_i[i];
        expr_res_feat += y_i[i];
    }

    /**if(max_sens < preselected_features.rows()){
        max_sens = max_sens + preselected_features.rows();
    }
    if(max_res < preselected_features.rows()){
        max_res = max_res + preselected_features.rows();
    }**/ //TODO: must be fixed otherwise

    model.add(expr_sens_feat + expr_res_feat <= max_feat);
    expr_sens_feat.end();
    expr_res_feat.end();


    //Special constraints as imposed by new matrix
    for(int i = 0; i< preselected_features.cols(); i++){

        for(int j = 0; j< preselected_features.rows(); j++){
            if(i == 0){
                if( doubleEpsilonEqual(preselected_features(j,i), 1.0)){
                    for(int g = 0; g< G.size(); g++){
                        if(x_i[g].getName() == preselected_features.rowName(j)){
                            IloExpr exprx(env);
                            exprx += x_i[g];
                            model.add(exprx == 1);
                            exprx.end();
                        }
                    }
                }
            }
            if(i == 1){
                if( doubleEpsilonEqual(preselected_features(j,i), 1.0)){
                    for(int g = 0; g< G.size(); g++){
                        if(x_i[g].getName() == preselected_features.rowName(j)){
                            IloExpr expry(env);
                            expry += x_i[g];
                            model.add(expry == 0);
                            expry.end();
                        }
                    }
                }
            }

            if(i == 2){
                if( doubleEpsilonEqual(preselected_features(j,i), 1.0)){
                    for(int g = 0; g< G.size(); g++){
                        if(y_i[g].getName() == preselected_features.rowName(j)){
                            IloExpr exprz(env);
                            exprz += y_i[g];
                            model.add(exprz == 1);
                            exprz.end();
                        }
                    }
                }
            }
            if(i == 3){
                if( doubleEpsilonEqual(preselected_features(j,i), 1.0)){
                    for(int g = 0; g< G.size(); g++){
                        if(y_i[g].getName() == preselected_features.rowName(j)){
                            IloExpr exprw(env);
                            exprw += y_i[g];
                            model.add(exprw == 0);
                            exprw.end();
                        }
                    }
                }
            }

        }
    }

    // std::cout << "WTF" << std::endl;
    //s_k and r_k
    s_k = IloNumVarArray(env, number_of_samples);
    r_k = IloNumVarArray(env, number_of_samples);

    for(int k = 0; k< number_of_samples; k++){
        s_k[k] = IloBoolVar(env, data_mtx.rowName(k).c_str());
        r_k[k] = IloBoolVar(env, data_mtx.rowName(k).c_str());
    }

    //Constraints for s_k and r_k
    for(int k = 0; k < number_of_samples; k++){
        IloExpr expr2(env);


        //s_k
        for(int i = 0; i< G.size(); i++){

            std::string elem(G[i]);
            std::vector<std::string>::iterator pos_elem;

            pos_elem = std::find(samples.at(k).begin(), samples.at(k).end(), elem);
            if(!(pos_elem == samples.at(k).end())){
                expr2 += x_i[i];
            }

        }
        model.add(s_k[k] <= expr2);
        model.add(expr2 <= s_k[k] * (int) samples.at(k).size() );
        expr2.end();

        //r_k
        IloExpr expr3(env);
        for(int i = 0; i< G.size(); i++){

            std::string elem(G[i]);
            std::vector<std::string>::iterator pos_elem;

            pos_elem = std::find(samples.at(k).begin(), samples.at(k).end(), elem);
            if(!(pos_elem == samples.at(k).end())){
                expr3 += y_i[i];
            }

        }
        model.add(r_k[k] <= expr3);
        model.add(expr3 <= r_k[k] * (int) samples.at(k).size());
        expr3.end();

    }



    //z_k
    z_k = IloNumVarArray(env, number_of_samples);

    for(int k = 0; k< number_of_samples; k++){
        z_k[k] = IloBoolVar(env, data_mtx.rowName(k).c_str());
    }
    //Constraints for z_k
    for(int i = 0; i < number_of_samples; i++){
        IloExpr expr4(env);

        model.add(0 <= s_k[i] + 1- r_k[i] -2*z_k[i]);
        model.add(s_k[i] +1 -r_k[i] - 2*z_k[i] <=1);
        expr4.end();

    }

    //Objective function
    IloExpr obj(env);


    for(int i = 0; i< number_of_samples; i++){
        if(doubleEpsilonEqual(data_mtx(i, data_mtx.cols()-1), 1.0)){
            obj += data_mtx(i, data_mtx.cols()-2)*z_k[i] - data_mtx(i, data_mtx.cols()-2) + data_mtx(i, data_mtx.cols()-2) * z_k[i] ; //weighted (TP-FN)
            //obj += z_k[i] * (1.0/double(n_sens)) - (1.0/double(n_sens)) + z_k[i] * (1.0/ double(n_sens)); // (TP - FN)/n_sens
        }
        else{
            obj += data_mtx(i, data_mtx.cols()-2) - data_mtx(i, data_mtx.cols()-2)*z_k[i] - data_mtx(i, data_mtx.cols()-2)*z_k[i]; //weighted (TN - FP)
            //obj += 1 *(1.0/double(n_res))-z_k[i]*(1.0/ double(n_res)) - z_k[i] *(1.0/double(n_res)); //(TN - FP)/n_res
        }
    }


    model.add(IloMaximize(env, obj));

    return true;
}

auto MERIDA_ILP::solveModel() -> void {


    cplex.extract(model);

    //Running on testosterone?

    cplex.setParam(IloCplex::Threads, 32); //max number of threads
    cplex.setParam(IloCplex::ParallelMode, 1); //deterministic calculations enforced
    //cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.05); //needs only to be used if calculation takes too long


    try{
        //cplex.setParam(IloCplex::ClockType, 2);//wall clock time, cplex does not work ...
        time_t timer;
        timer = time(NULL);
        //double first = cplex.getTime();
        feasible = cplex.solve();

        std::cout << "ILP successfully solved" << std::endl;
        time_t timer2;
        timer2 = time(NULL);
        //double second = cplex.getTime();
        elapsed_time = difftime(timer2, timer);
    }
    catch(IloException & e){
        std::cout << "CPLEX raised an exception" << std::endl;
        std::cout << e << std::endl;
        e.end();

    }
    catch(std::exception & e){
        std::cout << "Exception occurred" << std::endl;
        std::cout << e.what() << std::endl;
    }

    if(!feasible){
        std::cout << "ILP was not feasible, additional information is written to Not_feasible.txt" << std::endl;
        std::ofstream output("Not_feasible.txt");
        output << "ILP was not feasible" << std::endl;
        output << "Status\t" << cplex.getStatus() << std::endl;
        output << "Solver Status\t" << cplex.getCplexStatus() << std::endl;
        output.close();
    }
    else{
        std::cout << "ILP was feasible, solution and solution information have been written to separate files" << std::endl;
        std::ofstream sol_info("Feasible.txt");
        sol_info << "ILP was feasible" << std::endl;
        sol_info << "Status\t" << cplex.getStatus() << std::endl;
        sol_info << "Solver Status\t" << cplex.getCplexStatus() << std::endl;
        sol_info.close();
        //Solution as is by cplex
        cplex.writeSolution("Solution.sol");
        //Solution mapped back to features
        printSol();

    }
}

auto MERIDA_ILP::isFeasible() -> bool {
    return feasible;
}

auto MERIDA_ILP::printSol() -> void {
    //Request solution information
    requestSol();



    std::string delim("_");
    std::string filename(result_directory + "Result" +
                         delim + "M_" + std::to_string(max_feat) +
                         delim + "Feat_" + std::to_string(number_of_features) +
                         delim + "Samp_" + std::to_string(number_of_samples) + ".txt");

    std::cout << "Results have been written to " << filename << std::endl;
    std::ofstream output(filename);

    printHeader(output);
    printShortForm(output);

    print_x_i(output);
    print_y_i(output);
    print_s_k(output);
    print_r_k(output);
    print_y_k(output);

    output.close();
}

auto MERIDA_ILP::printHeader(std::ofstream & output) -> void {

    output << "MERIDA_ILP analysis" << std::endl;
    output << "Parameters"<< std::endl;
    output << "M: " << std::to_string(max_feat) << std::endl;
    double obj_value = cplex.getObjValue(); //gecracked...
    output <<  "Objective_value_of_best_solution: " << std::to_string(obj_value) << std::endl;
    output << "\n";


}

auto MERIDA_ILP::print_x_i(std::ofstream & output) -> void {

    for(int i = 0; i < x_i.getSize(); i++){
        output << x_i[i].getName() << "\t" << sol_x_i[i] << "\n";
    }
}

auto MERIDA_ILP::print_y_i(std::ofstream & output) -> void {

    for(int i = 0; i< y_i.getSize(); i++){
        output << y_i[i].getName() << "\t" << sol_y_i[i] << "\n";
    }
}

auto MERIDA_ILP::print_s_k(std::ofstream & output) -> void {

    for(int k = 0; k< s_k.getSize(); k++){
        output << "s_" << std::to_string(k) << "\t" << s_k[k].getName() << "\t" << sol_s_k[k] << "\n";
    }
}

auto MERIDA_ILP::print_r_k(std::ofstream & output) -> void {

    for(int k = 0; k< r_k.getSize(); k++){
        output << "r_" << std::to_string(k) << "\t" << r_k[k].getName() << "\t" << sol_r_k[k] << "\n";
    }
}

auto MERIDA_ILP::mytime() -> double{
    return elapsed_time;
}

//How much memory was allocated for future use
auto MERIDA_ILP::totalMemory() -> int{
    return env.getTotalMemoryUsage();
}

//How much memory is actually used
auto MERIDA_ILP::usedMemory() -> int{
    return env.getMemoryUsage();
}

auto MERIDA_ILP::setServer(bool b) -> void{
    testo = b;
}

auto MERIDA_ILP::setDebugstrm(bool append) ->void{

    std::string name{result_directory + "Debug" +
                     "M_" + std::to_string(max_feat) + 
                     "Feat_" + std::to_string(number_of_features) +
                     "Samp_" + std::to_string(number_of_samples) +
                     ".txt"};
    if(append){
        out.open(name, std::ios_base::app);
    }
    else{
        out.open(name);
    }

    if(!out.good()){
        throw invalid_argument("Debug file could not be opened, std::cout will be used instead");
    }
    else{

        cplex.setOut(out);

    }
}

auto MERIDA_ILP::setResultDirectory(std::string dir) -> void {
    result_directory = dir;
}



auto MERIDA_ILP::requestSol()-> void {
    //Request solution information

    if(sol_requested){ //solution has already been generated
        return;
    }

    cplex.getValues(sol_x_i, x_i);
    cplex.getValues(sol_y_i, y_i);
    cplex.getValues(sol_s_k, s_k);
    cplex.getValues(sol_r_k, r_k);

    //Create model in short form
    for(int i = 0; i< sol_x_i.getSize(); i++){
        if(doubleEpsilonEqual(sol_x_i[i], 1.0)){
            sensitive_features.emplace_back(x_i[i].getName());
        }
    }
    for(int i = 0; i < sol_y_i.getSize(); i++){
        if(doubleEpsilonEqual(sol_y_i[i], 1.0)){
            resistant_features.emplace_back(y_i[i].getName());
        }
    }

    sol_requested = true;

}
auto MERIDA_ILP::printShortForm(std::ofstream & output) -> void{

    requestSol();

    std::string logical_model{};
    logical_model += "Sensitivity_features: ";
    for(auto feature: sensitive_features){
        logical_model += feature + "\t";

    }
    logical_model += "Resistance features: ";
    for(auto feature: resistant_features){
        logical_model += feature + "\t";
    }
    logical_model += "\n";

    output << "Resulting model: " << logical_model << std::endl;

}

auto MERIDA_ILP::predict(GeneTrail::DenseMatrixSubset& mtx) -> std::vector<int> {
    requestSol();

    std::vector<int> prediction(mtx.rows());

    for (unsigned int i = 0; i < mtx.rows(); ++i) {
        bool one_sens{false};
        bool one_res{false};
        for(auto s_feature: sensitive_features){
            if(doubleEpsilonEqual(mtx(i, mtx.colIndex(s_feature)), 1.0)){
                one_sens = true;
            }
        }
        for(auto r_feature: resistant_features){

            if(doubleEpsilonEqual(mtx(i, mtx.colIndex(r_feature)), 1.0)){
                one_res = true;
            }
        }
        if(one_res == true){
            prediction[i] = 0;
        }
        else{
            if(one_sens == true){
                prediction[i] = 1;
            }
            else{
                prediction[i] = 0;
            }
        }

    }


    return prediction;
}


auto MERIDA_ILP::calculateObjectiveFunction(GeneTrail::DenseMatrixSubset& mtx_response, std::vector<int>& prediction) -> double{


    double objf_value{0.0};
    //Calculate objective function value
    for(int i = 0; i< mtx_response.rows(); i++){
     
        if(doubleEpsilonEqual(mtx_response(i, mtx_response.cols()-1), 1.0)){ //sensitive samples
            objf_value += prediction[i] * mtx_response(i, mtx_response.cols()-2)- (1-prediction[i]) * mtx_response(i, mtx_response.cols()-2); // (TP - FN)/n_sens
        }
        else{//resistant samples
            objf_value += (1- prediction[i])* mtx_response(i, mtx_response.cols()-2) - (prediction[i]) * mtx_response(i, mtx_response.cols()-2); // (TN - FP)/n_res
        }
    }

    return objf_value;

}

auto MERIDA_ILP::print_y_k(std::ofstream & output) -> void{

    std::vector<int> prediction(std::move(predict(data_mtx)));
    output << "Prediction" << "\n";

    for(int i = 0; i< data_mtx.rows(); i++){
        output << "sample_" << i << " " << data_mtx.rowName(i) << " " << prediction[i] << std::endl;
    }
}


