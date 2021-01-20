//
// Created by klenhof on 19.08.19.
//

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#ifndef Pathway_Pathway_H
#define Pathway_Pathway_H

class Pathway{

    std::string pathway_definition_filename;
    std::map<std::string, std::vector<std::string>> pw_name_to_gene;

public:

    Pathway(std::string pathway_filename);
    auto getNumberOfPathways() -> size_t;
    auto getPathwaysOfGene(std::string gene_name) -> std::vector<std::string>;
    auto getPathways() -> std::vector<std::string>; //is sorted, hence order is preserved
    auto getPathwayElements(std::string pathway_name) -> std::vector<std::string>;
    auto isElementOfPathway(std::string pathway_name, std::string gene_name) -> bool;
    auto getPathwayDefinitionFilename() -> std::string;

};

#endif //Pathway_Pathway_H