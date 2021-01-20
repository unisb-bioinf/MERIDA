
//
// Created by klenhof on 19.08.19.
//

#include <Pathway.h>

Pathway::Pathway(std::string pathway_filename): pathway_definition_filename(pathway_filename){

    std::ifstream pw_is(pathway_filename);

    std::string line;

    if(!pw_is){
        std::cout << "Pathway definition filename was invalid." << std::endl;
    }
    else{
        while(getline(pw_is, line)){

            std::istringstream liss(line);
            std::string element;

            int counter(-1);
            std::string pw_name;
            std::vector<std::string> genes;
            while(std::getline(liss, element, '\t')){

                counter++;
                if(counter == 0){

                    pw_name = element;
                }
                else if(counter == 1){
                    continue;
                }

                else{
                    genes.push_back(element);
                }

            }

            if(!pw_name.empty()){
                pw_name_to_gene.insert(std::pair<std::string, std::vector<std::string>>(pw_name, genes));
            }

        }
    }


    pw_is.close();

}

auto Pathway::getNumberOfPathways() -> size_t {

    return pw_name_to_gene.size();
}

auto Pathway::getPathwaysOfGene(std::string gene_name) -> std::vector<std::string> {

    std::vector<std::string> pathways;

    for(auto &it : pw_name_to_gene){

        std::vector<std::string> genes(it.second);
        for(size_t i = 0; i< genes.size(); i++){
            if(boost::iequals(genes[i], gene_name)){//should be case insensitive
                pathways.push_back(it.first);
                break;
            }
        }

    }

    return pathways;

}

auto Pathway::getPathways() -> std::vector<std::string>{

    std::vector<std::string> pathways;
    for(auto &it : pw_name_to_gene){

        pathways.push_back(it.first);
    }

    std::sort(pathways.begin(), pathways.end());

    return pathways;
}

auto Pathway::getPathwayElements(std::string pathway_name) -> std::vector<std::string>{

    try{
        return pw_name_to_gene.at(pathway_name);
    }
    catch(const std::out_of_range& not_available){

        return std::vector<std::string>();
    }
}

auto Pathway::getPathwayDefinitionFilename() -> std::string{
    return boost::filesystem::path(pathway_definition_filename).stem().string();
}

auto Pathway::isElementOfPathway(std::string pathway_name, std::string gene_name) -> bool{


    for (size_t i = 0; i < pw_name_to_gene.at(pathway_name).size(); i++) {

        if (boost::iequals(pw_name_to_gene.at(pathway_name)[i], gene_name)) {
            return true;
        }
    }

    return false;


}