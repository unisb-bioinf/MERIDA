//
// Created by klenhof on 08.03.18.
//

#ifndef ACT_LOMO_CROSSVALIDATOR_H
#define ACT_LOMO_CROSSVALIDATOR_H
#include <string>

/*Actually wanted a pure virtual function with a template argument, needed to use CRTP instead.*/
template<class T>
class CrossValidator {
protected:
    std::string paramfile;
    int foldsize;
public:
    CrossValidator(std::string paramfile, int foldsize): paramfile(paramfile), foldsize(foldsize){}
    /*This function needs to be implemented by a derived class*/
    template <typename Logmod>
    auto crossvalidate() -> void {
        static_cast<T*>(this)->crossvalidate();
    }
};


#endif //ACT_LOMO_CROSSVALIDATOR_H
