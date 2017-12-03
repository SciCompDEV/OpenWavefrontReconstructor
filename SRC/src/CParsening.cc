#include "CParsening.h" 



CParsening::CParsening(std::string file) {
    cfg  = new ConfigFile(file);
    
};
CParsening::~CParsening() {

    delete cfg;

}


std::string CParsening::get_lib_type() {

    exists = cfg->keyExists("libtype");
    if(!exists)  error_arguments("libtype is missing");
    return  cfg->getValueOfKey<std::string>("libtype", "armadillo");

}

std::string CParsening::get_pl_basis() {

    exists = cfg->keyExists("polynomialbasis");
    if(!exists)  error_arguments("polynomialbasis is missing");
    return cfg->getValueOfKey<std::string>("polynomialbasis", "Zernike");

}

int  CParsening::get_pl_order() {

        exists = cfg->keyExists("polynomialorder");
        if(!exists)  error_arguments("polynomialorder is missing");
        return cfg->getValueOfKey<int>("polynomialorder",10);

}
std::string CParsening::get_mock_type_wavefront() {

        exists = cfg->keyExists("mock_type_wavefront");
        if(!exists)  error_arguments("mock_type_wavefrontis missing");
        return cfg->getValueOfKey<string>("mock_type_wavefront","CircularTestF2");

}


