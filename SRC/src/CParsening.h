


#ifndef __CPARSENING_H__
#define __CPARSENING_H__

#include "CParser.h"
#include "CError.h"
using namespace std;

class CParsening
{
private:
    bool exists;
    ConfigFile* cfg;

public:
    CParsening(std::string file);
    virtual ~CParsening();
    // 
    std::string get_lib_type();
    std::string get_pl_basis();
    int get_pl_order();
    std::string get_mock_type_wavefront();
};

#endif 



