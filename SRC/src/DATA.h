
#ifndef __DATA_H__
#define __DATA_H__
#include <string>
#include <vector>
#include "CLogging.h"


struct DATA {
private:
    std::string source; // source of data: mock_wave_front vs camera_stream
    CLogging* LDATA;
    
    
public:
    int nbr; // TODO: what is this number
    std::vector<double> x, y, z, dx, dy;
    DATA(std::string _source, int _nbr):source(_source),nbr(_nbr) {

        LDATA =new CLogging((char*)"Log4cxxConfig.xml","DATA");
		LDATA->INFO((char *)"Initializing DATA");
    }
    ~DATA(){ delete LDATA;}
    
    //
    std::string getSource() {
        return source;
    }
    
};
#endif

