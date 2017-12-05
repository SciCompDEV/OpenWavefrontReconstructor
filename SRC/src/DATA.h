#ifndef __DATA_H__
#define __DATA_H__
#include <string>
#include <vector>
#include "CLogging.h"
//#include "UTILs.h"
struct DATA {
private:
    CLogging* LDATA;
    std::string input_field; // input_field name of data: mock_wave_front vs external_stream
public:
    // User input field
    std::vector<double> x,
        y,
        z,
        dx,
        dy;
public:
    DATA(std::string _input_field): 
        input_field(_input_field) {
        LDATA = new CLogging((char*)"Log4cxxConfig.xml","DATA");
        LDATA->INFO((char *)"Initializing wavefront DATA from mockup");
    }
    DATA(std::vector<double>& _x,
         std::vector<double>& _y,
         std::vector<double>& _z,
         std::vector<double>& _dx,
         std::vector<double>& _dy)  {
        // getter/setter
        input_field = "user";
        LDATA = new CLogging((char*)"Log4cxxConfig.xml","DATA");
        LDATA->INFO((char *)"Initializing wafefront DATA from external input_field");
        // TODO: autocheck/sanitize
        x=_x;
        y=_y;
        z=_z;
        dx=_dx;
        dy=_dy;
    }
    ~DATA() {
        delete LDATA;
    }
    //
    std::string getSource() {
        return input_field;
    }
};
#endif
