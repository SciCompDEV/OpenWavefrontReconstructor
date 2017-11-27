
#ifndef __CERROR
#define __CERROR

#include <string>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

// TODO: check thread safe
template<int _error_type>
class CError : public std::exception
{
  private:
     string msg;
     int error_type;
     const char * file;
     int line;

  public:

     CError(string _msg, const char *_file, int _line) : exception(),  msg(_msg), error_type(_error_type), file(_file), line(_line) {
        stringstream ss;
        ss.str("");
        ss	<< "EXCEPTION type: "<< getType() << " " << msg <<  " ["<<file<<  ","<< line << "]";
        msg=ss.str();
     }   
    const char* what() const throw (){
       return msg.c_str();
    }
    string getType() {
        switch (error_type) {
            case 1:
                return "PARAMETERS";
            case 2:
                return "ARGUMENTS";
            case 3:
                return "UNDEFINED";
            default:
                return "UNDEFINED";
        }
    
    }

};


#define error_parameters(arg) throw CError<1>(arg, __FILE__, __LINE__);
#define error_arguments(arg) throw CError<2>(arg, __FILE__, __LINE__);
#define error_undefined(arg) throw CError<3>(arg, __FILE__, __LINE__);

#endif



