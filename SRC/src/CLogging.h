


#ifndef __CLOGGING_H_
#define __CLOGGING_H_

//hola
#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>
#include <log4cxx/helpers/exception.h>
#include <log4cxx/basicconfigurator.h>

#include <iostream>
#include <string>
#include <cstdarg>
#include <cstdio>
#include <sys/stat.h>

using namespace std;

using namespace log4cxx;
using namespace log4cxx::xml;
using namespace log4cxx::helpers;

#define TRACE_ LOG4CXX_TRACE
#define DEB_   LOG4CXX_DEBUG
#define INFO_  LOG4CXX_INFO
#define WARN_  LOG4CXX_WARN
#define ERROR_ LOG4CXX_ERROR
#define FATAL_ LOG4CXX_FATAL

class CLogging {

private:

    string cfgFile;
    string ptrLoggerName;

    LoggerPtr ptrLogger;

    char buffer[200];
    string str;

    bool fileExists(const std::string& filename) {
        struct stat buf;
        if (stat(filename.c_str(), &buf) != -1) {
            return true;
        }
        return false;
    }

public:

    CLogging(string , string );

    void TRACE(char* fmt ,...);
    void DEB(char* fmt ,...);
    void INFO(char* fmt ,...);
    void WARN(char* fmt ,...);
    void ERROR(char* fmt ,...);
    void FATAL(char* fmt ,...);

};

#endif
