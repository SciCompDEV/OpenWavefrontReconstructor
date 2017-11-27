#include "CLogging.h"


CLogging::CLogging(string _cfgFile, string _ptrLoggerName):
    cfgFile(_cfgFile),
    ptrLoggerName(_ptrLoggerName) {

    if (fileExists(cfgFile)) {

        try {
            DOMConfigurator::configure(cfgFile);
        } catch (log4cxx::helpers::Exception & e) {
            std::cerr << "Error : log4cxx::DOMConfigurator:configure" << std::endl;
        }

    } else {

        try {
            // Set up a simple configuration that logs on the console.
            std::cout << "Log4Cxx: BasicConfigurator (cfg file not found, using console for logging)" << std::endl;
            BasicConfigurator::configure();
        } catch (log4cxx::helpers::Exception & e) {
            std::cerr << "Error : log4cxx::BasicConfigurator" << std::endl;
        }
    }

    ptrLogger=Logger::getLogger(ptrLoggerName);

}




void CLogging::TRACE(char* format, ...) {
    va_list args;
    va_start(args, format);
    vsprintf(buffer,format, args);
    str=buffer;
    va_end (args);  
    TRACE_(ptrLogger,str);
}

void CLogging::DEB(char* format, ...) {
    va_list args;
    va_start(args, format);
    vsprintf(buffer,format, args);
    str=buffer;
    va_end (args);  
    DEB_(ptrLogger,str);
}

void CLogging::INFO(char* format, ...) {
    va_list args;
    va_start(args, format);
    vsprintf(buffer,format, args);
    str=buffer;
    va_end (args);  
    INFO_(ptrLogger,str);
}

void CLogging::WARN(char* format, ...) {
    va_list args;
    va_start(args, format);
    vsprintf(buffer,format, args);
    str=buffer;
    va_end (args);  
    WARN_(ptrLogger,str);
}

void CLogging::ERROR(char* format, ...) {
    va_list args;
    va_start(args, format);
    vsprintf(buffer,format, args);
    str=buffer;
    va_end (args);  
    ERROR_(ptrLogger,str);
}

void CLogging::FATAL(char* format, ...) {
    va_list args;
    va_start(args, format);
    vsprintf(buffer,format, args);
    str=buffer;
    va_end (args);  
    FATAL_(ptrLogger,str);
}
