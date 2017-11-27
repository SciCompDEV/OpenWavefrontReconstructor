
#include "CParser.h"



        // Constructor
ConfigFile::ConfigFile(const std::string &fName) {
    this->fName = fName;
    ExtractKeys();
}

