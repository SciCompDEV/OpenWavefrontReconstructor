#ifndef __CFactory_H_
#define __CFactory_H_
#include "WaveFrontReconstructor.h"
#include "MockWavefrontGenerator.h"
#include "ZernikeReconstructor.h"
#include "LegendreReconstructor.h"
#include "HalfCircularHarmonicsReconstructor.h"
#include "CError.h"
//
#include "CLogging.h"
//
#include <memory>
#include <iostream>
#include <vector>
// -----------------------------------------------------------------------------
// make_unique not yet implemented in C++11
using namespace std;
	template<typename T, typename... Args>
	std::unique_ptr<T> make_unique(Args&&... args) {
		return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
	}

// -----------------------------------------------------------------------------
CLogging CFACTORY((char*)"Log4cxxConfig.xml","CFACTORY");
class CFactory {
public:
    // TODO: pl-... wrapper && typedef
    static unique_ptr<WavefrontReconstructor> make(string pl_type, int pl_order, MockWavefrontGenerator* wf) {
        CFACTORY.INFO((char *)"Trying to build up reconstructor type: %s", const_cast<char*>(pl_type.c_str()));
		unique_ptr<WavefrontReconstructor> product; 
        // ---------------------------------------------------------------------
        // TODO: put a map for selection here
		if( pl_type.compare("Zernike")==0 ) {
			product = make_unique<ZernikeReconstructor>(pl_order,(wf->data).release());
		} else if ( pl_type.compare("HalfCircularHarmonics")==0 ) {
			product = make_unique<HalfCircularHarmonicsReconstructor>(pl_order,(wf->data).release());
		} else if( pl_type.compare("Legendre")==0 ) {
			product = make_unique<LegendreReconstructor>(pl_order,(wf->data).release());
		}  else {
			error_parameters("Polynomial basis have NOT been provided?");
		}
        // ---------------------------------------------------------------------
        if(product == nullptr) {
            error_undefined("Reconstructor cannot be created");
        } else {
            return std::move(product); // copy elision?
        }
    } // END WaveFrontReconstructor
};
#endif

