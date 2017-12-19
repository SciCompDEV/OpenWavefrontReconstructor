// -----------------------------------------------------------------------------
#include <string>
#include <memory> 
#include <iomanip>
#include <iostream>
#include <sstream>
// -----------------------------------------------------------------------------
#include "CLogging.h"
#include "CError.h"
#include "DATA.h"
#include "WaveFrontReconstructor.h"
#include "MockWavefrontGenerator.h"
#include "CFactory.h"
#include "CParsening.h"
#include "AuxiliaryTemporaryFunctions.h"
// -----------------------------------------------------------------------------
using namespace std; // TODO
CLogging MAIN((char*)"Log4cxxConfig.xml","MAIN");
//
int main() {
    int ret =-1;
    try {
        CParsening parser("config.cfg"); // TODO: move class to static functions?
        
        unique_ptr<DATA> input_data(new DATA(parser.get_mock_type_wavefront())); // TODO: confirm the 30 is not neccesary anymore
        // OR:
        //std::vector<double> x, y, z, dx, dy;
        //unique_ptr<DATA> input_data(new DATA(x, y, z, dx, dz )); // TODO: check if  thi sis correct:in runtime by the sensor's slopes retrieving function.
        
        auto reconstructor = CFactory::make(parser.get_pl_basis(), parser.get_pl_order(), new MockWavefrontGenerator(std::move(input_data)));

        MAIN.INFO((char *)"Computing image reconstruction...");
        // TODO: compute instead of ComputeReconstructedWaveFront
        // TODO: why do you need parsing the wf twice? here and in the factory#make?
        vector<double> zReconstructed = reconstructor->ComputeReconstructedWaveFront(reconstructor->data->x,reconstructor->data->y,reconstructor->_coeffs);

        reconstructor->print_CPU_gathered_info();

        // TODO: Solano, i want to move the args to a return value, why are calling this function twice?
        reconstructor->CenterWavefrontAlongZ(reconstructor->data->z);
        reconstructor->CenterWavefrontAlongZ(zReconstructed); // TODO: Move this as a return values

        AuxiliaryTemporaryFunctions::SaveData3D_gathered(reconstructor->data->x, 
                reconstructor->data->y, 
                reconstructor->data->z, 
                zReconstructed, "generated.tsv", 
                "reconstructed.tsv", 
                "difference.tsv", 
                "error.dat"); 

        MAIN.DEB((char *)"Finished Saving...");


        // TODO: Solano sumtot, and sumres are not used!
        double sumtot,sumres;
        AuxiliaryTemporaryFunctions::print_simu_param_gather(sumtot, 
                sumres, 
                reconstructor->data->z, 
                zReconstructed,
                reconstructor->PolynomialType(),
                reconstructor->PolynomialOrder(),
                reconstructor->NumberOfPolynomialTerms() );

        // TODO: encapsulate in higher function + parse as an option
        AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string("wf"),string("generated"),string("generated.tsv"),string("generated.pdf"));
        AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string("wf"),string("reconstructed"),string("reconstructed.tsv"),string("reconstructed.pdf"));
        AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string("wf"),string("difference"),string("difference.tsv"),string("difference.pdf"));
        AuxiliaryTemporaryFunctions::MakePlot2DGnuplot(string("error [perc.]"),string("error"),string("error.dat"),string("error.pdf"));
        MAIN.DEB((char *)"Finished plotting.");
        
        MAIN.DEB((char *)"Starting RA and VSUG products loop...");
        double cpuRATime=0.0e0,cpuVSUGTime=0.0e0;
        int III=10; // TODO why 10??
        for ( int i=0 ; i<III ; ++i ) {
            // TODO Solano WavefrontReconstructor should have a DATA object.There is no point at passing reconstructor->data->dx
            // to the reconstructor itself
            reconstructor->SetSlopes(reconstructor->data->dx,reconstructor->data->dy);
            reconstructor->ComputeCoefficients();
            zReconstructed = reconstructor->ComputeReconstructedWaveFront(reconstructor->data->x,reconstructor->data->y,reconstructor->_coeffs);
            cpuVSUGTime+=reconstructor->GetCPUTimeCoefficientEstimation();
            cpuRATime+=reconstructor->GetCPUTimeSimpleReconstruction();
        }
        cout << "AverageCPUTimeVSUGProduct: " << (cpuVSUGTime/double(III)) << endl;
        cout << "AverageCPUTimeRAProduct: " << (cpuRATime/double(III)) << endl;
        MAIN.DEB((char *)"RA and VSUG products loop finished...");
        ret = 1;
    } catch (std::exception & e) {
        std::cout << e.what() << std::endl;
    }
    if(ret==1) MAIN.INFO((char*)"DONE!");
    return ret;
}
