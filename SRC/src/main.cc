// -----------------------------------------------------------------------------
#include <string>
#include <memory> // unique
#include <iomanip>
#include <iostream>
//#include <exception>
//#include <stdexcept>
#include <sstream>
// -----------------------------------------------------------------------------
#include "CLogging.h"
#include "CError.h"
#include "DATA.h"
#include "WaveFrontReconstructor.h"
#include "MockWavefrontGenerator.h"
#include "CFactory.h"
#include "CParser.h"
#include "CParsening.h"
#include "AuxiliaryTemporaryFunctions.h"
// -----------------------------------------------------------------------------
using namespace std;
CLogging MAIN((char*)"Log4cxxConfig.xml","MAIN");
//
//int main(int argc, char** argv) {
int main() {
    // TODO: CParsening parser;
    // TODO: Decouple Mock from the contructor of WaveFrontReconstructor.
    // The MockWavefrontGenerator is an independent class, which
    // should be chosen by the user, and it may be replaced
    // in runtime by the sensor's slopes retrieving function.
    //bool generate_only=true;
    int ret =-1;
    try {
        CParsening parser("config.cfg"); // TODO: move class to static functions?
        
        // Building data -------------------------------------------------------
        unique_ptr<DATA> input_data(new DATA(parser.get_mock_type_wavefront()));
        //unique_ptr<DATA> input_data(new DATA(parser.get_mock_type_wavefront(), 30));
        // OR:
        //std::vector<double> x, y, z, dx, dy;
        //unique_ptr<DATA> input_data(new DATA(x, y, z, dx, dz ));
        // Buiding reconstructor -----------------------------------------------
        // MockWavefrontGenerator should return only the values in case DATA is a user input 
        auto reconstructor = CFactory::make(parser.get_pl_basis(), parser.get_pl_order(), new MockWavefrontGenerator(std::move(input_data)));
        // ---------------------------------------------------------------------
        // Testing zone (used for checking numerical procedure, please
        // do not erase until numerical algorithms are fixed.
        MAIN.INFO((char *)"Computing image reconstruction...");
        // TODO: compute instead of ComputeReconstructedWaveFront
        // TODO: why do you need parsing the wf twice? here and in the factory#make?
        vector<double> zReconstructed = reconstructor->ComputeReconstructedWaveFront(reconstructor->data->x,reconstructor->data->y,reconstructor->_coeffs);
        cout << scientific << setprecision(10) << endl;
        // TODO: encapsulate all this stuff
        cout << "NumberOfNodes: " << reconstructor->data->dx.size() << endl;
        cout << "CPUTimeSVD: " << reconstructor->GetCPUTimeMatrixGeneration() << endl;
        cout << "CPUTimePureSVD: " << reconstructor->GetCPUTimePureSVD() << endl;
        cout << "CPUTimeGenMatrixM: " << reconstructor->GetCPUTimeGenerationOfMatrixM() << endl;
        cout << "CPUTimeGenMatrixR: " << reconstructor->GetCPUTimeGenerationOfMatrixR() << endl;
        cout << "SingleCPUTimeVSUGProduct: " << reconstructor->GetCPUTimeCoefficientEstimation() << endl;
        cout << "SingleCPUTimeRAProduct: " << reconstructor->GetCPUTimeSimpleReconstruction() << endl;
        // ---------------------------------------------------------------------
        // TODO: twice the same call??
        reconstructor->CenterWavefrontAlongZ(reconstructor->data->z);
        reconstructor->CenterWavefrontAlongZ(zReconstructed);
        // ---------------------------------------------------------------------
        AuxiliaryTemporaryFunctions::SaveData3D(reconstructor->data->x, reconstructor->data->y,reconstructor->data->z,string("generated.tsv"));
        AuxiliaryTemporaryFunctions::SaveData3D(reconstructor->data->x,reconstructor->data->y,zReconstructed,string("reconstructed.tsv"));
        AuxiliaryTemporaryFunctions::SaveComparedData3D(reconstructor->data->x, reconstructor->data->y,reconstructor->data->z,zReconstructed,string("difference.tsv"));
        AuxiliaryTemporaryFunctions::SaveDiffData2D(reconstructor->data->x, reconstructor->data->y,reconstructor->data->z,zReconstructed,string("error.dat"));
        // ---------------------------------------------------------------------
        MAIN.DEB((char *)"Finished Saving...");
        double sumtot,sumres;
        // TODO: the scientific has to be encapsulted in reconstructor + encapsulte + parse as an option
        cout << std::scientific << std::setprecision(12);
        cout << "CoefficientOfDetermination: " << AuxiliaryTemporaryFunctions::GetCoefficientOfDetermination(reconstructor->data->z,zReconstructed,sumtot,sumres) << endl;
        cout << "NormalizedRMS: " << AuxiliaryTemporaryFunctions::GetNormalizedRMS(reconstructor->data->z,zReconstructed) << endl;
        cout << "PolynomialType: " << reconstructor->PolynomialType() << endl;
        cout << "PolynomialOrder: " << reconstructor->PolynomialOrder() << endl;
        cout << "NumberOfPolynomialTerms: " << reconstructor->NumberOfPolynomialTerms() << endl;
        // TODO: encapsulate in higher function + parse as an option
        AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string("wf"),string("generated"),string("generated.tsv"),string("generated.pdf"));
        AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string("wf"),string("reconstructed"),string("reconstructed.tsv"),string("reconstructed.pdf"));
        AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string("wf"),string("difference"),string("difference.tsv"),string("difference.pdf"));
        AuxiliaryTemporaryFunctions::MakePlot2DGnuplot(string("error [perc.]"),string("error"),string("error.dat"),string("error.pdf"));
        MAIN.DEB((char *)"Finished plotting.");
        MAIN.DEB((char *)"Starting RA and VSUG products loop...");
        double cpuRATime=0.0e0,cpuVSUGTime=0.0e0;
        int III=10;
        for ( int i=0 ; i<III ; ++i ) {
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
