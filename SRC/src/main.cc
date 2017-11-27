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
#include "AuxiliaryTemporaryFunctions.h"
// -----------------------------------------------------------------------------
using namespace std;
CLogging MAIN((char*)"Log4cxxConfig.xml","MAIN");
//
//int main(int argc, char** argv) {
int main() {
    
    // TODO: CParsening parser;  
    bool exists;
    // TODO: Decouple Mock from the contructor of WaveFrontReconstructor.
    // The MockWavefrontGenerator is an independent class, which
    // should be chosen by the user, and it may be replaced
    // in runtime by the sensor's slopes retrieving function.
    //bool generate_only=true;
    string ttt;
    int pl_order;
    string pl_basis;
    string mock_wave_front;
    int ret =-1;

    try {
        // Parsening data ------------------------------------------------------
        // TODO: ifnoexist file
        ConfigFile cfg("config.cfg");
       
        exists = cfg.keyExists("libtype");
        if(!exists)  error_arguments("libtype is missing");
        std::string lib_type= cfg.getValueOfKey<std::string>("libtype", "armadillo");
       
        //exists = cfg.keyExists("generate_only");
        //if(!exists) {  error_arguments("genarate_only is missing"); generate_only = false; }
        //ttt=cfg.getValueOfKey<std::string>("generate_only","false");
        
        exists = cfg.keyExists("polynomialbasis");
        if(!exists)  error_arguments("polynomialbasis is missing");
        pl_basis= cfg.getValueOfKey<std::string>("polynomialbasis", "Zernike"); 
        
        exists = cfg.keyExists("polynomialorder");
        if(!exists)  error_arguments("polynomialorder is missing");
        pl_order=cfg.getValueOfKey<int>("polynomialorder",10);
        
        exists = cfg.keyExists("mock_wave_front");
        if(!exists)  error_arguments("mock_wave_frontis missing");
        mock_wave_front=cfg.getValueOfKey<string>("mock_wave_front","CircularTestF2");

        // Building data -------------------------------------------------------
        unique_ptr<DATA> input_data(new DATA(mock_wave_front, 30));

        // Buiding reconstructor -----------------------------------------------
        auto reconstructor = CFactory::make(pl_basis, pl_order, new MockWavefrontGenerator(std::move(input_data)));

		// TODO: vector as a return value
		vector<double> zReconstructed;
		// ---------------------------------------------------------------------
		// Testing zone (used for checking numerical procedure, please
		// do not erase until numerical algorithms are fixed.
		MAIN.INFO((char *)"Computing image reconstruction...");
		// TODO: compute instead of ComputeReconstructedWaveFront
		// TODO: why do you need parsing the wf twice? here and in the factory#make?
		reconstructor->ComputeReconstructedWaveFront(reconstructor->data->x,reconstructor->data->y,reconstructor->_coeffs,zReconstructed);
		cout << scientific << setprecision(10) << endl;
		cout << "NumberOfNodes: " << reconstructor->data->dx.size() << endl;
		cout << "CPUTimeSVD: " << reconstructor->GetCPUTimeMatrixGeneration() << endl;
        cout << "CPUTimePureSVD: " << reconstructor->GetCPUTimePureSVD() << endl;
        cout << "CPUTimeGenMatrixM: " << reconstructor->GetCPUTimeGenerationOfMatrixM() << endl;
        cout << "CPUTimeGenMatrixR: " << reconstructor->GetCPUTimeGenerationOfMatrixR() << endl;
		cout << "SingleCPUTimeVSUGProduct: " << reconstructor->GetCPUTimeCoefficientEstimation() << endl;
		cout << "SingleCPUTimeRAProduct: " << reconstructor->GetCPUTimeSimpleReconstruction() << endl;
        // ---------------------------------------------------------------------
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
		cout << std::scientific << std::setprecision(12);
		cout << "CoefficientOfDetermination: " << AuxiliaryTemporaryFunctions::GetCoefficientOfDetermination(
				reconstructor->data->z,zReconstructed,sumtot,sumres) << endl;
        cout << "NormalizedRMS: " << AuxiliaryTemporaryFunctions::GetNormalizedRMS(reconstructor->data->z,zReconstructed) << endl;
		cout << "PolynomialType: " << reconstructor->PolynomialType() << endl;
		cout << "PolynomialOrder: " << reconstructor->PolynomialOrder() << endl;
		cout << "NumberOfPolynomialTerms: " << reconstructor->NumberOfPolynomialTerms() << endl;

		AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string("wf"),string("generated"),string("generated.tsv"),string("generated.pdf"));
		AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string("wf"),string("reconstructed"),string("reconstructed.tsv"),
				string("reconstructed.pdf"));
		AuxiliaryTemporaryFunctions::MakePlot3DGnuplot(string("wf"),string("difference"),string("difference.tsv"),
				string("difference.pdf"));
		AuxiliaryTemporaryFunctions::MakePlot2DGnuplot(string("error [perc.]"),string("error"),string("error.dat"),string("error.pdf"));
		MAIN.DEB((char *)"Finished plotting.");

		MAIN.DEB((char *)"Starting RA and VSUG products loop...");
		double cpuRATime=0.0e0,cpuVSUGTime=0.0e0;
		int III=10;
		for ( int i=0 ; i<III ; ++i ) {
			reconstructor->SetSlopes(reconstructor->data->dx,reconstructor->data->dy);
			reconstructor->ComputeCoefficients();
			reconstructor->ComputeReconstructedWaveFront(reconstructor->data->x,reconstructor->data->y,reconstructor->_coeffs,zReconstructed);
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

