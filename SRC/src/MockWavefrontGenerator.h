
#ifndef MOCKWAVEFRONTGENERATOR_H_
#define MOCKWAVEFRONTGENERATOR_H_
#include <memory>
#include <vector>
using std::vector;
#include <string>
#include <iostream>
#include "DATA.h" 
using namespace std;

#include "DATA.h"

class MockWavefrontGenerator {

    public:

       unique_ptr<DATA> data;

    protected:

        bool onSquare;
        double rMin,rMax;

    public:
       // ======================================================================
       MockWavefrontGenerator(unique_ptr<DATA> _data);
       // ======================================================================
       // If the grid is nxn, this function returns n.  
       int Size(void);
       // If false, it only saves the coordinates of points:
       // { (x,y) | rMin <= x^2+y^2 <= rMax  } 
       void OnSquare(bool os) { onSquare=os; }
       // ======================================================================
       /** Generates a simple tilted wave front with slopes mx and my.  */
       void GenerateXYTiltWaveFront(int n,double mx,double my);
       void GenerateXTiltWaveFront(int n,double mx) {return GenerateXYTiltWaveFront(n,mx,0.0e0);}
       void GenerateYTiltWaveFront(int n,double my) {return GenerateXYTiltWaveFront(n,0.0e0,my);}
       /** Returns a Gaussian-shaped wavefront, centered at (x0,y0) with width=width and
        * height=height. n is the number of points per direction of the grid (nxn).  */
       void GenerateGaussianWaveFront(int n, double x0,double y0,double width,double height=1.0e0);
       /** Generates a wavefront that has two off-centered Gaussian functions. Gaussian
        * 1[2] is centered at (x1,y1)[(x2,y2)], and have widht=width1[width2] and 
        * height=height1[height2]. n is the number of points per direction of the
        * grid (nxn), and it only produces points for which \f$\sqrt{x^2+y^2}<1\f$  */
       void GenerateCircularDoubleOffCenteredGaussians(int n,double x1,double y1,double width1,
               double height1,double x2,double y2,double width2,double height2);
       /** Generates a Super-Gaussian-shaped wavefront, centered at (x0,y0), etc. The order
         of the superGaussian function is n (a Gaussian function is a super Gaussian function
         of order 2).  */
       void GenerateSuperGaussianWaveFront(int n, double x0,double y0,int ord,double width,\
             double height=1.0e0);
       /** Generates z=height*exp(-(x^{2N}/width^{2N}+y^{2N}/width^{2N})), and  */
       void GenerateSquaredSuperGaussianWaveFront(int n, double x0,double y0,int ord,double width,\
             double height=1.0e0);
       // TODO: what is n?
       void GenerateTestF1WaveFront(int n);
       void GenerateTestF2WaveFront(int n);
       // ======================================================================
       void CenterWavefrontAlongZ(void);
       void Save(const char *fileName);
       void Save(std::string fileName) {Save(fileName.c_str());}
       // ======================================================================
       void SetRMin(double rr) { rMin=rr; }
       void SetRMax(double rr) { rMax=rr; }
       void Print(void);
       // ======================================================================

    protected:
       void ComputeCoordinatesCircle(int n);
       void ComputeCoordinatesSquare(int n);
       void ComputeXYTiltWaveFront(double mx,double my);
       void ComputeGaussianWaveFront(double x0, double y0, double width,double height);
       void ComputeDoubleOffCenteredGaussiansWaveFront(double x01,double y01,\
               double x02,double y02,double width1,double height1,\
               double width2,double height2);
       void ComputeTestF1WaveFront(void);
       void ComputeTestF2WaveFront(void);
       void ComputeSuperGaussianWaveFront(double x0,double y0,int order,\
             double width,double height);
       void ComputeSquaredSuperGaussianWaveFront(double x0,double y0,int order,\
             double width,double height);
       void SetupCoordinates(int n);
};


#endif  /* MOCKWAVEFRONTGENERATOR_H_ */

