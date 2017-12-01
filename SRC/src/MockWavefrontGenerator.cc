#ifndef _MOCKWAVEFRONTGENERATOR_CPP_
#define _MOCKWAVEFRONTGENERATOR_CPP_
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <iomanip>
#include <string>
using std::string;
using std::to_string;
#include <fstream>
using std::ofstream;
#include <cmath>
#include "MockWavefrontGenerator.h"
#define MIN_NON_ZERO_EPS_FOR_R (1.0e-06)
//Test-only class. Removed/replaced in user applications.
// CONSTRUCTOR =================================================================
MockWavefrontGenerator::MockWavefrontGenerator(unique_ptr<DATA> _data) {
    data=move(_data);
    data->x.clear();
    data->y.clear();
    data->z.clear();
    data->dx.clear();
    data->dy.clear();
    onSquare=false;
    rMin=0.0e0;
    rMax=1.0e0;
    int gridDim=30; //here 30 is used to form the 30x30 grid from which it selects points inside a circle
                    
    if(data->getSource().compare("CircularTestF1")==0) {
        cout << "Setting up Test F1 wavefront, circular." << endl;
        GenerateTestF1WaveFront(gridDim);
    } else if ( data->getSource().compare("CircularXYTilted")==0 ) {
        cout << "Setting up Tilted (on X and Y) wavefront, circular." << endl;
        GenerateXYTiltWaveFront(gridDim,0.5e0,0.2e0);
    } else if ( data->getSource().compare("CircularCenteredGaussian")==0 ) {
        cout << "Setting up a Centered Gaussian wavefront, circular." << endl;
        GenerateGaussianWaveFront(gridDim,0.0e0,0.0e0,0.5e0);
    } else if ( data->getSource().compare("CircularOffCenteredGaussian")==0 ) {
        cout << "Setting up an Off-centered Gaussian wavefront, circular." << endl;
        GenerateGaussianWaveFront(gridDim,0.1e0,0.3e0,0.6e0,0.6);
    } else if ( data->getSource().compare("CircularDoubleOffCenteredGaussians")==0 ) {
        cout << "Setting up Two off-centered Gaussians wavefront, circular." << endl;
        GenerateCircularDoubleOffCenteredGaussians(gridDim,0.3e0,0.5e0,0.25e0,0.5e0,
                -0.5e0,-0.3e0,0.10e0,0.3e0);
    } else if ( data->getSource().compare("CircularSuperGaussian4")==0 ) {
        cout << "Setting up a Super-Gaussian wavefront, circular, of order 4" << endl;
        GenerateSuperGaussianWaveFront(gridDim,0.0e0,0.0e0,4,0.8e0,0.7e0);
    } else if ( data->getSource().compare("CircularSuperGaussian6")==0 ) {
        cout << "Setting up a Super-Gaussian wavefront, circular, of order 6" << endl;
        GenerateSuperGaussianWaveFront(gridDim,0.0e0,0.0e0,6,0.6e0,0.5e0);
    } else if ( data->getSource().compare("CircularSuperGaussian8")==0 ) {
        cout << "Setting up a Super-Gaussian wavefront, circular, of order 8" << endl;
        GenerateSuperGaussianWaveFront(gridDim,0.00e0,0.00e0,8,0.6e0,0.6e0);
    } else if ( data->getSource().compare("SquareXYTilted")==0 ) {
        cout << "Setting up Tilted (on X and Y) wavefront, square." << endl;
        onSquare=true;
        rMax=0.99e0;
        GenerateXYTiltWaveFront(gridDim,0.5e0,0.2e0);
    } else if ( data->getSource().compare("SquareTestF1")==0 ) {
        cout << "Setting up a Test F1 wavefront, square." << endl;
        onSquare=true;
        rMax=0.99e0;
        GenerateTestF1WaveFront(gridDim);
    } else if ( data->getSource().compare("SquareCenteredGaussian")==0 ) {
        cout << "Setting up a Centered Gaussian wavefront, square." << endl;
        onSquare=true;
        rMax=0.99e0;
        GenerateGaussianWaveFront(gridDim,0.0e0,0.0e0,0.5e0);
    } else if ( data->getSource().compare("SquareOffCenteredGaussian")==0 ) {
        cout << "Setting up an Off-centered Gaussian wavefront, square." << endl;
        onSquare=true;
        rMax=0.99e0;
        GenerateGaussianWaveFront(gridDim,0.1e0,0.3e0,0.6e0,0.6);
    } else if ( data->getSource().compare("SquareDoubleOffCenteredGaussians")==0 ) {
        cout << "Setting up Two off-centered Gaussians wavefront, square." << endl;
        onSquare=true;
        rMax=0.99e0;
        GenerateCircularDoubleOffCenteredGaussians(gridDim,0.3e0,0.5e0,0.25e0,0.5e0,
                -0.5e0,-0.3e0,0.10e0,0.3e0);
    } else if ( data->getSource().compare("SquareSuperGaussian4")==0 ) {
        cout << "Setting up a Super-Gaussian wavefront, square, of order 4" << endl;
        onSquare=true;
        rMax=0.99e0;
        GenerateSuperGaussianWaveFront(gridDim,0.0e0,0.0e0,4,0.8e0,0.7e0);
    } else if ( data->getSource().compare("SquareSuperGaussian6")==0 ) {
        cout << "Setting up a Super-Gaussian wavefront, square, of order 6" << endl;
        onSquare=true;
        rMax=0.99e0;
        GenerateSuperGaussianWaveFront(gridDim,0.0e0,0.0e0,6,0.6e0,0.5e0);
    } else if ( data->getSource().compare("SquareSuperGaussian8")==0 ) {
        cout << "Setting up a Super-Gaussian wavefront, square, of order 8" << endl;
        onSquare=true;
        rMax=0.99e0;
        GenerateSuperGaussianWaveFront(gridDim,0.00e0,0.00e0,8,0.6e0,0.6e0);
    } else {
        cout << "Setting up the default xy tilted wavefront." << endl;
        GenerateXYTiltWaveFront(gridDim,0.5e0,0.2e0);
    }
}
// =============================================================================
// COORDINATES =================================================================
void MockWavefrontGenerator::SetupCoordinates(int n) {
    data->x.clear();
    data->y.clear();
    data->z.clear();
    data->dx.clear();
    data->dy.clear();
    if ( onSquare ) {
        ComputeCoordinatesSquare(n);
    } else {
        ComputeCoordinatesCircle(n);
    }
    data->z.resize((data->x.size()),0.0e0);
    data->dx.resize((data->x.size()),0.0e0);
    data->dy.resize((data->x.size()),0.0e0);
} // END -----------------------------------------------------------------------
void MockWavefrontGenerator::ComputeCoordinatesCircle(int n) {
    if ( n<=0 ) {
        throw (string("You should use n>0!\n" )\
               +string(__FILE__)+string(": ")+to_string(__LINE__));
        return;
    }
    double dx=2.0e0/double(n-1);
    double cur2x=-1.0e0,cur2y,r2;
    double r2Max=rMax*rMax;
    double r2Min=rMin*rMin;
    for ( int i=0 ; i<n ; ++i ) {
        cur2y=-1.0e0;
        for ( int j=0 ; j<n ; ++j ) {
            r2=cur2x*cur2x+cur2y*cur2y;
            if ( (r2<=r2Max) && (r2>=r2Min) ) {
                data->x.push_back(cur2x);
                data->y.push_back(cur2y);
            }
            cur2y+=dx;
        }
        cur2x+=dx;
    }
}
void MockWavefrontGenerator::ComputeCoordinatesSquare(int n) {
    if ( n<=0 ) {
        throw (string("You should use n>0!\n" )\
               +string(__FILE__)+string(": ")+to_string(__LINE__));
        return;
    }
    double dx=2.0e0*rMax/double(n-1);
    double cur2x=-rMax,cur2y;
    for ( int i=0 ; i<n ; ++i ) {
        cur2y=-rMax;
        for ( int j=0 ; j<n ; ++j ) {
            data->x.push_back(cur2x);
            data->y.push_back(cur2y);
            cur2y+=dx;
        }
        cur2x+=dx;
    }
}
// =============================================================================
//  MOCK GENERATORS ============================================================
//  Generates a grid of points:
/*
 1- The initial grid is set to be 30x30 points, from -1 to 1 in x and y directions.
 Then the class removes the points outside a circle of radius 1 reducing the actual
 number of points that are used which reduces the computational cost.
 The library assumes the user provide arrays with coordinates and slopes in cartesian coordinates, thus
 MockWavefrontGenerator object is not needed in the final user code.
 The class MockWavefrontGenerator is, a Mock wavefront generator, and the wavefront is known.
 In a real problem, the function (wavefront or image) is not known!
*/
void MockWavefrontGenerator::GenerateTestF2WaveFront(int n) {
    SetupCoordinates(n);
    ComputeTestF2WaveFront();
}
void MockWavefrontGenerator::GenerateGaussianWaveFront(int n, double x0,double y0,double width,double height) {
    SetupCoordinates(n);
    ComputeGaussianWaveFront(x0,y0,width,height);
}
void MockWavefrontGenerator::GenerateCircularDoubleOffCenteredGaussians(int n,double x1,double y1,double width1,
               double height1,double x2,double y2,double width2,double height2) {
    SetupCoordinates(n);
    ComputeDoubleOffCenteredGaussiansWaveFront(x1,y1,x2,y2,width1,height1,width2,height2);
}
void MockWavefrontGenerator::ComputeDoubleOffCenteredGaussiansWaveFront(double x01,double y01,\
        double x02,double y02,double width1,double height1,\
        double width2,double height2) {
    ComputeGaussianWaveFront(x01,y01,width1,height1);
    int nn=int(data->x.size());
    vector<double> dxt(nn),dyt(nn),zt(nn);//,xt(nn),yt(nn)
    for ( int i=0 ; i<nn ; ++i ) {
        //xt[i]=data->x[i];
        //yt[i]=data->y[i];
        zt[i]=data->z[i];
        dxt[i]=data->dx[i];
        dyt[i]=data->dy[i];
    }
    ComputeGaussianWaveFront(x02,y02,width2,height2);
    for ( int i=0 ; i<nn ; ++i ) {
        //data->x[i]+=xt[i];
        //data->y[i]+=yt[i];
        data->z[i]+=zt[i];
        data->dx[i]+=dxt[i];
        data->dy[i]+=dyt[i];
    }
}
void MockWavefrontGenerator::GenerateSuperGaussianWaveFront(int n, double x0,double y0,int ord,double width,double height) {
    SetupCoordinates(n);
    ComputeSuperGaussianWaveFront(x0,y0,ord,width,height);
}
void MockWavefrontGenerator::GenerateSquaredSuperGaussianWaveFront(int n, double x0,double y0,int ord,double width,double height) {
    SetupCoordinates(n);
    ComputeSquaredSuperGaussianWaveFront(x0,y0,ord,width,height);
}
void MockWavefrontGenerator::GenerateTestF1WaveFront(int n) {
    SetupCoordinates(n);
    ComputeTestF1WaveFront();
}
void MockWavefrontGenerator::GenerateXYTiltWaveFront(int n,double mx,double my) {
    SetupCoordinates(n);
    ComputeXYTiltWaveFront(mx,my);
}
// =============================================================================
// COMPUTE MOCK WAVEFRONT ======================================================
void MockWavefrontGenerator::ComputeXYTiltWaveFront(double mx,double my) {
    int nn=data->x.size();
    for ( int i=0 ; i<nn ; ++i ) {
        data->z[i]=mx*data->x[i]+my*data->y[i];
        data->dx[i]=mx;
        data->dy[i]=my;
    }
}
void MockWavefrontGenerator::ComputeGaussianWaveFront(double x0,double y0, double width,double height) {
    int nn=data->x.size();
    double r2;
    double oow2=1.0e0/(width*width);
    double twoow2=2.0e0*oow2;
    for ( int i=0 ; i<nn ; ++i ) {
        r2=((data->x[i]-x0)*(data->x[i]-x0));
        r2+=((data->y[i]-y0)*(data->y[i]-y0));
        data->z[i]=height*exp(-(r2*oow2));
        data->dx[i]=-(twoow2*(data->x[i]-x0)*data->z[i]);
        data->dy[i]=-(twoow2*(data->y[i]-y0)*data->z[i]);
    }
}
void MockWavefrontGenerator::ComputeSuperGaussianWaveFront(double x0, double y0,int order,double width,double height) {
    int nn=data->x.size();
    double r,rn,ee,row;
    for ( int i=0 ; i<nn ; ++i ) {
        r=((data->x[i]-x0)*(data->x[i]-x0));
        r+=((data->y[i]-y0)*(data->y[i]-y0));
        r=sqrt(r);
        if ( r<MIN_NON_ZERO_EPS_FOR_R ) {
            data->z[i]=height;
            data->dx[i]=0.0e0;
            data->dy[i]=0.0e0;
        } else {
            row=r/width;
            rn=1.0e0;
            for ( int j=0 ; j<order ; ++j ) {
                rn*=row;
            }
            ee=exp(-(rn*0.5e0))*height;
            r=1.0e0/r;
            r*=r; //1/r^2
            data->z[i]=ee;
            ee*=(-0.5e0*double(order)*rn*r);
            data->dx[i]=data->x[i]*ee;
            data->dy[i]=data->y[i]*ee;
        }
    }
}
void MockWavefrontGenerator::ComputeSquaredSuperGaussianWaveFront(double x0, double y0,int order,double width,double height) {
    int nn=data->x.size();
    int efford=order;
    if ( order<1 ) {
        cerr << "Non valid order for square Gaussian beam!" << endl
             << "using N=2..." << endl;
#if DEBUG
        cerr << "file: " << __FILE__ << ", line: " << __LINE__ << endl;
#endif /* ( DEBUG ) */
        efford=2;
    }
    double xx,yy,ee,xow,yow,xn,yn;
    for ( int i=0 ; i<nn ; ++i ) {
        xx=((data->x[i]-x0)*(data->x[i]-x0));
        yy=((data->y[i]-y0)*(data->y[i]-y0));
        if ( sqrt(xx+yy)<MIN_NON_ZERO_EPS_FOR_R ) {
            data->z[i]=height;
            data->dx[i]=0.0e0;
            data->dy[i]=0.0e0;
        } else {
            xow=xx/(width*width);
            xn=1.0e0;
            for ( int j=0 ; j<efford ; ++j ) {
                xn*=xow;
            }
            yow=yy/(width*width);
            yn=1.0e0;
            for ( int j=0 ; j<efford ; ++j ) {
                yn*=yow;
            }
            ee=exp(-(xn+yn))*height;
            data->z[i]=ee;
            ee*=(-double(2*efford));
            data->dx[i]=xn*ee/(data->x[i]-x0);
            data->dy[i]=yn*ee/(data->y[i]-y0);
        }
    }
}
// TODO: what the hell is that
void MockWavefrontGenerator::ComputeTestF2WaveFront() {
    // TODO: documnent
    size_t nn=data->x.size();
    double tmp1,tmp2,tmp3;
    double height=0.02e0;
    for ( size_t i=0 ; i<nn ; ++i ) {
        tmp1=data->x[i]*data->x[i]*data->x[i]; //x^3
        tmp2=data->y[i]*data->y[i];
        tmp2*=tmp2;
        tmp2*=data->y[i]; //y^5
        data->z[i]=(-data->x[i]+20.0e0*tmp1+80.0e0*tmp2);
        data->z[i]*=(12.0e0*exp(1.0e0+4.0e0*data->x[i]));
        data->z[i]-=1.0e0;
        data->z[i]*=(exp(4.0e0*data->y[i]));
        tmp3=1.0e0-2.0e0*data->x[i];
        tmp3*=tmp3; //(1-2x)^2
        data->z[i]+=(9.0e0*exp(4.0e0*data->x[i])*tmp3);
        tmp3=1.0e0+2.0e0*data->y[i];
        tmp3*=tmp3; //(1+2y)^2
        data->z[i]*=(exp(-4.0e0*data->x[i]*(1.0e0+data->x[i])-tmp3)/15.0e0);
        data->z[i]*=height;
    }
    for ( size_t i=0 ; i<nn ; ++i ) {
        tmp1=data->x[i]*data->x[i]; //x^2
        tmp2=tmp1*tmp1; //x^4
        tmp3=data->y[i]*data->y[i];
        tmp3*=tmp3;
        tmp3*=data->y[i]; //y^5
        data->dx[i]=(1.0e0-68.0e0*tmp1+160.0e0*tmp2+640.0e0*data->x[i]*tmp3);
        data->dx[i]*=(-3.0e0*exp(1.0e0+4.0e0*data->x[i]));
        data->dx[i]+=(1.0e0+2.0e0*data->x[i]);
        data->dx[i]*=(4.0e0*exp(4.0e0*data->y[i]));
        data->dx[i]-=(36.0e0*exp(4.0e0*data->x[i])*(1.0e0+8.0e0*(data->x[i]-1.0e0)*tmp1));
        tmp3=1.0e0+2.0e0*data->y[i];
        tmp3*=tmp3; //(1+2y)^2
        data->dx[i]*=(exp(-4.0e0*data->x[i]*(1.0e0+data->x[i])-tmp3)/15.0e0);
        data->dx[i]*=height;
    }
    for ( size_t i=0 ; i<nn ; ++i ) {
        tmp1=data->x[i]*data->x[i]*data->x[i]; //x^3
        tmp2=data->y[i]*data->y[i];
        tmp3=tmp2*data->y[i]; //y^3
        tmp2*=tmp3; //y^5
        data->dy[i]=(-data->x[i]+20.0e0*tmp1-50.0e0*tmp3+80.0e0*tmp2);
        data->dy[i]*=(-12.0e0*exp(1.0e0+4.0e0*data->x[i]));
        data->dy[i]+=1.0e0;
        data->dy[i]*=(8.0e0*data->y[i]*exp(4.0e0*data->y[i]));
        tmp1=1.0e0-2.0e0*data->x[i];
        tmp1*=tmp1; //(1-2x)^2
        tmp2=1.0e0+2.0e0*data->y[i];
        data->dy[i]-=(36.0e0*exp(4.0e0*data->x[i])*tmp1*tmp2);
        data->dy[i]*=(exp(-4.0e0*data->x[i]*(1.0e0+data->x[i])-tmp2*tmp2)/15.0e0);
        data->dy[i]*=height;
    }
    // TODO: move internal x,y,z to the DATA object:
    //data->x=x;
    //data->y=y;
    //data->dx=dx;
    //data->dy=dy;
    //data->z=z;
}
void MockWavefrontGenerator::ComputeTestF1WaveFront() {
    size_t nn=data->x.size();
    double tmp1,tmp2,tmp3;
    for ( size_t i=0 ; i<nn ; ++i ) {
        /*
        tmp1=1.0e0-2.0e0*data->x[i]; tmp1*=tmp1; //(1-2x)^2
        tmp2=1.0e0+2.0e0*data->y[i]; tmp2*=tmp2; //(1+2y)^2
        tmp3=4.0e0*data->x[i]*data->x[i]; //(2x)^2
        data->z[i]=3.0e0*tmp1*exp(-tmp3-tmp2);
        tmp1=tmp3; //(2x)^2
        tmp2=2.0e0*tmp1*data->x[i]; //(2x)^3
        tmp3=4.0e0*data->y[i]*data->y[i]; tmp3*=tmp3; tmp3*=(2.0e0*data->y[i]); //(2y)^5
        data->z[i]-=(10.0e0*(0.4e0*data->x[i]-tmp2-tmp3)*exp(-tmp1-4.0e0*data->y[i]*data->y[i]));
        tmp1=2.0e0*data->x[i]+1.0e0; tmp1*=tmp1; //(2x+1)^2
        tmp2=4.0e0*data->y[i]*data->y[i];
        data->z[i]-=(exp(-tmp1-tmp2)/3.0e0);
        data->z[i]*=0.2e0;
        // */
        //*
        tmp1=data->x[i]*data->x[i]*data->x[i]; //x^3
        tmp2=data->y[i]*data->y[i];
        tmp2*=tmp2;
        tmp2*=data->y[i]; //y^5
        data->z[i]=(-data->x[i]+20.0e0*tmp1+80.0e0*tmp2);
        data->z[i]*=(12.0e0*exp(1.0e0+4.0e0*data->x[i]));
        data->z[i]-=1.0e0;
        data->z[i]*=(exp(4.0e0*data->y[i]));
        tmp3=1.0e0-2.0e0*data->x[i];
        tmp3*=tmp3; //(1-2x)^2
        data->z[i]+=(9.0e0*exp(4.0e0*data->x[i])*tmp3);
        tmp3=1.0e0+2.0e0*data->y[i];
        tmp3*=tmp3; //(1+2y)^2
        data->z[i]*=(exp(-4.0e0*data->x[i]*(1.0e0+data->x[i])-tmp3)/15.0e0);
        // */
    }
    for ( size_t i=0 ; i<nn ; ++i ) {
        tmp1=data->x[i]*data->x[i]; //x^2
        tmp2=tmp1*tmp1; //x^4
        tmp3=data->y[i]*data->y[i];
        tmp3*=tmp3;
        tmp3*=data->y[i]; //y^5
        data->dx[i]=(1.0e0-68.0e0*tmp1+160.0e0*tmp2+640.0e0*data->x[i]*tmp3);
        data->dx[i]*=(-3.0e0*exp(1.0e0+4.0e0*data->x[i]));
        data->dx[i]+=(1.0e0+2.0e0*data->x[i]);
        data->dx[i]*=(4.0e0*exp(4.0e0*data->y[i]));
        data->dx[i]-=(36.0e0*exp(4.0e0*data->x[i])*(1.0e0+8.0e0*(data->x[i]-1.0e0)*tmp1));
        tmp3=1.0e0+2.0e0*data->y[i];
        tmp3*=tmp3; //(1+2y)^2
        data->dx[i]*=(exp(-4.0e0*data->x[i]*(1.0e0+data->x[i])-tmp3)/15.0e0);
    }
    for ( size_t i=0 ; i<nn ; ++i ) {
        tmp1=data->x[i]*data->x[i]*data->x[i]; //x^3
        tmp2=data->y[i]*data->y[i];
        tmp3=tmp2*data->y[i]; //y^3
        tmp2*=tmp3; //y^5
        data->dy[i]=(-data->x[i]+20.0e0*tmp1-50.0e0*tmp3+80.0e0*tmp2);
        data->dy[i]*=(-12.0e0*exp(1.0e0+4.0e0*data->x[i]));
        data->dy[i]+=1.0e0;
        data->dy[i]*=(8.0e0*data->y[i]*exp(4.0e0*data->y[i]));
        tmp1=1.0e0-2.0e0*data->x[i];
        tmp1*=tmp1; //(1-2x)^2
        tmp2=1.0e0+2.0e0*data->y[i];
        data->dy[i]-=(36.0e0*exp(4.0e0*data->x[i])*tmp1*tmp2);
        data->dy[i]*=(exp(-4.0e0*data->x[i]*(1.0e0+data->x[i])-tmp2*tmp2)/15.0e0);
    }
}
// =============================================================================
// AUX =========================================================================
int MockWavefrontGenerator::Size(void) {
    return data->x.size();
}
void MockWavefrontGenerator::Save(const char *fileName) {
    ofstream ofil(fileName);
    int nn=data->x.size();
    ofil << std::scientific << std::setprecision(12);
    double tmp=data->x[0];
    for ( int i=0 ; i<nn ; ++i ) {
        if ( data->x[i]!=tmp ) {
            ofil << endl;
            tmp=data->x[i];
        }
        ofil << data->x[i] << " " << data->y[i] << " " << data->z[i];
        ofil << " " << data->dx[i] << " " << data->dy[i];
        ofil  << endl;
    }
    ofil.close();
}
void MockWavefrontGenerator::Print(void) {
    int nn=data->x.size();
    cout << "     x             y                z                dx              dy" << endl;
    for ( int i=0 ; i<80 ; ++i ) {
        putchar('*');
    }
    cout << endl;
    cout << std::scientific << std::setprecision(8);
    for ( int i =0 ; i<nn ; ++i ) {
        cout << std::setw(16)
             << data->x[i] << std::setw(16)
             << data->y[i] << std::setw(16)
             << data->z[i] << std::setw(16)
             << data->dx[i] << std::setw(16)
             << data->dy[i] << endl;
    }
}
//  Centers the result along the z axis.
void MockWavefrontGenerator::CenterWavefrontAlongZ(void) {
    int nn=data->z.size();
    if ( nn<=0 ) {
        cerr << "Error: non valid size of array!" << endl;
    }
    double sum=0.0e0;
    for ( int i=0 ; i<nn ; ++i ) {
        sum+=data->z[i];
    }
    sum/=double(nn);
    for ( int i=0 ; i<nn ; ++i ) {
        data->z[i]-=sum;
    }
}
// =============================================================================
#endif  /* _MOCKWAVEFRONTGENERATOR_CPP_ */
