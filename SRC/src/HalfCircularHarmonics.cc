#ifndef _HALFCIRCULARHARMONICS_CPP_
#define _HALFCIRCULARHARMONICS_CPP_

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "HalfCircularHarmonics.h"
#include <cmath>
#ifndef DEBUG
#define DEBUG 0
#endif


#define CIRCULARHARMONICS_EPS_TO_ZERO 1.0e-4

double HalfCircularHarmonics::NormalPnm(int n,int m,double x) {
#if DEBUG
    if (n<0) {
        throw(string("Requesting n<0!"));
    }
    if (m<0) {
        throw(string("Requesting m<0!"));
    }
    if (m>n) {
        throw(string("Requesting m>n!"));
    }
    if (fabs(x)>1.0e0) {
        throw(string("Requesting P(x), |x|>1."));
    }
#endif
    double omx2,fct,Pmm;
    Pmm=1.0e0;
    if (m>0) {
        omx2=(-x*x);
        omx2+=1.0e0;
        fct=1.0e0;
        for (int i=1;i<=m;i++) {
            Pmm*=(omx2*fct/(fct+1.0e0));
            fct+=2.0e0;
        }
    }
    Pmm=sqrt(double(2*m+1)*Pmm*0.5e0);
    if (m&1) {Pmm=-Pmm;}
    double Pmmp1,prevfct,k2,Pk=0.0e0;
    if (n==m){
        return Pmm;
    } else {
        Pmmp1=x*sqrt(2.0e0*double(m)+3.0e0)*Pmm;
        if (n==(m+1)){
            return Pmmp1;
        } else {
            prevfct=sqrt(2.0e0*double(m)+3.0e0);
            for (int k=(m+2);k<=n;k++) {
                k2=double(k*k);
                fct=sqrt((4.0e0*k2-1.0e0)/(k2-double(m*m)));
                Pk=(x*Pmmp1-Pmm/prevfct)*fct;
                prevfct=fct;
                Pmm=Pmmp1;
                Pmmp1=Pk;
            }
            return Pk;
        }
    }
}
double HalfCircularHarmonics::FirstDerivativeNormalPnm(int n,int m,double x) {
    double pnm=NormalPnm(n,m,x);
    double pnp1m=NormalPnm((n+1),m,x);
    double np1=double(n+1);
    double coeff=sqrt(double(2*n+1)*(np1*np1-double(m*m))/double(2*n+3));
    return ((np1*x*pnm-coeff*pnp1m)/(1.0e0-x*x));
}

double HalfCircularHarmonics::FirstDerivativeNormalPnmFromPnmAndPnplus1m(int n,int m,\
        double x,double pnm,double pnp1m) {
    double np1=double(n+1);
    double coeff=sqrt(double(2*n+1)*(np1*np1-double(m*m))/double(2*n+3));
    return ((np1*x*pnm-coeff*pnp1m)/(1.0e0-x*x));
}

double HalfCircularHarmonics::FirstDerivativeNormalPnmFromPnmAndPnminus1m(int n,int m,\
        double x,double pnm,double pnm1m) {
    if ( n==0 ) {return 0.0e0;}
    double nr=double(n);
    double coeff=sqrt((2.0e0*nr+1.0e0)*(nr*nr-double(m*m))/double(2*n-1));
    return ((-nr*x*pnm+coeff*pnm1m)/(1.0e0-x*x));
}
int HalfCircularHarmonics::ComputeNumberOfPolynomialsFromPolynomialOrder(int n) {
    return (n*(n+2)+1);
}
double HalfCircularHarmonics::CircularHarmonic(const int n,const int m,const double mu,const double phi) {
    int am=(m<0?(-m):m);
    double pnm=NormalPnm(n,am,mu);
    double ang=((m<0)?(sin(double(am)*phi)):cos(double(am)*phi));
    double normFact=((m==0)?cOneOverSqrtTwoPi:cOneOverSqrtPi);
    return (normFact*pnm*ang);
}
double HalfCircularHarmonics::HalfCircularHarmonicFromCartesian(int n,int m,double x,double y) {
    double theta,mu,phi;
    ComputeThetaMuAndPhiFromXiEta(x,y,theta,mu,phi);
    return CircularHarmonic(n,m,mu,phi);
}
void HalfCircularHarmonics::ComputeThetaMuAndPhiFromXiEta(const double xi,const double eta,
        double &theta,double &mu,double &vphi) {
    double r=sqrt(xi*xi+eta*eta);
    theta=cHalfPi*r;
    mu=cos(theta);
    vphi=((r<CIRCULARHARMONICS_EPS_TO_ZERO)?0.0e0:atan2(eta,xi));
}
double HalfCircularHarmonics::DCircularHarmonicDxiFromCartesian(const int n,const int m,const double xi,const double eta) {
    double theta,mu,phi;
    ComputeThetaMuAndPhiFromXiEta(xi,eta,theta,mu,phi);
    return DCircularHarmonicDxi(n,m,theta,mu,phi);
}
double HalfCircularHarmonics::DCircularHarmonicDxi(const int n,const int m,
        const double theta,const double mu,const double phi) {
    int am=((m<0)?(-m):m);
    double pnm=NormalPnm(n,am,mu);
    double dpnm=FirstDerivativeNormalPnm(n,am,mu);
    double sqrt1mmu2=sqrt(1.0e0-mu*mu);
    double ang,dang;
    double normfact; //stemming from the sin(m phi) or cos(m phi) part.
    if ( m==0 ) {
        ang=1.0e0;
        dang=0.0e0;
        normfact=cOneOverSqrtTwoPi;
    } else {
        ang=(m<0)?(sin(double(am)*phi)):cos(double(am)*phi);
        dang=am*((m<0)?(cos(double(am)*phi)):(-sin(double(am)*phi)));
        normfact=cOneOverSqrtPi;
    }
    return normfact*(-cHalfPi*sqrt1mmu2*cos(phi)*ang*dpnm-(cHalfPi*sin(phi)*dang*pnm/theta));
}
double HalfCircularHarmonics::DCircularHarmonicDetaFromCartesian(const int n,const int m,const double xi,const double eta) {
    double theta,mu,phi;
    ComputeThetaMuAndPhiFromXiEta(xi,eta,theta,mu,phi);
    return DCircularHarmonicDeta(n,m,theta,mu,phi);
}
double HalfCircularHarmonics::DCircularHarmonicDeta(const int n,const int m,const double theta,
        const double mu,const double phi) {
    int am=((m<0)?(-m):m);
    double pnm=NormalPnm(n,am,mu);
    double dpnm=FirstDerivativeNormalPnm(n,am,mu);
    double sqrt1mmu2=sqrt(1.0e0-mu*mu);
    double ang,dang;
    double normfact; //stemming from the sin(m phi) or cos(m phi) part.
    if ( m==0 ) {
        ang=1.0e0;
        dang=0.0e0;
        normfact=cOneOverSqrtTwoPi;
    } else {
        ang=(m<0)?(sin(double(am)*phi)):cos(double(am)*phi);
        dang=am*((m<0)?(cos(double(am)*phi)):(-sin(double(am)*phi)));
        normfact=cOneOverSqrtPi;
    }
    return normfact*(-cHalfPi*sqrt1mmu2*sin(phi)*ang*dpnm+(cHalfPi*cos(phi)*dang*pnm/theta));
}


#endif  /* _HALFCIRCULARHARMONICS_CPP_ */

