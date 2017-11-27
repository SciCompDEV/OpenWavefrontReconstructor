#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include "LegendrePolynomials.h"
#include <cmath>

#ifndef DEBUG
#define DEBUG 0
#endif


double LegendrePolynomials::NormalPnm(int n,int m,double x)
{
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



double LegendrePolynomials::FirstDerivativeNormalPn(int n,double x)
{
   double np1=double(n+1);
   if ( x==1.0e0 ) {
      return 0.5e0*double(n)*np1/sqrt(double(2*n+1));
   }
   if ( x==-1.0e0 ) {
      if ( n&1 ) {
         return 0.5e0*double(n)*np1/sqrt(double(2*n+1));
      } else {
         return -0.5e0*double(n)*np1/sqrt(double(2*n+1));
      }
   }
   double pn=NormalPnm(n,0,x);
   double pnp1=NormalPnm((n+1),0,x);
   double coeff=sqrt(double(2*n+1)*(np1*np1)/double(2*n+3));
   return ((np1*x*pn-coeff*pnp1)/(1.0e0-x*x));
}



double LegendrePolynomials::FirstDerivativeNormalPnFromPnAndPnplus1(int n,\
      double x,double pn,double pnp1)
{
   double np1=double(n+1);
   double coeff=sqrt(double(2*n+1)*(np1*np1)/double(2*n+3));
   return ((np1*x*pn-coeff*pnp1)/(1.0e0-x*x));
}



double LegendrePolynomials::FirstDerivativeNormalPnFromPnAndPnminus1(int n,\
      double x,double pn,double pnm1)
{
   if ( n==0 ) {return 0.0e0;}
   double nr=double(n);
   double coeff=sqrt((2.0e0*nr+1.0e0)*(nr*nr)/double(2*n-1));
   return ((-nr*x*pn+coeff*pnm1)/(1.0e0-x*x));
}



