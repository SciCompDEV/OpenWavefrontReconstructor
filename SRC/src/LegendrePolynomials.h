#ifndef _LEGENDRE_POLYNOMIALS_H_
#define _LEGENDRE_POLYNOMIALS_H_

class LegendrePolynomials {
public:
    /** Returns the normalized Associated Legendre polynomial Pnm evaluated at x.  */
    static double NormalPnm(int n,int m,double x);
    /** Returns the normalized Legendre polynomial Pn evaluated at x.  */
    static inline double NormalPn(int n,double x) {return NormalPnm(n,0,x);}
    /** Returns dPn/dt, evaluated at t=x  */
    static double FirstDerivativeNormalPn(int n,double x);
    /** Auxiliar function (used for saving computing time)  */
    static double FirstDerivativeNormalPnFromPnAndPnplus1(int n,\
            double x,double pn,double pnp1);
    /** Auxiliar function (used for saving computing time)  */
    static double FirstDerivativeNormalPnFromPnAndPnminus1(int n,\
            double x,double pn,double pnm1);
};

#endif /* defined(_LEGENDRE_POLYNOMIALS_H_) */
