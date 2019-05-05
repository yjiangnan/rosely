//
//  integral.hpp
//  
//
//  Created by Jiang-Nan Yang on 9/10/18.
//
//

#ifndef integral_hpp
#define integral_hpp

#include <stdio.h>

double fVls2 (double V, double s2, double n, double ElogV, double VlogV0, double Vmax);
double EVfV  (double V, double s2, double n, double ElogV, double VlogV0, double Vmax, double M);
double EsgmfV(double V, double s2, double n, double ElogV, double VlogV0, double Vmax, double M);
double VsgmfV(double V, double s2, double n, double ElogV, double VlogV0, double Vmax, double M, double Esgm);

double fs2lV (double s2,  double n, double V);
double fs2alV(double s2a, double a, double n, double V);

double  fVls2a(double V, double s2, double n, double EVa, double Vga, double a, double n0, double Vmax);

double EVfVa  (double V, double s2, double n, double EVa, double Vga, double a, double n0, double Vmax, double M);
double EsgmfVa(double V, double s2, double n, double EVa, double Vga, double a, double n0, double Vmax, double M);
double VsgmfVa(double V, double s2, double n, double EVa, double Vga, double a, double n0, double Vmax, double M, double Esgm);

double ms2fs2a  (double s2a, double a, double n, double V);
double Vsmplfs2a(double s2a, double a, double n, double V, double ms2);

double z_score(double p);

#endif /* integral_hpp */
