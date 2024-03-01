#ifndef D_BASIS_F_H
#define D_BASIS_F_H
#include <math.h>

// GAUSSIANS
inline double G(double r, double eta, double miu) {
    return exp(-eta*(r-miu)*(r-miu));
}
inline double G(double r, double eta, double miu, double f) {
    return exp(-eta*(r-miu)*(r-miu))*f;
}
inline double dG(double r, double eta, double miu) {
    return -2.0*eta*(r-miu)*exp(-eta*(r-miu)*(r-miu));
}
inline double dG(double r, double eta, double miu, double f, double fp) {
    return exp(-eta*(r-miu)*(r-miu) )*(-2.0*f*eta*(r-miu)+fp);
}

// BLIPS
inline double B(double x_c)
    // x_c = eta*(rij-r_s)
    // return B(x_c)
{
    if (-2.0 < x_c && x_c <= -1.0) {
        double b = 2.0 + x_c;
        return 0.25 * b*b*b;
    }
    else if (-1.0 < x_c && x_c <= 1.0) {
        double x2 = x_c*x_c;
        double x3 = x2*x_c;
        return 1.0 - 1.5*x2 + 0.75*std::fabs(x3);
    }
    else if (1.0 < x_c && x_c < 2.0) {
        double b = 2.0 - x_c;
        return 0.25 * b*b*b;
    }
    else {
        return 0.0;
    }
}
inline double B(double x_c, double f)
    // x_c = eta*(rij-r_s)
    // return B(x_c)*f
{
    if (-2.0 < x_c && x_c <= -1.0) {
        double b = 2.0 + x_c;
        return 0.25 * b*b*b *f;
    }
    else if (-1.0 < x_c && x_c <= 1.0) {
        double x2 = x_c*x_c;
        double x3 = x2*x_c;
        return f*(1.0 - 1.5*x2 + 0.75*std::fabs(x3));
    }
    else if (1.0 < x_c && x_c < 2.0) {
        double b = 2.0 - x_c;
        return f*(0.25 * b*b*b);
    }
    else {
        return 0.0;
    }
}
inline double dB(double x_c, double eta)
    // def: x_c = eta*(rij-r_s)
    // return d/dr_ij B(x_c) = eta*d/dx_c B(x_c)
{
    if (-2.0 < x_c && x_c <= -1.0) {
        double d = 2.0+x_c;
        return 0.75*eta*d*d;
    }
    else if (-1.0 < x_c && x_c <= 1.0) {
        return -3.0*eta*x_c + 2.25*eta*x_c*std::fabs(x_c);
    }
    else if (1.0 < x_c && x_c < 2.0) {
        double d = 2.0-x_c;
        return -0.75*eta*d*d;
    }
    else {
        return 0.0;
    }
}
inline double dB(double x_c, double eta, double f, double fp)
    // def: x_c = eta*(rij-r_s)
    // return d/dr_ij f(r_ij)*B(x_c)
{
    if (-2.0 < x_c && x_c <= -1.0) {
        double d = 2.0+x_c;
        return f*(0.75*eta*d*d)+fp*(0.25 * d*d*d);
    }
    else if (-1.0 < x_c && x_c <= 1.0) {
        return f*(-3.0*eta*x_c + 2.25*eta*x_c*std::fabs(x_c))
            +fp*(1.0 - 1.5*x_c*x_c + 0.75*std::fabs(x_c*x_c*x_c));
    }
    else if (1.0 < x_c && x_c < 2.0) {
        double d = 2.0-x_c;
        return f*(-0.75*eta*d*d)+fp*(0.25 * d*d*d);
    }
    else {
        return 0.0;
    }
}
#endif
