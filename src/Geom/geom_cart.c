#include "geom.h"

double geom_dV_cart(double *x, double *dx)
{
    return dx[0]*dx[1]*dx[3];
}

double geom_dA1_cart(double *x, double dx2, double dx3)
{
    return dx2*dx3;
}

double geom_dA2_cart(double *x, double dx1, double dx3)
{
    return dx1*dx3;
}
double geom_dA3_cart(double *x, double dx1, double dx2)
{
    return dx1*dx2;
}
