#ifndef RLYEH_GEOM
#define RLYEH_GEOM

double (*geom_dV)(double *, double *);
double (*geom_dA1)(double *, double, double);
double (*geom_dA2)(double *, double, double);
double (*geom_dA3)(double *, double, double);

double geom_dV_cart(double *x, double *dx);
double geom_dA1_cart(double *x, double dx2, double dx3);
double geom_dA2_cart(double *x, double dx1, double dx3);
double geom_dA3_cart(double *x, double dx1, double dx2);

#endif
