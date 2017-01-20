#ifndef RLYEH_HYDRO
#define RLYEH_HYDRO

int prim2cons(double *prim, double *cons, double *x, double dV,
                double a, double *b, double *g, double t);
int cons2prim(double *prim, double *cons, double *x, double dV,
                double a, double *b, double *g, double t);
int flux(double *prim, double *f, int dim, double *x, 
            double a, double *b, double *g, double t);
int source(double *prim, double *s, double *x1, double *x2, double *x3,
            double a, double *b, double *g, double t);
int wavespeeds(double *prim, double *vel, int dim, double *x, 
                double a, double *b, double *g, double t);
int maxv(double *prim, double *vel, double *x, 
            double a, double *b, double *g, double t);

int cons2prim(double *cons, double *prim, double *x, double dV,
                double a, double *b, double *g, double t);

#endif
