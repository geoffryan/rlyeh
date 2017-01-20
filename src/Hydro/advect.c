#include "../rlyeh.h"
#include "../hydro.h"

double c[] = {1.0, 0.0, 0.0};

int prim2cons(double *prim, double *cons, double *x, double dV, double a, double *b, double *g, double t)
{
    int q;
    for(q=0; q<NUMQ; q++)
        cons[q] = prim[q]*dV;

    return 0;
}

int flux(double *prim, double *f, int dim, double *x, double a, double *b, double *g, double t)
{
    int q;
    for(q=0; q<NUMQ; q++)
        f[q] = c[dim] * prim[q];

    return 0;
}

int source(double *prim, double *s, double *x1, double *x2, double *x3, double a, double *b, double *g, double t)
{
    int q;
    for(q=0; q<NUMQ; q++)
        s[q] = 0.0;

    return 0;
}

int wavespeeds(double *prim, double *vel, int dim, double *x, double a, double *b, double *g, double t)
{
    vel[0] = c[dim];
    vel[1] = c[dim];

    return 0;
}

int maxv(double *prim, double *vel, double *x, double a, double *b, double *g, double t)
{
    vel[2*0+0] = c[0];
    vel[2*0+1] = c[0];
    vel[2*1+0] = c[1];
    vel[2*1+1] = c[1];
    vel[2*2+0] = c[2];
    vel[2*2+1] = c[2];

    return 0;
}

int cons2prim(double *prim, double *cons, double *x, double dV,
                double a, double *b, double *g, double t)
{
    int q;
    for(q=0; q<NUMQ; q++)
        prim[q] = cons[q]/dV;

    return 0;
}
