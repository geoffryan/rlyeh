#ifndef RLYEH
#define RLYEH

#define NUMC 4
#define NUMP 0
#define NUMQ (NUMC+NUMP)

#include "par.h"

enum{X1,X2,X3};

struct Domain{
    double ***prim;
    double ***cons;
    double ***consRK;
    double ***grad;
    double ***x1f;
    double ***x1fRK;
    double **x2f;
    double **x2fRK;
    double *x3f;
    double *x3fRK;
    double ***w1;
    double **w2;
    double *w3;
    int **N1;
    int *N2;
    int N3;

    int periodic1;
    int periodic2;
    int periodic3;

    double t;
    double cfl;
    double plm;

    struct ParList *pars;
};


#endif
