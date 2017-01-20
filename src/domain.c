#include <stdlib.h>
#include "rlyeh.h"

int setup_domain(struct Domain *dom, struct ParList *pars)
{
    int N1 = pars->nx1;
    int N2 = pars->nx2;
    int N3 = pars->nx3;
    dom->pars = pars;

    double x1a = pars->x1min;
    double x1b = pars->x1max;
    double x2a = pars->x2min;
    double x2b = pars->x2max;
    double x3a = pars->x3min;
    double x3b = pars->x3max;

    dom->t = pars->tmin;
    dom->cfl = pars->cfl;
    dom->plm = pars->plm;

    dom->N3 = pars->nx3;
    dom->N2 = (int *) malloc(N3 * sizeof(int));
    dom->N1 = (int **) malloc(N3 * sizeof(int *));
    dom->prim = (double ***) malloc(N3 * sizeof(double **));
    dom->cons = (double ***) malloc(N3 * sizeof(double **));
    dom->consRK = (double ***) malloc(N3 * sizeof(double **));
    dom->grad = (double ***) malloc(N3 * sizeof(double **));
    dom->x1f = (double ***) malloc(N3 * sizeof(double **));
    dom->x1fRK = (double ***) malloc(N3 * sizeof(double **));
    dom->w1 = (double ***) malloc(N3 * sizeof(double **));
    dom->x2f = (double **) malloc(N3 * sizeof(double *));
    dom->x2fRK = (double **) malloc(N3 * sizeof(double *));
    dom->w2 = (double **) malloc(N3 * sizeof(double *));
    dom->x3f = (double *) malloc((N3+1) * sizeof(double));
    dom->x3fRK = (double *) malloc((N3+1) * sizeof(double));
    dom->w3 = (double *) malloc((N3+1) * sizeof(double));

    int i,j,k,q;

    for(k=0; k<N3; k++)
    {
        dom->N2[k] = N2;
        dom->N1[k] = (int *) malloc(N2 * sizeof(int));
        dom->prim[k] = (double **) malloc(N2 * sizeof(double *));
        dom->cons[k] = (double **) malloc(N2 * sizeof(double *));
        dom->consRK[k] = (double **) malloc(N2 * sizeof(double *));
        dom->grad[k] = (double **) malloc(N2 * sizeof(double *));
        dom->x1f[k] = (double **) malloc(N2 * sizeof(double *));
        dom->x1fRK[k] = (double **) malloc(N2 * sizeof(double *));
        dom->w1[k] = (double **) malloc(N2 * sizeof(double *));
        dom->x2f[k] = (double *) malloc((N2+1) * sizeof(double));
        dom->x2fRK[k] = (double *) malloc((N2+1) * sizeof(double));
        dom->w2[k] = (double *) malloc((N2+1) * sizeof(double));
        dom->x3f[k] = x3a + k*(x3b-x3a)/N3;
        dom->x3fRK[k] = dom->x3f[k];
        dom->w3[k] = 0.0;
        
        for(j=0; j<dom->N2[k]; j++)
        {
            dom->N1[k][j] = N1;
            dom->prim[k][j] = (double *) malloc(N1*NUMQ * sizeof(double));
            dom->cons[k][j] = (double *) malloc(N1*NUMQ * sizeof(double));
            dom->consRK[k][j] = (double *) malloc(N1*NUMQ * sizeof(double));
            dom->grad[k][j] = (double *) malloc(3*N1*NUMQ * sizeof(double));
            dom->x1f[k][j] = (double *) malloc((N1+1) * sizeof(double));
            dom->x1fRK[k][j] = (double *) malloc((N1+1) * sizeof(double));
            dom->w1[k][j] = (double *) malloc((N1+1) * sizeof(double));
            dom->x2f[k][j] = x2a + j*(x2b-x2a)/N2;
            dom->x2fRK[k][j] = dom->x2f[k][j];
            dom->w2[k][j] = 0.0;

            for(i=0; i<dom->N1[k][j]; i++)
            {
                dom->x1f[k][j][i] = x1a + i*(x1b-x1a)/N1;
                dom->x1fRK[k][j][i] = dom->x1f[k][j][i];
                dom->w1[k][j][i] = 0.0;

                for(q=0; q<NUMQ; q++)
                {
                    dom->prim[k][j][NUMQ*i+q] = 0.0;
                    dom->cons[k][j][NUMQ*i+q] = 0.0;
                    dom->consRK[k][j][NUMQ*i+q] = 0.0;
                    dom->grad[k][j][NUMQ*(3*i+0)+q] = 0.0;
                    dom->grad[k][j][NUMQ*(3*i+1)+q] = 0.0;
                    dom->grad[k][j][NUMQ*(3*i+2)+q] = 0.0;
                }
            }
            dom->x1f[k][j][i] = x1b;
            dom->x1fRK[k][j][i] = dom->x1f[k][j][i];
            dom->w1[k][j][i] = 0.0;
        }
        dom->x2f[k][j] = x2b;
        dom->x2fRK[k][j] = dom->x2f[k][j];
        dom->w2[k][j] = 0.0;
    }
    dom->x3f[k] = x1a;
    dom->x3fRK[k] = dom->x3f[k];
    dom->w3[k] = 0.0;

    return 0;
}

int free_domain(struct Domain *dom)
{
    int j,k;

    for(k=0; k<dom->N3; k++)
    {
        for(j=0; j<dom->N2[k]; j++)
        {
            free(dom->prim[k][j]);
            free(dom->cons[k][j]);
            free(dom->consRK[k][j]);
            free(dom->grad[k][j]);
            free(dom->x1f[k][j]);
            free(dom->x1fRK[k][j]);
            free(dom->w1[k][j]);
        }
        free(dom->prim[k]);
        free(dom->cons[k]);
        free(dom->consRK[k]);
        free(dom->grad[k]);
        free(dom->x1f[k]);
        free(dom->x1fRK[k]);
        free(dom->w1[k]);
        free(dom->x2f[k]);
        free(dom->x2fRK[k]);
        free(dom->w2[k]);
        free(dom->N1[k]);
    }
    free(dom->prim);
    free(dom->cons);
    free(dom->consRK);
    free(dom->grad);
    free(dom->x1f);
    free(dom->x1fRK);
    free(dom->w1);
    free(dom->x2f);
    free(dom->x2fRK);
    free(dom->w2);
    free(dom->x3f);
    free(dom->x3fRK);
    free(dom->w3);
    free(dom->N1);
    free(dom->N2);

    return 0;
}
