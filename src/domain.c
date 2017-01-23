#include <stdlib.h>
#include "rlyeh.h"
#include "comm.h"
#include "domain.h"

double face_loc(double a, double b, int i, int N, int scale, 
                    struct ParList *pars)
{
    if(scale == 1)
    {
        if(i < N)            
            return a * pow(b/a, ((float) i)/N);
        else
            return b * pow(b/a, ((float) i-N)/N);
    }
    else
    {
        if(i < N)            
            return a + i*(b-a)/N;
        else
            return b + (i-N)*(b-a)/N;
    }
}

int setup_domain(struct Domain *dom, struct ParList *pars)
{
    int N1 = pars->nx1;
    int N2 = pars->nx2;
    int N3 = pars->nx3;
    dom->pars = pars;

    double X1a = pars->x1min;
    double X1b = pars->x1max;
    double X2a = pars->x2min;
    double X2b = pars->x2max;
    double X3a = pars->x3min;
    double X3b = pars->x3max;
    int Ng = pars->nghost;

    dom->t = pars->tmin;
    dom->cfl = pars->cfl;
    dom->plm = pars->plm;

    dom->Ng3 = Ng;
    dom->Ng1[0] = Ng;
    dom->Ng1[1] = Ng;

    dom->N3glob = N3;
    dom->N2glob = N2;
    dom->N1glob = N1;

    dom->periodic1 = (pars->bc1a == 1 || pars->bc1b == 1) ? 1 : 0;
    dom->periodic2 = (pars->bc2a == 1 || pars->bc2b == 1) ? 1 : 0;
    dom->periodic3 = (pars->bc3a == 1 || pars->bc3b == 1) ? 1 : 0;

    setup_comm(dom);

    int n1, n2, n3, ia, ib, ja, jb, ka, kb;
    int *dim = dom->mpidim;
    int *coord = dom->mpicoord;
    
    n1 = N1;
    ia = 0;
    ib = n1;
    
    n2 = N2 / dim[1];
    ja = coord[1] * (N2/dim[1]);
    if(coord[1] < N2 % dim[1])
    {
        n2++;
        ja += coord[1];
    }
    else
        ja += N2 % dim[1];
    jb = ja + n2;

    n3 = N3 / dim[0];
    ka = coord[0] * (N3/dim[0]);
    if(coord[0] < N3 % dim[0])
    {
        n3++;
        ka += coord[0];
    }
    else
        ka += N3 % dim[0];
    kb = ka + n3;

    n1 += 2*Ng;
    n2 += 2*Ng;
    n3 += 2*Ng;
    ia -= Ng;
    ja -= Ng;
    ka -= Ng;
    ib += Ng;
    jb += Ng;
    kb += Ng;

    dom->N3 = n3;
    dom->N2 = (int *) malloc(n3 * sizeof(int));
    dom->N1 = (int **) malloc(n3 * sizeof(int *));
    dom->Ng2 = (int *) malloc(2*n3 * sizeof(int));
    dom->prim = (double ***) malloc(n3 * sizeof(double **));
    dom->cons = (double ***) malloc(n3 * sizeof(double **));
    dom->consRK = (double ***) malloc(n3 * sizeof(double **));
    dom->grad = (double ***) malloc(n3 * sizeof(double **));
    dom->x1f = (double ***) malloc(n3 * sizeof(double **));
    dom->x1fRK = (double ***) malloc(n3 * sizeof(double **));
    dom->w1 = (double ***) malloc(n3 * sizeof(double **));
    dom->x2f = (double **) malloc(n3 * sizeof(double *));
    dom->x2fRK = (double **) malloc(n3 * sizeof(double *));
    dom->w2 = (double **) malloc(n3 * sizeof(double *));
    dom->x3f = (double *) malloc((n3+1) * sizeof(double));
    dom->x3fRK = (double *) malloc((n3+1) * sizeof(double));
    dom->w3 = (double *) malloc((n3+1) * sizeof(double));

    int i,j,k,q;

    for(k=0; k<dom->N3; k++)
    {
        dom->Ng2[2*k] = Ng;
        dom->Ng2[2*k+1] = Ng;
        dom->N2[k] = n2;
        dom->N1[k] = (int *) malloc(n2 * sizeof(int));
        dom->prim[k] = (double **) malloc(n2 * sizeof(double *));
        dom->cons[k] = (double **) malloc(n2 * sizeof(double *));
        dom->consRK[k] = (double **) malloc(n2 * sizeof(double *));
        dom->grad[k] = (double **) malloc(n2 * sizeof(double *));
        dom->x1f[k] = (double **) malloc(n2 * sizeof(double *));
        dom->x1fRK[k] = (double **) malloc(n2 * sizeof(double *));
        dom->w1[k] = (double **) malloc(n2 * sizeof(double *));
        dom->x2f[k] = (double *) malloc((n2+1) * sizeof(double));
        dom->x2fRK[k] = (double *) malloc((n2+1) * sizeof(double));
        dom->w2[k] = (double *) malloc((n2+1) * sizeof(double));

        dom->x3f[k] = face_loc(X3a, X3b, ka+k, N3, pars->x3type, pars);
        dom->x3fRK[k] = dom->x3f[k];
        dom->w3[k] = 0.0;
        
        for(j=0; j<dom->N2[k]; j++)
        {
            dom->N1[k][j] = n1;
            dom->prim[k][j] = (double *) malloc(n1*NUMQ * sizeof(double));
            dom->cons[k][j] = (double *) malloc(n1*NUMQ * sizeof(double));
            dom->consRK[k][j] = (double *) malloc(n1*NUMQ * sizeof(double));
            dom->grad[k][j] = (double *) malloc(3*n1*NUMQ * sizeof(double));
            dom->x1f[k][j] = (double *) malloc((n1+1) * sizeof(double));
            dom->x1fRK[k][j] = (double *) malloc((n1+1) * sizeof(double));
            dom->w1[k][j] = (double *) malloc((n1+1) * sizeof(double));
            dom->x2f[k][j] = face_loc(X2a, X2b, ja+j, N2, pars->x2type, pars);
            dom->x2fRK[k][j] = dom->x2f[k][j];
            dom->w2[k][j] = 0.0;

            for(i=0; i<dom->N1[k][j]; i++)
            {
                int N1int = dom->N1[k][j]-dom->Ng1[0]-dom->Ng1[1];
                dom->x1f[k][j][i] = face_loc(X1a, X1b, ia+i, N1int,
                                                pars->x1type, pars);
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
            dom->x1f[k][j][i] = face_loc(X1a, X1b, ia+i, dom->N1[k][j], 
                                                pars->x1type, pars);
            dom->x1fRK[k][j][i] = dom->x1f[k][j][i];
            dom->w1[k][j][i] = 0.0;
        }
        dom->x2f[k][j] = face_loc(X2a, X2b, ja+j, N2, pars->x2type, pars);
        dom->x2fRK[k][j] = dom->x2f[k][j];
        dom->w2[k][j] = 0.0;
    }
    dom->x3f[k] = face_loc(X3a, X3b, ka+k, N3, pars->x3type, pars);
    dom->x3fRK[k] = dom->x3f[k];
    dom->w3[k] = 0.0;

    for(k=0; k<=dom->N3; k++)
        printf("%.3g ", dom->x3f[k]);
    printf("\n\n");
    for(k=0; k<dom->N3; k++)
    {
        for(j=0; j<=dom->N2[k]; j++)
            printf("%.3g ", dom->x2f[k][j]);
        printf("\n");
    }
    printf("\n\n");
    for(k=0; k<dom->N3; k++)
    {
        for(j=0; j<dom->N2[k]; j++)
        {
            for(i=0; i<=dom->N1[k][j]; i++)
                printf("%.3g ", dom->x1f[k][j][i]);
            printf("\n");
        }
        printf("\n");
    }

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
