#include <math.h>
#include "rlyeh.h"
#include "timestep.h"
#include "hydro.h"

int forward_euler(struct Domain *dom, double dt)
{
    substep(dom, dt);

    return 0;
}

double calc_dt(struct Domain *dom)
{
    int i,j,k;

    double dtmin = 1.0e100;

    for(k=0; k<dom->N3; k++)
    {
        double w3 = 0.5*(dom->w3[k] + dom->w3[k+1]);
        double x3a = dom->x3f[k];
        double x3b = dom->x3f[k+1];
        double dx3 = x3b-x3a;

        for(j=0; j<dom->N2[k]; j++)
        {
            double w2 = 0.5*(dom->w2[k][j] + dom->w2[k][j+1]);
            double x2a = dom->x2f[k][j];
            double x2b = dom->x2f[k][j+1];
            double dx2 = x2b-x2a;

            for(i=0; i<dom->N1[k][j]; i++)
            {
                double w1= 0.5*(dom->w1[k][j][i] + dom->w1[k][j][i+1]);
                double x1a = dom->x1f[k][j][i];
                double x1b = dom->x1f[k][j][i+1];
                double dx1 = x1b-x1a;

                double x[3] = {0.5*(x1a+x1b),0.5*(x2a+x2b),0.5*(x3a+x3b)};
                double a = 1.0;
                double b[3] = {0.0, 0.0, 0.0};
                double g[9] = {1.0, 0.0, 0.0,
                               0.0, 1.0, 0.0,
                               0.0, 0.0, 1.0};
                double v[6];
                maxv(&(dom->prim[k][j][NUMQ*i]), v, x, a, b, g, dom->t);
                double dv[3][2];
                dv[0][0] = fabs(v[0] - w1);
                dv[0][1] = fabs(v[1] - w1);
                dv[1][0] = fabs(v[2] - w2);
                dv[1][1] = fabs(v[3] - w2);
                dv[2][0] = fabs(v[4] - w3);
                dv[2][1] = fabs(v[5] - w3);

                double idt1, idt2, idt3;

                idt1 = dv[0][0] > dv[0][1] ? dv[0][0]/dx1 : dv[0][1]/dx1;
                idt2 = dv[1][0] > dv[1][1] ? dv[1][0]/dx2 : dv[1][1]/dx2;
                idt3 = dv[2][0] > dv[2][1] ? dv[2][0]/dx3 : dv[2][1]/dx3;
                
                double dt = 1.0 / (idt1 + idt2 + idt3);
                if (dt < dtmin)
                    dtmin = dt;
            }
        }
    }

    return dtmin;
}

