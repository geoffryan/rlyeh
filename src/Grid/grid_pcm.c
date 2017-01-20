#include "../rlyeh.h"
#include "grid.h"

int recon_pcm(struct Grid *theGrid, struct Face *theFaces1, struct Face *theFaces2, struct Face *theFaces3, int nf1, int nf2, int nf3)
{
    int i,j,k,q;
    struct Cell *c;

    for(k=0; k<theGrid->Nx3; k++)
        for(i=0; i<theGrid->Nx1[k]; i++)
            for(j=0; j<theGrid->Nx2[k][i]; j++)
            {
                c = &(theGrid->theCells[k][i][j]);
                for(q=0; q<NUMQ; q++)
                {
                    c->grad1[q] = 0.0;
                    c->grad2[q] = 0.0;
                    c->grad3[q] = 0.0;
                }
            }
    return 0;
}
