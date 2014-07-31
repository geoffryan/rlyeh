#include "../rlyeh.h"
#include "hydro.h"

int calc_cons(struct Grid *theGrid)
{
    int i,j,k;

    double x[3];
    struct Cell *c;

    for(k=0; k<theGrid->Nx3; k++)
        for(i=0; i<theGrid->Nx1[k]; i++)
            for(j=0; j<theGrid->Nx2[k][i]; j++)
            {
                c = &(theGrid->theCells[k][i][j]);
                x[X1_DIR] = c->xiph[X1_DIR] - c->dx[X1_DIR];
                x[X2_DIR] = c->xiph[X2_DIR] - c->dx[X2_DIR];
                x[X3_DIR] = c->xiph[X3_DIR] - c->dx[X3_DIR];
                prim2cons(theGrid->theCells[k][i][j].prim, theGrid->theCells[k][i][j].cons, x);
            }

}

int calc_prim(struct Grid *theGrid)
{
    int i,j,k;

    double x[3];
    struct Cell *c;

    for(k=0; k<theGrid->Nx3; k++)
        for(i=0; i<theGrid->Nx1[k]; i++)
            for(j=0; j<theGrid->Nx2[k][i]; j++)
            {
                c = &(theGrid->theCells[k][i][j]);
                x[X1_DIR] = c->xiph[X1_DIR] - c->dx[X1_DIR];
                x[X2_DIR] = c->xiph[X2_DIR] - c->dx[X2_DIR];
                x[X3_DIR] = c->xiph[X3_DIR] - c->dx[X3_DIR];
                cons2prim(theGrid->theCells[k][i][j].cons, theGrid->theCells[k][i][j].prim, x);
            }

}
