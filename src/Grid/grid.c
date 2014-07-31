#include <stdlib.h>
#include "../rlyeh.h"
#include "grid.h"
#include "../Geom/geom.h"

int face_initialize(struct Face *myface, struct Cell *cL, struct Cell *cR, double dx1, double dx2, int dir)
{
    myface->cL = cL;
    myface->cR = cR;
    myface->cm[X1_DIR] = cL->xiph[X1_DIR]-0.5*cL->dx[X1_DIR];
    myface->cm[X2_DIR] = cL->xiph[X2_DIR]-0.5*cL->dx[X2_DIR];
    myface->cm[X3_DIR] = cL->xiph[X3_DIR]-0.5*cL->dx[X3_DIR];
    myface->cm[dir] = cL->xiph[dir];

    if(dir == X1_DIR)
        myface->dA = geom_dA1(myface->cm, dx2, cL->dx[X3_DIR]);
    else if(dir == X2_DIR)
        myface->dA = geom_dA2(myface->cm, dx1, cL->dx[X3_DIR]);
    else
        myface->dA = geom_dA3(myface->cm, dx1, dx2);

    return 0;
}

// Name lovingly stolen from Paul's code.
void build_jloop(struct Grid *theGrid, struct Face *theFaces, int *nf, int i1, int i2, int k)
{
    int j1, j2;
    int nf1;
    int flag, neighbour;
    int **nx2 = theGrid->Nx2;
    struct Cell ***theCells = theGrid->theCells;

    double xLp, xLm, xRp, xRm, dx;

    j2 = 0;
    nf1 = 0;
    for(j1=0; j1<nx2[k][i1]; j1++)
    {
        flag = 0;
        xLp = theCells[k][i1][j1].xiph[X1_DIR];
        xLm = theCells[k][i1][j1].xiph[X1_DIR]-theCells[k][i1][j1].dx[X1_DIR];
        xRp = theCells[k][i2][j2].xiph[X1_DIR];
        xRm = theCells[k][i2][j2].xiph[X1_DIR]-theCells[k][i2][j2].dx[X1_DIR];
        neighbour = (xRp > xLm && xRm < xLp);
        while(!flag || neighbour)
        {
            if(neighbour)
            {
                if(theFaces != NULL)
                {
                    flag = 1;
                    if(xLp > xRp && xLm > xRm)
                        dx = xRp-xLm;
                    else if(xLp > xRp && xLm < xRm)
                        dx = xRp-xRm;
                    else if(xLp < xRp && xLm > xRm)
                        dx = xLp-xLm;
                    else
                        dx = xLp-xRm;
                    face_initialize(&(theFaces[*nf+nf1]), &(theCells[k][i1][j1]), &(theCells[k][i2][j2]), 0.0, dx, X1_DIR);
                }
                nf1++;
            }
            j2++;
            if(j2 == nx2[k][i2])
                j2 = 0;
            xRp = theCells[k][i2][j2].xiph[X1_DIR];
            xRm = theCells[k][i2][j2].xiph[X1_DIR]-theCells[k][i2][j2].dx[X1_DIR];
            neighbour = (xRp > xLm && xRm < xLp);
        }
        j2--;
        if(j2==-1)
            j2 = nx2[k][i2]-1;
    }

    *nf += nf1;
}

int setup_faces_x1(struct Grid *theGrid, struct Face **theFaces, int *nf)
{
    //Allocate and initialize all constant-x1 faces.
    
    int i,k;

    int nf1;
    int *nx1 = theGrid->Nx1;
    int nx3 = theGrid->Nx3;

    //Count the faces.
    
    nf1 = 0;
    for(k=0; k<nx3; k++)
        for(i=0; i<nx1[k]-1; i++)
            build_jloop(theGrid, NULL, &nf1, i, i+1, k); 

    //Allocate faces.
    *nf = nf1;
    *theFaces = (struct Face *) malloc(nf1 * sizeof(struct Face));

    //Initialize the faces.
    nf1 = 0;
    for(k=0; k<nx3; k++)
        for(i=0; i<nx1[k]-1; i++)
            build_jloop(theGrid, *theFaces, &nf1, i, i+1, k); 

    return 0;
}

int setup_faces_x2(struct Grid *theGrid, struct Face **theFaces, int *nf)
{
    //Allocate and initialize all constant-x2 faces.

    int i,j,k;
    int nf2 = 0;
    int fnum = 0;
    int *nx1 = theGrid->Nx1;
    int **nx2 = theGrid->Nx2;
    int nx3 = theGrid->Nx3;
    struct Face *myface;

    //Count the faces.
    if(theGrid->periodic & 2) //We're periodic in x2
        for(k=0; k<nx3; k++)
            for(i=0; i<nx1[k]; i++)
                nf2 += nx2[k][i];
    else
        for(k=0; k<nx3; k++)
            for(i=0; i<nx1[k]; i++)
                nf2 += nx2[k][i] - 1;

    //Allocate the faces.
    *nf = nf2;
    *theFaces = (struct Face *) malloc(nf2 * sizeof(struct Face));

    //Initialize the faces.
    fnum = 0;
    for(k=0; k<nx3; k++)
        for(i=0; i<nx1[k]; i++)
        {
            for(j=0; j<nx2[k][i]; j++)
            {
                face_initialize(&(*theFaces[fnum]), 
                        &(theGrid->theCells[k][i][j]), 
                        &(theGrid->theCells[k][i][j+1]), 
                        theGrid->theCells[k][i][j].dx[X1_DIR], 0.0, X2_DIR);
                fnum++;
            }

            if(theGrid->periodic & 2) //We're periodic in x2
            {
                face_initialize(&(*theFaces[fnum]), 
                        &(theGrid->theCells[k][i][j]), 
                        &(theGrid->theCells[k][i][0]), 
                        theGrid->theCells[k][i][j].dx[X1_DIR], 0.0, X2_DIR);
                fnum++;
            }
        }

    return 0;
}

int setup_faces_x3(struct Grid *theGrid, struct Face **theFaces, int *nf)
{
    *nf = 0;
    return 0;
}
