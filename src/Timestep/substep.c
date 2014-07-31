#include <stdlib.h>
#include "../rlyeh.h"
#include "../Boundary/boundary.h"
#include "../Grid/grid.h"
#include "../Hydro/hydro.h"

int substep(struct Grid *theGrid, double dt)
{
    int i,k;
    int nf1 = 0;
    int nf2 = 0;
    int nf3 = 0;
    struct Face *theFaces1 = NULL;
    struct Face *theFaces2 = NULL;
    struct Face *theFaces3 = NULL;

    for(k=0; k<theGrid->Nx3; k++)
        if (theGrid->Nx1[k] > 1)
        {
            setup_faces_x1(theGrid, &theFaces1, &nf1);
            break;
        }
    
    for(k=0; k<theGrid->Nx3; k++)
        for(i=0; i<theGrid->Nx1[k]; i++)
            if (theGrid->Nx2[k][i] > 1)
            {
                setup_faces_x2(theGrid, &theFaces2, &nf2);
                break;
            }
    if (theGrid->Nx3 > 1)
        setup_faces_x3(theGrid, &theFaces3, &nf3);

    reconstruct(theGrid, theFaces1, theFaces2, theFaces3, nf1, nf2, nf3);

    if(theFaces1 != NULL)
    {
        hydro_flux(theGrid, theFaces1, nf1, X1_DIR, dt);
        free(theFaces1);
    }
    if(theFaces2 != NULL)
    {
        hydro_flux(theGrid, theFaces2, nf2, X2_DIR, dt);
        free(theFaces2);
    }
    if(theFaces3 != NULL)
    {
        hydro_flux(theGrid, theFaces3, nf3, X3_DIR, dt);
        free(theFaces3);
    }

    hydro_source(theGrid, dt);

    calc_prim(theGrid);

    synchronize(theGrid);

    apply_boundary(theGrid);
}
