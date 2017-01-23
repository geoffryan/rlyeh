#include <stdio.h>
#include "par.h"
#include "rlyeh.h"
#include "domain.h"
#include "timestep.h"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    if(argc < 2)
    {
        printf("\nPlease provide a parameter file:\n");
        printf("\n   $ rlyeh <parameter file>\n\n");
        return 0;
    }

    struct ParList pars = PAR_DEFAULT;
    read_pars(&pars, argv[1]);
    print_pars(&pars, NULL);

    struct Domain dom;
    setup_domain(&dom, &pars);

    double dt = calc_dt(&dom);
    printf("t=%.8lg dt=%.8lg\n", dom.t, dt);

    free_domain(&dom);

    MPI_Finalize();

    return 0;
}
