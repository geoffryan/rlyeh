#include <stdio.h>
#include "par.h"
#include "rlyeh.h"
#include "domain.h"

int main(int argc, char *argv[])
{
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
    free_domain(&dom);

    return 0;
}
