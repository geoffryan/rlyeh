#ifndef RLYEH_TIMESTEP
#define RLYEH_TIMESTEP

#include "rlyeh.h"

int substep(struct Domain *dom, double dt);

int forward_euler(struct Domain *dom, double dt);
double calc_dt(struct Domain *dom);


#endif
