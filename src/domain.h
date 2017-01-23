#ifndef RLYEH_DOMAIN
#define RLYEH_DOMAIN

#include "rlyeh.h"

double face_loc(double a, double b, int i, int N, int scale, 
                    struct ParList *pars);
int setup_domain(struct Domain *dom, struct ParList *pars);
int free_domain(struct Domain *dom);

#endif
