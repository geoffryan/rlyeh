#ifndef RLYEH_HYDRO
#define RHYEH_HYDRO

int (*prim2cons)(double *, double *, double *);
int (*cons2prim)(double *, double *, double *);
int (*hydro_flux)(struct Grid *, struct Face *, int, int, double);
int (*hydro_source)(struct Grid *, double);
double (*max_speed)(double *, double *);

int calc_prim(struct Grid *theGrid);
int calc_cons(struct Grid *theGrid);

int prim2cons_newt(double *prim, double *cons, double *x);
int cons2prim_newt(double *cons, double *prim, double *x);
int hydro_flux_newt(struct Grid *theGrid, struct Face *theFaces, int nf, int dir, double dt);
int hydro_source_newt(struct Grid *theGrid, double dt);
double max_speed_newt(double *prim, double *x);

#endif
