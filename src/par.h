#ifndef RLYEH_PAR
#define RLYEH_PAR

enum{VAR_DBL, VAR_INT, VAR_LON, VAR_STR};

struct ParList
{
    int hydro;
    int geom;
    int recon;
    int riemann;
    int step;
    int bc1a;
    int bc1b;
    int bc2a;
    int bc2b;
    int bc3a;
    int bc3b;
    int nx1;
    int nx2;
    int nx3;
    int nghost;
    int nc;
    int np;
    double tmin;
    double tmax;
    double x1min;
    double x1max;
    int x1type;
    int x1pifac;
    double x2min;
    double x2max;
    int x2type;
    int x2pifac;
    double x3min;
    double x3max;
    int x3type;
    int x3pifac;
    double plm;
    double cfl;
    int eos;
    int cool;
    double gammalaw;

    int io;
    int nChkpt;

    double geomPar1;
    double geomPar2;
    double geomPar3;
    double geomPar4;

    int init;
    double initPar1;
    double initPar2;
    double initPar3;
    double initPar4;
    double initPar5;
    double initPar6;
    double initPar7;
    double initPar8;
};

const static struct ParList PAR_DEFAULT = {
    .hydro = 0,
    .geom = 0,
    .recon = 0,
    .riemann = 0,
    .step = 0,
    .bc1a = 0,
    .bc1b = 0,
    .bc2a = 0,
    .bc2b = 0,
    .bc3a = 0,
    .bc3b = 0,
    .nx1 = 1,
    .nx2 = 1,
    .nx3 = 1,
    .nghost = 1,
    .nc = 3,
    .np = 0,
    .tmin = 0.0,
    .tmax = 1.0,
    .x1min = 0.0,
    .x1max = 1.0,
    .x1type = 0,
    .x1pifac = 0,
    .x2min = 0.0,
    .x2max = 1.0,
    .x2type = 0,
    .x2pifac = 0,
    .x3min = 0.0,
    .x3max = 1.0,
    .x3type = 0,
    .x3pifac = 0,
    .plm = 1.5,
    .cfl = 0.5,
    .eos = 0,
    .cool = 0,
    .gammalaw = 1.4,

    .io = 0,
    .nChkpt = 0,

    .geomPar1 = 0.0,
    .geomPar2 = 0.0,
    .geomPar3 = 0.0,
    .geomPar4 = 0.0,

    .init = 0,
    .initPar1 = 0.0,
    .initPar2 = 0.0,
    .initPar3 = 0.0,
    .initPar4 = 0.0,
    .initPar5 = 0.0,
    .initPar6 = 0.0,
    .initPar7 = 0.0,
    .initPar8 = 0.0
};

int readvar(char filename[], char key[], int vtype, void *ptr);
void read_pars(struct ParList *theParList, char filename[]);
void print_pars(struct ParList *theParList, char filename[]);

#endif
