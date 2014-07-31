#ifndef RLYEH
#define RLYEH

#define NUMC 4
#define NUMP 0
#define NUMQ (NUMC+NUMP)

enum{X1_DIR,X2_DIR,X3_DIR};


struct Cell{
    double prim[NUMQ];
    double cons[NUMQ];
    double RKcons[NUMQ];

    double grad1[NUMQ];
    double grad2[NUMQ];
    double grad3[NUMQ];

    double xiph[3];
    double dx[3];
    double dS;
};

struct Face{
    struct Cell *cL;
    struct Cell *cR;
    
    double cm[3];
    double dA;
};

struct Grid{
    struct Cell ***theCells;
    int ng;
    int *Nx1;
    int **Nx2;
    int Nx3;
    int periodic;

    double *x1ph;
    double *x3ph;
};

int (*synchronize)(struct Grid *);

#endif
