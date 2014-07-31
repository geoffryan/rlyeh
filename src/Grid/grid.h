#ifndef RLYEH_GRID
#define RLYEH_GRID

int (*reconstruct)(struct Grid *, struct Face *, struct Face *, struct Face *, int, int, int);

int face_initialize(struct Face *myface, struct Cell *cL, struct Cell *cR, double dx1, double dx2, int dir);
int setup_faces_x1(struct Grid *theGrid, struct Face **theFaces, int *nf);
int setup_faces_x2(struct Grid *theGrid, struct Face **theFaces, int *nf);
int setup_faces_x3(struct Grid *theGrid, struct Face **theFaces, int *nf);

int recon_plm(struct Grid *, struct Face *, struct Face *, struct Face *, int, int, int);

#endif
