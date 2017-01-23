#include <mpi.h>
#include "rlyeh.h"

int setup_comm(struct Domain *dom)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int ndim = 2;
    int dim[] = {0,0};
    if(dom->N3glob == 1)
        dim[0] = 1;
    MPI_Dims_create(size, ndim, dim);

    int reorder = 1;
    int periods[] = {dom->periodic3, dom->periodic2};
    MPI_Comm *comm = &(dom->comm);
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dim, periods, reorder, comm);

    MPI_Comm_size(*comm, &size);
    MPI_Comm_rank(*comm, &rank);

    int coord[2];
    int neighbour[4];
    MPI_Cart_coords(*comm, rank, ndim, coord);

    int coord3a[2], coord3b[2], coord2a[2], coord2b[2];
    coord3a[0] = coord[0]-1;
    coord3a[1] = coord[1];
    coord3b[0] = coord[0]+1;
    coord3b[1] = coord[1];
    coord2a[0] = coord[0];
    coord2a[1] = coord[1]-1;
    coord2b[0] = coord[0];
    coord2b[1] = coord[1]+1;

    if(periods[0] || coord[0] > 0)
        MPI_Cart_rank(*comm, coord3a, &(neighbour[0])); 
    else
        neighbour[0] = -1;
    if(periods[0] || coord[0] < dim[0]-1)
        MPI_Cart_rank(*comm, coord3b, &(neighbour[1])); 
    else
        neighbour[1] = -1;
    if(periods[1] || coord[1] > 0)
        MPI_Cart_rank(*comm, coord2a, &(neighbour[2])); 
    else
        neighbour[2] = -1;
    if(periods[1] || coord[1] < dim[1]-1)
        MPI_Cart_rank(*comm, coord2b, &(neighbour[3])); 
    else
        neighbour[3] = -1;

    dom->mpirank = rank;
    dom->mpisize = size;
    dom->mpidim[0] = dim[0];
    dom->mpidim[1] = dim[1];
    dom->mpicoord[0] = coord[0];
    dom->mpicoord[1] = coord[1];
    dom->mpineighbour[0] = neighbour[0];
    dom->mpineighbour[1] = neighbour[1];
    dom->mpineighbour[2] = neighbour[2];
    dom->mpineighbour[3] = neighbour[3];

    return 0;
}

int sync(struct Domain *dom)
{
    

    return 0;
}
