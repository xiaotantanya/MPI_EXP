#include <stdio.h>
#include <mpi.h>
#include <limits.h>
#include <stdlib.h>


int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int *a;
    int *b;
    int num = INT_MAX;
    MPI_Status status;
    // for(int i = 0; i < size; i++){
    //     if(rank != i){
    //         MPI_Send(&a, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    //     }
    // }

    // for(int i = 0; i < size; i++){
    //     if(rank != i){
    //         MPI_Recv(&b, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
    //     }
    // }
    a = (int*)malloc(sizeof(int)*num);
    b = (int*)malloc(sizeof(int)*num);
    if(rank==0){
        MPI_Send(&a, num, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&b, num, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
    } else{
        
        MPI_Send(&a, num, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(&b, num, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        printf("b[0] = %d\n", b[0]);
    }
    free(a);
    free(b);
    MPI_Finalize();
}