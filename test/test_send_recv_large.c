#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>


int MPI_Send_Large(void *buf, uint64_t count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, int element_size) {
    uint64_t int_max = (uint64_t)(INT_MAX);
    uint64_t block_num = count / int_max;
    #ifdef DEBUG
        printf("send_block_num = %lu\n", block_num);
    #endif
    for (uint64_t block_id = 0; block_id < block_num; block_id++){
        MPI_Send(&(((char*)buf)[block_id * int_max * (uint64_t)element_size]), INT_MAX, datatype, dest, tag + (int)block_id, comm);
        // printf("address of %lu * int_max position: %p\n", block_id, &(((char*)buf)[block_id * int_max]));
    }
    if(count- int_max*block_num != 0){
        MPI_Send(&(((char*)buf)[block_num*int_max * (uint64_t)element_size]), (int)(count - int_max*block_num), datatype, dest, tag + (int)block_num, comm);
        #ifdef DEBUG
            printf(" send extend last elements\n");
        #endif
    }
    return 0;
}

int MPI_Recv_Large(void *buf, uint64_t count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status, int element_size) {
    uint64_t int_max = (uint64_t)(INT_MAX);
    uint64_t block_num = count / int_max;
    #ifdef DEBUG
        printf("recv_block_num = %lu\n", block_num);
    #endif
    for (uint64_t block_id = 0; block_id < block_num; block_id++){
        MPI_Recv(&(((char*)buf)[block_id * int_max * (uint64_t)element_size]), INT_MAX, datatype, source, tag + (int)block_id, comm, status);
    }
    if(count- int_max*block_num != 0){
        MPI_Recv(&(((char*)buf)[block_num * int_max * (uint64_t)element_size]), (int)(count- int_max*block_num), datatype, source, tag + (int)block_num, comm, status);
        #ifdef DEBUG
            printf(" recv extend last elements\n");
        #endif 
    }   
    return 0;
}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // uint64_t large = ((uint64_t)INT_MAX) * ((uint64_t)2);

    // uint64_t large = (uint64_t)INT_MAX * (uint64_t)2;
    uint64_t large = (uint64_t)4 * (uint64_t)1024 * (uint64_t)1024 * (uint64_t)1024;
    if(rank == 0){
        printf("large = %lu\n", (unsigned long)large);
    }
    float *send_large_data;
    float *recv_large_data;
    MPI_Status status;

    if(rank==0){
        send_large_data = (float*)malloc(sizeof(float)*large);
        printf("address of start position: %p\n", send_large_data);
        printf("address of int_max position: %p\n", &(send_large_data[INT_MAX]));
        printf("address of int_max position: %p\n", &(send_large_data[(uint64_t)INT_MAX]));
        // send_large_data[large - 1] = 0.314;
        for(uint64_t j = 10; j > 0; j--){
            send_large_data[large - j] = (float)j;
        }
        // send_large_data[large/2] = 0.314;
        // printf("send last float number: %f\n", send_large_data[large - 1]);
        MPI_Send_Large(send_large_data, large, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, sizeof(float));
        free(send_large_data);
    } else{
        recv_large_data = (float*)malloc(sizeof(float)*large);
        MPI_Recv_Large(recv_large_data, large, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status, sizeof(float));
        for(uint64_t j = 10; j > 0; j--){
            printf("recv last %d float number: %f\n", (int)j, recv_large_data[large - j]);
        }
        
        // printf("recv large / 2 float number: %f\n", recv_large_data[large / 2]);
        free(recv_large_data);
    }
    
    MPI_Finalize();
    return 0;
}