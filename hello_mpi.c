#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(){
    MPI_Init(NULL,NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Process %d out of %d is on %s\n", rank, size, processor_name);
    float* arr = NULL;
    int send_num = 1024*1024*1024/size;
    float* local_arr = (float*)malloc(sizeof(float)*send_num);

    if(rank==0){
        arr = (float*)malloc(sizeof(float)*1024*1024*1024);
        for(int i = 0; i < 1024*1024*1024; i++){
            arr[i] = (float)i;
        }
    }
    // printf("%d\n", send_num);

    MPI_Scatter(arr, send_num, MPI_FLOAT, local_arr, send_num, MPI_FLOAT, 0, MPI_COMM_WORLD);
    printf("%d scatter completed!\n", rank);
    
    if(rank==0){
        free(arr);
    }
    free(local_arr);
    MPI_Finalize();
    return 0;
}