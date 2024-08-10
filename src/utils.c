#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include "utils.h"

float rand_float(float s){
    return 4*s*(1-s);
}

void matrix_gen(float *a, uint64_t N, float seed){
    float s = seed;
    for(uint64_t i=0; i<N; i++){
        s = rand_float(s);
        a[i] = s;
    }
}

uint64_t get_Num_from_cmd(int argc, char** argv, int rank, int is_test){
    char * Num_str = argv[1];
    const uint64_t M_num = 1024 * 1024;
    uint64_t N = 0;
    if(strcmp(Num_str, "256M")==0){
        N = 256 * M_num;
    } else if(strcmp(Num_str, "1G")==0){
        N = 1024 * M_num;
    } else if(strcmp(Num_str, "4G")==0){
        N = 4 * 1024 * M_num;
    } else if(is_test){
        N = atoi(Num_str);
    } else{
        printf("Number Error !\n");
        return 1;
    }
    if (rank ==0) {
        printf("N = %lu\n", N);
    }
    return N;
}

// 顺序
int compare (const void * a, const void * b)
{
    return ( *(float*)a - *(float*)b > 0 ? 1 :- 1);
}

// 逆序
int re_compare (const void * a, const void * b)
{
    return ( *(float*)a - *(float*)b > 0 ? -1 : 1);
}

void gen_correct_file(char* check_file, uint64_t N, int rank, int force, int is_test, float seed){
    if(rank == 0){
        if(force || access(check_file, 0)){
            float* arr_data = (float*)malloc(sizeof(float) * N);
            matrix_gen(arr_data, N, seed);
            qsort(arr_data, N, sizeof(float), compare);

            FILE *pd = NULL;
            pd = fopen(check_file, "wb");
            fwrite(arr_data, sizeof(float), N, pd);
            fclose(pd);
            free(arr_data);
        }
    }
}

void check_correct(float* sort_arr, char* check_file, uint64_t N){
    FILE *pd = NULL;
    pd = fopen(check_file, "rb");
    float *read_arr = (float*)malloc(sizeof(float) * N);
    fread(read_arr, sizeof(float), N, pd);

    int correct = memcmp(sort_arr, read_arr, sizeof(float) * N);
    char* correct_str = (correct == 0 ? "true" : "false");
    printf("rank[%d]: sort correct: %s\n", 0, correct_str);
    free(read_arr);
}

int MPI_Send_Large(void *buf, uint64_t count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, int element_size) {
    uint64_t int_max = (uint64_t)(INT_MAX);
    uint64_t block_num = count / int_max;
    #ifdef DEBUG
        printf("send_block_num = %lu\n", block_num);
    #endif
    for (uint64_t block_id = 0; block_id < block_num; block_id++){
        MPI_Send(&(buf[block_id * int_max * (uint64_t)element_size]), INT_MAX, datatype, dest, tag + (int)block_id, comm);
    }
    if(count- int_max * block_num != 0){
        MPI_Send(&(buf[block_num*int_max * (uint64_t)element_size]), (int)(count - int_max*block_num), datatype, dest, tag + (int)block_num, comm);
        #ifdef DEBUG
            printf(" send extend last elements\n");
        #endif
    }
    
}

int MPI_Recv_Large(void *buf, uint64_t count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status, int element_size) {
    uint64_t int_max = (uint64_t)(INT_MAX);
    uint64_t block_num = count / int_max;
    #ifdef DEBUG
        printf("recv_block_num = %lu\n", block_num);
    #endif
    for (uint64_t block_id = 0; block_id < block_num; block_id++){
        MPI_Recv(&(buf[block_id * int_max * (uint64_t)element_size]), INT_MAX, datatype, source, tag + (int)block_id, comm, status);
    }
    if(count- int_max * block_num != 0){
        MPI_Recv(&(buf[block_num * int_max * (uint64_t)element_size]), (int)(count- int_max * block_num), datatype, source, tag + (int)block_num, comm, status);
        #ifdef DEBUG
            printf(" recv extend last elements\n");
        #endif
    }   
    
}

// 目前该函数只能用于两个进程互换信息，不要多个进程交叉在一起。
// 目前该函数还有bug ！！！
int MPI_Sendrecv_Large(void *sendbuf, uint64_t sendcount, MPI_Datatype sendtype,
                        int dest, int sendtag, void *recvbuf, uint64_t recvcount,
                        MPI_Datatype recvtype, int source, int recvtag,
                        MPI_Comm comm, MPI_Status *status, int element_size)
{
        uint64_t int_max = (uint64_t)(INT_MAX);
        // 注意：这里我们默认sendcount等于recvcount
        uint64_t send_block_num = sendcount / int_max;
        uint64_t recv_block_num = recvcount / int_max;
        uint64_t min_block_num = send_block_num < recv_block_num ? send_block_num : recv_block_num;
        int min_index = send_block_num < recv_block_num ? 0 : 1;
        printf("send_block_num = %lu\n", send_block_num);
        printf("recv_block_num = %lu\n", recv_block_num);
        for (uint64_t block_id = 0; block_id < min_block_num; block_id++){
            MPI_Sendrecv(&(sendbuf[block_id * int_max * (uint64_t)element_size]), INT_MAX, sendtype, 
                            dest, sendtag + (int)block_id, &(recvbuf[block_id * int_max * (uint64_t)element_size]), INT_MAX,
                            MPI_FLOAT, source,  recvtag + (int)block_id,
                            comm, status);
        }
        if(send_block_num - recv_block_num != 0){
            if(min_index == 0){
                for(uint64_t block_id = 0; block_id < recv_block_num - min_block_num; block_id++){
                    MPI_Sendrecv(&(sendbuf), 0, sendtype, dest, sendtag + (int)block_id + (int)min_block_num, &(recvbuf[(block_id + min_block_num) * int_max * (uint64_t)element_size]), INT_MAX,
                            recvtype, source,  recvtag + (int)block_id + (int)min_block_num, comm, status);
                }
            } else{
                for(uint64_t block_id = 0; block_id < send_block_num - min_block_num; block_id++){
                    MPI_Sendrecv(&(sendbuf[(block_id + min_block_num) * int_max * (uint64_t)element_size]), int_max, sendtype, dest, sendtag + (int)block_id + (int)min_block_num, &(recvbuf[(block_id + min_block_num) * int_max * (uint64_t)element_size]), 0,
                            recvtype, source,  recvtag + (int)block_id + (int)(min_block_num), comm, status);
            }
        }
        MPI_Sendrecv(&(sendbuf[send_block_num * int_max * (uint64_t)element_size]), (int)(sendcount- int_max * send_block_num), sendtype, 
                            dest, sendtag + (int)(send_block_num) + (int)(recv_block_num), &(recvbuf[recv_block_num * int_max * (uint64_t)element_size]), (int)(recvcount- int_max * recv_block_num),
                            recvtype, source,  recvtag + (int)(send_block_num) + (int)(recv_block_num),
                            comm, status);
        // printf(" sendrecv extend last elements\n");
    }   
}