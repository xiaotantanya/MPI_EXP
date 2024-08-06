#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>
#include <string.h>

void matrix_gen(float *a, uint64_t N, float seed);
uint64_t get_Num_from_cmd(int argc, char** argv, int rank, int is_test);

int compare (const void * a, const void * b)
{
    return ( *(float*)a - *(float*)b > 0 ? 1 :- 1);
}

int MPI_Send_Large(void *buf, uint64_t count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, int element_size) {
    uint64_t int_max = (uint64_t)(INT_MAX);
    uint64_t block_num = count / int_max;
    printf("send_block_num = %lu\n", block_num);
    for (uint64_t block_id = 0; block_id < block_num; block_id++){
        MPI_Send(&(buf[block_id * int_max * (uint64_t)element_size]), INT_MAX, datatype, dest, tag + (int)block_id, comm);
        printf("address of %lu * int_max position: %p\n", block_id, &(buf[block_id * int_max]));
    }
    if(count- int_max*block_num != 0){
        MPI_Send(&(buf[block_num*int_max * (uint64_t)element_size]), (int)(count - int_max*block_num), datatype, dest, tag + (int)block_num, comm);
        printf(" send extend last elements\n");
    }
    
}

int MPI_Recv_Large(void *buf, uint64_t count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status, int element_size) {
    uint64_t int_max = (uint64_t)(INT_MAX);
    uint64_t block_num = count / int_max;
    printf("recv_block_num = %lu\n", block_num);
    for (uint64_t block_id = 0; block_id < block_num; block_id++){
        MPI_Recv(&(buf[block_id * int_max * (uint64_t)element_size]), INT_MAX, datatype, source, tag + (int)block_id, comm, status);
    }
    if(count- int_max*block_num != 0){
        MPI_Recv(&(buf[block_num * int_max * (uint64_t)element_size]), (int)(count- int_max*block_num), datatype, source, tag + (int)block_num, comm, status);
        printf(" recv extend last elements\n");
    }   
    
}

int MPI_Sendrecv_Large(void *sendbuf, uint64_t sendcount, MPI_Datatype sendtype,
                        int dest, int sendtag, void *recvbuf, uint64_t recvcount,
                        MPI_Datatype recvtype, int source, int recvtag,
                        MPI_Comm comm, MPI_Status *status, int element_size)
    {
        uint64_t int_max = (uint64_t)(INT_MAX);
        // 注意：这里我们默认sendcount等于recvcount
        uint64_t block_num = sendcount / int_max;
        printf("send_recv_block_num = %lu\n", block_num);
        for (uint64_t block_id = 0; block_id < block_num; block_id++){
            MPI_Sendrecv(&(sendbuf[block_id * int_max * (uint64_t)element_size]), INT_MAX, sendtype, 
                            dest, sendtag + (int)block_id, &(recvbuf[block_id * int_max * (uint64_t)element_size]), INT_MAX,
                            MPI_FLOAT, source,  recvtag + (int)block_id,
                            comm, status);
        }

        if(sendcount- int_max*block_num != 0){
        MPI_Sendrecv(&(sendbuf[block_num * int_max * (uint64_t)element_size]), (int)(sendcount- int_max*block_num), sendtype, 
                            dest, sendtag + (int)block_num, &(recvbuf[block_num * int_max * (uint64_t)element_size]), (int)(sendcount- int_max*block_num),
                            MPI_FLOAT, source,  recvtag + (int)block_num,
                            comm, status);
        printf(" sendrecv extend last elements\n");
    }   
    }

void odd_even_sort_parallel(float *arr, uint64_t n, int rank, int size) {
    uint64_t local_n = n / size;
    // printf("n: %ld\n", n);
    // printf("size: %d\n", size);
    printf("local_n: %lu\n", local_n);

    // 检查是否会丢失数据
    if (local_n > UINT_MAX) {
        printf("Warning: data loss may occur during conversion.\n");
    }
    // 强制转换
    // int send_num = (int)local_n;

    float *local_arr = (float *)malloc(local_n * sizeof(float));
    MPI_Status status;
    
    // MPI_Scatter(arr, send_num, MPI_FLOAT, local_arr, send_num, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(local_n <= (uint64_t)INT_MAX){
        MPI_Scatter(arr, (int)local_n, MPI_FLOAT, local_arr, (int)local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    } else{
        if(rank==0){
            // rank==0的时候，我们直接将数据复制到local_arr中
            memcpy(local_arr, arr, local_n * sizeof(float));
            for(int i=1;i<size;i++){
                MPI_Send_Large(&(arr[local_n * i]), local_n, MPI_FLOAT, i, i * 1000, MPI_COMM_WORLD, sizeof(float));
            }
        } else{
            for(int i = 1; i < size; i++){
                MPI_Recv_Large(local_arr, local_n, MPI_FLOAT, 0, i * 1000, MPI_COMM_WORLD, &status, sizeof(float));
            }
        }
    }
    // printf("%d: scatter compelete!\n", rank);
    //局部排序(快排)
    qsort(local_arr, local_n, sizeof(float), compare);
    printf("%d: local quick sort compelete!\n", rank);
    MPI_Barrier(MPI_COMM_WORLD);

        // 奇数阶段和偶数阶段切换
    for (int phase = 0; phase < size; phase++) {
        int partner = (phase + rank) % 2 == 0 ? rank + 1 : rank - 1;
        if (partner >= 0 && partner < size) {
            float *recv_arr = (float *)malloc(local_n * sizeof(float));
            if(local_n <= (uint64_t)INT_MAX){
                MPI_Sendrecv(local_arr, (int)local_n, MPI_FLOAT, partner, 0,
                            recv_arr, (int)local_n, MPI_FLOAT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else{
                MPI_Sendrecv_Large(local_arr, local_n, MPI_FLOAT, partner, 2024,
                            recv_arr, local_n, MPI_FLOAT, partner, 2024, MPI_COMM_WORLD, MPI_STATUS_IGNORE, sizeof(float));
            }
            
            // printf("38\n");
            float *temp_arr = (float*)malloc(2 * local_n * sizeof(float));
            if (rank < partner) {
                for (uint64_t i = 0; i < local_n; i++) temp_arr[i] = local_arr[i];
                for (uint64_t i = 0; i < local_n; i++) temp_arr[i + local_n] = recv_arr[i];
                //这里rank进程只需要前local_n个小的，partner只需要后面local_n个大的
                float* first_point = temp_arr;
                float* second_point = temp_arr + local_n;
                uint64_t  local_index = 0;
                // printf("init local index: %lu\n", local_index);
                while(local_index < local_n){

                    if (*first_point < *second_point){
                        local_arr[local_index] = *first_point;
                        first_point = first_point + 1;
                    } else{
                        local_arr[local_index] = *second_point;
                        second_point = second_point + 1;
                    }
                    local_index = local_index + 1;
                    // printf("rank < partner local_index: %lu\n", local_index);
                }
                //两个指针，偏移
                //加一步local_arr赋值操作
            } else {
                for (uint64_t i = 0; i < local_n; i++) temp_arr[i] = recv_arr[i];
                for (uint64_t i = 0; i < local_n; i++) temp_arr[i + local_n] = local_arr[i];
                //这里rank进程只需要后local_n个大的，partner只需要签名local_n个小的
                float* first_point = temp_arr +  local_n - 1;
                float* second_point = temp_arr + 2 * local_n - 1;
                uint64_t  local_index = local_n - 1;
                uint64_t local_sum = 0;
                // printf("init local index: %lu\n", local_index);
                // printf("69\n");
                while(local_sum < local_n){
                    if (*first_point < *second_point){
                        local_arr[local_index] = *second_point;
                        second_point = second_point - 1;
                    } else{
                        local_arr[local_index] = *first_point;
                        first_point = first_point - 1;
                    }
                    local_index = local_index - 1;
                    local_sum = local_sum + 1;
                    // printf("rank > partner local_index: %llu\n", local_index);
                    // printf("rank > partner local_sum: %llu\n", local_sum);
                    // printf("rank > partner local_n: %llu\n", local_n);
                }
            }
            
            free(recv_arr);
            free(temp_arr);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // int global_sorted;
        // MPI_Allreduce(&sorted_local, &global_sorted, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        // sorted_local = global_sorted;
    }
    if(local_n <= (uint64_t)INT_MAX){
        MPI_Gather(local_arr, (int)local_n, MPI_FLOAT, arr, (int)local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    } else{
        if(rank==0){
            memcpy(arr, local_arr, local_n);
            for(int j = 1; j < size; j++){
                MPI_Recv_Large(&(arr[local_n * j]), local_n, MPI_FLOAT, j, j * 1000, MPI_COMM_WORLD, &status, sizeof(float));
            }
        } else{
            for(int j = 1; j < size; j++){
                MPI_Send_Large(local_arr, local_n, MPI_FLOAT, 0, j * 1000, MPI_COMM_WORLD, sizeof(float));
            }
        }
    }
    
    free(local_arr);
}

void test(int rank, int size){
    float test_arr[10] = {4.4,2.3,4.6,5.3,7.5,3.2,1.4,6.7,9.1,1.3};
    uint64_t n = 10;
    odd_even_sort_parallel(test_arr, n, rank, size);
    for(int i = 0; i < 10; i++){
        printf("%f ",test_arr[i]);
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    int is_test = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(rank==0){
        if (argc < 2){
            printf("Usage: %s <number_str>\n", argv[0]);
            return 1;
        }
    }
    uint64_t N = get_Num_from_cmd(argc, argv, rank, is_test);
    MPI_Barrier(MPI_COMM_WORLD);
    // uint64_t N = 16;
    float *arr = NULL;
    float seed = 0.344;
    clock_t start_time, end_time;
    if (rank == 0) {
        arr = (float *)malloc(N * sizeof(float));
        matrix_gen(arr, N, seed);
        // printf("Unsorted array:\n");
        // for (int i = 0; i < N; i++) printf("%f ", arr[i]);
        // printf("\n");
        start_time = clock();
    }

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Process %d out of %d is on %s\n", rank, size, processor_name);
    
    odd_even_sort_parallel(arr, N, rank, size);
    // test(rank, size);
    if (rank == 0) {
        // printf("Sorted array:\n");
        // for (int i = 0; i < N; i++) printf("%f ", arr[i]);
        // printf("\n");
        end_time = clock();
        // 计算排序时间
        double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        // 打印排序时间
        printf("qsort took %f seconds to sort %s elements.\n", time_taken, argv[1]);
        if(is_test) {
            for (uint64_t i = 0; i < N; i++){
                printf("%f ", arr[i]);
            }
            printf("\n");
        }
        free(arr);
    }

    MPI_Finalize();
    return 0;
}
