#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>
#include <string.h>

// 是否计算每个进程发送，接收数据量
#define CACULATE_TRANSFER_DATA 
// 是否进行排序结果验证
#define CHECK 
// 是否输出语句
#define DEBUG
// 是否计算并输出每段花费时间
#define STEPTIME

// 初始数据形态，0代表随机，1代表完全逆序，2代表完全顺序， 3代表近似顺序
const int initial_data = 0;

void matrix_gen(float *a, uint64_t N, float seed);
uint64_t get_Num_from_cmd(int argc, char** argv, int rank, int is_test);
void check_correct(float* sort_arr, char* check_file, uint64_t N);
void gen_correct_file(char* check_file, uint64_t N, int rank, int force, int is_test, float seed);

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
    // 定义两个变量来存储发送总量和接收总量
    uint64_t send_count = 0;
    uint64_t recv_count = 0;
    // 定义记录步骤时间的变量
    clock_t one_step_start, one_step_end;
    clock_t two_step_start, two_step_end;
    clock_t three_step_start, three_step_end;
    clock_t four_step_start, four_step_end;
    double time_taken;

    uint64_t local_n = n / size;

    #ifdef DEBUG
        if(rank == 0){
            printf("local_n: %lu\n", local_n);

            if(local_n > (uint64_t)INT_MAX){
                printf("Warning: data loss may occur during conversion.\n");
            }
        }
    #endif

    float *local_arr = (float *)malloc(local_n * sizeof(float));
    MPI_Status status;
    
    // MPI_Scatter(arr, send_num, MPI_FLOAT, local_arr, send_num, MPI_FLOAT, 0, MPI_COMM_WORLD);

    #ifdef STEPTIME
        one_step_start = clock();
    #endif
    if(local_n <= (uint64_t)INT_MAX){
        #ifdef DEBUG
            if(rank == 0){
                printf("local_n{%lu} <= INT_MAX\n", local_n);
            }
        #endif
        MPI_Scatter(arr, (int)local_n, MPI_FLOAT, local_arr, (int)local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    } else{
        if(rank==0){
            // rank==0的时候，我们直接将数据复制到local_arr中
            memcpy(local_arr, arr, local_n * sizeof(float));
            for(int i = 1; i < size; i++){
                MPI_Send_Large(&(arr[local_n * i]), local_n, MPI_FLOAT, i, i * 1000, MPI_COMM_WORLD, sizeof(float));
            }
        } else{
            MPI_Recv_Large(local_arr, local_n, MPI_FLOAT, 0, rank * 1000, MPI_COMM_WORLD, &status, sizeof(float));
        }
    }

    #ifdef STEPTIME
        MPI_Barrier(MPI_COMM_WORLD);
        one_step_end = clock();
        if(rank == 0){
            time_taken = (double)(one_step_end - one_step_start) / CLOCKS_PER_SEC;
            printf("first step time: %f\n", time_taken);
        }
    #endif

    #ifdef CACULATE_TRANSFER_DATA
        if(rank == 0){
            send_count += local_n * (size - 1) * sizeof(float);
        } else{
            recv_count += local_n * sizeof(float);
        }
    #endif
    
    #ifdef DEBUG
        if(rank == 0){
            printf("scatter compelete!\n");
        }
    #endif
    //局部排序(快排)
    #ifdef STEPTIME
        two_step_start = clock();
    #endif
    qsort(local_arr, local_n, sizeof(float), compare);
    
    // #ifdef DEBUG
    //     printf("[rank %d]: local quick sort compelete!\n", rank);
    //     printf("rank[%d]: ", rank);
    //     for(int i = 0; i < 32; i++){
    //         printf("%f ", local_arr[i * 1024 * 1024]);
    //     }
    //     printf("\n");
    // #endif
    // 保证所有进程都完成了局部排序
    MPI_Barrier(MPI_COMM_WORLD);
    #ifdef STEPTIME
        two_step_end = clock();
        time_taken = (double)(two_step_end - two_step_start) / CLOCKS_PER_SEC;
        printf("second step time: %f\n", time_taken);
        three_step_start = clock();
    #endif
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

            #ifdef CACULATE_TRANSFER_DATA
                send_count += local_n * sizeof(float);
                recv_count += local_n * sizeof(float);
            #endif
            
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
                // #ifdef DEBUG
                //     printf("rank[%d]: ", rank);
                //     for(int i = 0; i < 32; i++){
                //         printf("%f ", local_arr[i * 1024 * 1024]);
                //     }
                //     printf("\n");
                // #endif
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

                }
            }
            
            free(recv_arr);
            free(temp_arr);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
    }
    #ifdef STEPTIME
        three_step_end = clock();
        time_taken = (double)(three_step_end - three_step_start) / CLOCKS_PER_SEC;
        printf("third step time: %f\n", time_taken);
    #endif

    #ifdef DEBUG
        // printf("[rank %d]: local quick sort compelete!\n", rank);
        printf("rank[%d]: ", rank);
        for(int i = 0; i < 32; i++){
            printf("%f ", local_arr[i * 1024 * 1024]);
        }
        printf("\n");
    #endif

    #ifdef STEPTIME
        four_step_start = clock();
    #endif
    if(rank==0){
        // float* sorted_arr = (float*)malloc(sizeof(float) * n);
        memcpy(arr, local_arr, local_n * sizeof(float));
        for(int j = 1; j < size; j++){
            MPI_Recv_Large(&(arr[local_n * j]), local_n, MPI_FLOAT, j, j * 1000, MPI_COMM_WORLD, &status, sizeof(float));
            #ifdef CACULATE_TRANSFER_DATA
                recv_count += local_n * sizeof(float);
            #endif
        }

        // #ifdef CHECK
        //     qsort(arr, n, sizeof(float), compare);
        //     int correct = memcmp(arr, sorted_arr, sizeof(float) * n);
        //     char* correct_str = (correct == 0 ? "true" : "false");
        //     printf("rank[%d]: sort correct: %s\n", rank, correct_str);
        // #endif

        // free(sorted_arr);  
    } else{
        MPI_Send_Large(local_arr, local_n, MPI_FLOAT, 0, rank * 1000, MPI_COMM_WORLD, sizeof(float));
        #ifdef CACULATE_TRANSFER_DATA
            send_count += local_n * sizeof(float);
        #endif
    }

    #ifdef STEPTIME
        four_step_end = clock();
        time_taken = (double)(four_step_end - four_step_start) / CLOCKS_PER_SEC;
        printf("third step time: %f\n", time_taken);
    #endif

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
        if(initial_data == 1){
            qsort(arr, N, sizeof(float), re_compare);
        } else if(initial_data == 2){
            qsort(arr, N, sizeof(float), compare);
        } else if(initial_data == 3){
            qsort(arr, N, sizeof(float), compare);
            uint64_t exchange_num  = 0;
            uint64_t exchange_total = (uint64_t)N / (uint64_t)100;
            while(exchange_num < exchange_total){
                uint64_t rand_index = ((uint64_t)rand() / (uint64_t)RAND_MAX) * (N - 1);
                if(arr[rand_index] != arr[rand_index + 1]){
                    float temp = arr[rand_index];
                    arr[rand_index] = arr[rand_index + 1];
                    arr[rand_index + 1] = temp;
                    exchange_num += 1;
                }
            }
        }
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
        for(int i = 0; i < 20; i++){
            printf("%f ", arr[1024 * 1024 * i]);
        }
        printf("\n");
        printf("odd and even took %f seconds to sort %s elements.\n", time_taken, argv[1]);
        #ifdef CHECK
            char check_file[20] = "./sorted_";
            strcat(check_file, argv[1]);
            printf("check file: %s\n", check_file);
            gen_correct_file(check_file, N, 0, 0, 0, seed);
            check_correct(arr, check_file, N);
        #endif
        free(arr);
    }

    MPI_Finalize();
    return 0;
}
