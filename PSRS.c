#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>
#include <string.h>

// 定义一个发送字节数，接收字节数变量来完成。
#define CACULATE_TRANSFER_DATA 1
// 是否进行排序结果验证
#define CHECK 1

void matrix_gen(float *a, uint64_t N, float seed);
uint64_t get_Num_from_cmd(int argc, char **argv, int rank, int is_test);

int compare(const void *a, const void *b){
    return (*(float*)a - *(float*)b > 0 ? 1 : -1);
}

int check_sort(float* arr, float* local_sorted_arr, uint64_t start_index, uint64_t total_num, uint64_t local_num){
    qsort(arr, total_num, sizeof(float), compare);
    int is_correct = 0;
    for(uint64_t i = 0; i < local_num; i++){
        if(arr[i + start_index] != local_sorted_arr[i]){
            is_correct = 1;
            break;
        }
    }
    return is_correct;
}

int MPI_Send_Large(void *buf, uint64_t count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, int element_size) {
    uint64_t int_max = (uint64_t)(INT_MAX);
    uint64_t block_num = count / int_max;
    printf("send_block_num = %lu\n", block_num);
    for (uint64_t block_id = 0; block_id < block_num; block_id++){
        MPI_Send(&(buf[block_id * int_max * (uint64_t)element_size]), INT_MAX, datatype, dest, tag + (int)block_id, comm);
        printf("address of %lu * int_max position: %p\n", block_id, &(buf[block_id * int_max]));
    }
    if(count- int_max * block_num != 0){
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
    if(count- int_max * block_num != 0){
        MPI_Recv(&(buf[block_num * int_max * (uint64_t)element_size]), (int)(count- int_max * block_num), datatype, source, tag + (int)block_num, comm, status);
        printf(" recv extend last elements\n");
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

void PSRS(float *arr, uint64_t n, int rank, int size){
    uint64_t local_n = n / size;

    printf("local_n: %lu\n", local_n);

    if(local_n > (uint64_t)INT_MAX){
        printf("Warning: data loss may occur during conversion.\n");
    }

    float *local_arr = (float*)malloc(local_n * sizeof(float));
    MPI_Status status;

    if(local_n <= (uint64_t)INT_MAX){
        if(rank == 0){
            printf("local_n{%lu} <= INT_MAX\n", local_n);
        }
        MPI_Scatter(arr, (int)local_n, MPI_FLOAT, local_arr, (int)local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    } else{
        if(rank == 0){
            memcpy(local_arr, arr, local_n * sizeof(float));
            for(int i = 1; i < size; i++){
                MPI_Send_Large(&(arr[local_n * i]), local_n, MPI_FLOAT, i, i*1000, MPI_COMM_WORLD, sizeof(float));
            }
        } else{
            MPI_Recv_Large(local_arr, local_n, MPI_FLOAT, 0, rank * 1000, MPI_COMM_WORLD, &status, sizeof(float));
        }
    }

    if(rank == 0){
        printf("scatter compelete!\n");
    }
    qsort(local_arr, local_n, sizeof(float), compare);
    printf("[rank %d]: local quick sort compelete!\n", rank);
    MPI_Barrier(MPI_COMM_WORLD);

    // 选出局部数组中的关键元素，共 size-1 个元素。
    float *local_crucial_arr;
    if(local_n % (uint64_t)size != 0){
        printf("local_n / size must be 0 !!!\n");
        free(local_arr);
        exit -1;
    } else if(local_n <= (uint64_t)(size - 1)){
        printf("local_n must be larger than size - 1!!!\n");
        free(local_arr);
        exit -1;
    } else{
        local_crucial_arr = (float*)malloc(sizeof(float) * (size - 1));
        uint64_t w = local_n / size;
        for(int i = 1; i < size; i++){
            local_crucial_arr[i - 1] = local_arr[(uint64_t)i * w];
        }
    }

    if(rank == 0){
        for(int i = 0; i < (size - 1); i++){
            printf("%f ", local_crucial_arr[i]);
        }
        printf("\n");
        printf("local_arr[%lu] = %f\n", local_n / size, local_arr[local_n / size]);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // 将局部数组中的关键元素全发给0进程，然后选出所有进程的关键元素。
    float *total_local_crucial_arr;
    if(rank != 0){
        MPI_Send(local_crucial_arr, size - 1, MPI_FLOAT, 0, rank, MPI_COMM_WORLD);
    } else{
        total_local_crucial_arr = (float*)malloc(sizeof(float) * (size * (size - 1)));
        memcpy(total_local_crucial_arr, local_crucial_arr, (size - 1) * sizeof(float));
        for(int i = 0; i < (size - 1); i++){
            printf("%f ", total_local_crucial_arr[i]);
        }
        printf("\n");
        for(int i = 1; i < size; i++){
            MPI_Recv(&(total_local_crucial_arr[i * (size - 1)]), size - 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, &status);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        for(int i = 0; i < (size - 1) * size; i++){
            printf("%f ", total_local_crucial_arr[i]);
        }
        printf("\n");
    }
    
    // 这里我们需要从（size - 1) * size中选出 size - 1 个元素
    float * root_crucial_arr;
    root_crucial_arr = (float*)malloc(sizeof(float)*(size - 1));
    if(rank == 0){
        float *total_local_crucial_arr_tmp = (float*)malloc(sizeof(float) * (size * (size - 1)));
        uint64_t *now_index = (uint64_t*)malloc(sizeof(uint64_t) * (size));
        for(int i = 0; i < size; i++){
            now_index[i] = i * (size - 1);
        }
        int min_index = 0;
        int sort_num = 0;
        while(sort_num < size * (size - 1)){
            min_index = 0;
            for(int i = 1; i < size; i++){
                if(now_index[min_index] != local_n + 1 && now_index[i] != local_n + 1){
                    if(total_local_crucial_arr[now_index[i]] < total_local_crucial_arr[now_index[min_index]]){
                        min_index = i;
                    }
                } else if(now_index[i] != local_n + 1){
                    min_index = i;
                }
            }
            total_local_crucial_arr_tmp[sort_num] = total_local_crucial_arr[now_index[min_index]];
            now_index[min_index] = now_index[min_index] + 1;
            if(now_index[min_index] >= (uint64_t)(size - 1) * (min_index + 1)){
                now_index[min_index] = local_n + 1;
            }
            sort_num = sort_num + 1;
        }
        for(int i = 1; i < size; i++){
            root_crucial_arr[i - 1] = total_local_crucial_arr_tmp[i * (size - 1)];
        }
        for(int i = 0; i < size * (size - 1); i++){
            printf("%f ", total_local_crucial_arr_tmp[i]);
        }
        printf("\n");
        free(total_local_crucial_arr_tmp);
        free(now_index);
    }
    MPI_Bcast(root_crucial_arr, size - 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("root_crucial_arr be sure, size[%d]\n", size - 1);
    MPI_Barrier(MPI_COMM_WORLD);
    // if(rank == 0){
        printf("rank[%d]: ", rank);
        for(int i = 0; i < size - 1; i++){
            printf("%f ", root_crucial_arr[i]);
        }
        printf("\n");
    // }
    MPI_Barrier(MPI_COMM_WORLD);
    //将本地数据根据主元进行分段，这里保存小段最后一个元素索引来进行分段，第i段元素为[last_ele_index[i], last_ele_index[i+1])
    // 一共有size段
    uint64_t *last_ele_index = (uint64_t*)malloc(sizeof(uint64_t) * (size + 1));
    last_ele_index[0] = 0;
    for(uint64_t i = 1; i < size; i++){
        for(uint64_t j = last_ele_index[i - 1]; j < local_n; j++){
            if(local_arr[j] > root_crucial_arr[i - 1]){
                last_ele_index[i] = j;
                break;
            }
        }
    }
    last_ele_index[size] = local_n;
    if(rank == 0){
        printf("%f\n", local_arr[25165830]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("rank[%d]: ", rank);
    for(int i = 0; i < size + 1; i++){
        printf("%lu ", last_ele_index[i]);
    }
    printf("\n");
    MPI_Barrier(MPI_COMM_WORLD);
    // 将每个进程i第j段的数量发给j进程，并接收其他进程第i段的数量。
    uint64_t *ele_number = (uint64_t*)malloc(sizeof(uint64_t) * (size));
    for(int i = 0; i < size; i++){
        ele_number[i] = last_ele_index[i+1] - last_ele_index[i];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // if(rank == 0){
    printf("rank[%d]: ", rank);
    for(int i = 0; i < size; i++){
        printf("%lu ", ele_number[i]);
    }
    printf("\n");
    // }
    MPI_Barrier(MPI_COMM_WORLD);
    // 设置一个变量接收分配给进程的数据大小
    uint64_t *proc_ele_num = (uint64_t*)malloc(sizeof(uint64_t)*(size));
    float * proc_data;
    for(int i = 0; i < size; i++){
        if(i != rank){
            MPI_Sendrecv(&(ele_number[i]), 1, MPI_UINT64_T, i, 0,
                &(proc_ele_num[i]), 1, MPI_UINT64_T, i, 0, MPI_COMM_WORLD, &status);
        } else {
            proc_ele_num[i] = ele_number[i];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        for(int i = 0; i < size; i++){
            printf("%lu ", proc_ele_num[i]);
        }
        printf("\n");
    }
    uint64_t sum = 0;
    for(int i = 0; i < size; i++){
        sum += proc_ele_num[i];
    }
    // if(rank == 0){
    MPI_Barrier(MPI_COMM_WORLD);
    printf("rank[%d]: ", rank);
    printf("sum = %lu\n", sum);
    // }
    if(CHECK && rank != 0){
        MPI_Send(&sum, 1, MPI_UINT64_T, 0, 345, MPI_COMM_WORLD);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    uint64_t *proc_start_index = malloc(sizeof(uint64_t) * (size));
    proc_start_index[0] = 0;
    for(int i = 1; i < size; i++){
        proc_start_index[i] = proc_start_index[i-1] + proc_ele_num[i-1];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("rank[%d]: ", rank);
    // if(rank == 0){
        for(int i = 0; i < size; i++){
            printf("%lu ", proc_start_index[i]);
        }
        printf("\n");
    // }
    proc_data = (float*)malloc(sizeof(float) * sum);
    for(int i = 0; i < size; i++){
        if(ele_number[i] > (uint64_t)INT_MAX || proc_ele_num[i] > (uint64_t)INT_MAX){
            printf("ele_number[%d] or proc_ele_num[%d] > INT_MAX\n", i, i);
            free(local_arr);
            free(local_crucial_arr);
            free(root_crucial_arr);
            free(last_ele_index);
            free(ele_number);
            free(proc_ele_num);
            free(proc_start_index);
            free(proc_data);
            exit -1;
        }
    }
    for(int i = 0; i < size; i++){
        if(i != rank){
            MPI_Sendrecv(&(local_arr[last_ele_index[i]]), ele_number[i], MPI_FLOAT, i, 0,
                    &(proc_data[proc_start_index[i]]), proc_ele_num[i], MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status);
        } else{
            memcpy(&(proc_data[proc_start_index[i]]), &(local_arr[last_ele_index[i]]), ele_number[i] * sizeof(float));
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //每个进程进行局部归并排序
    float *sort_proc_data = (float*)malloc(sizeof(float) * sum);
    uint64_t *now_index = (uint64_t*)malloc(sizeof(uint64_t) * (size));
    for(int i = 0; i < size; i++){
        now_index[i] = proc_start_index[i];
    }
    int min_index = 0;
    float min_value = 0.;
    int sort_num = 0;
    for(int i = 0; i < size; i++){
        if(proc_ele_num[i] == 0){
            now_index[i] = -1;
        }
    }
    while(sort_num < sum){
        min_index = 0;
        for(int i = 1; i < size; i++){
            if(now_index[min_index] != sum + 1 && now_index[i] != sum + 1){
                if(proc_data[now_index[i]] < proc_data[now_index[min_index]]){
                    min_index = i;
                }
            } else if(now_index[i] != sum + 1){
                min_index = i;
            }
        }
        sort_proc_data[sort_num] = proc_data[now_index[min_index]];
        now_index[min_index] = now_index[min_index] + 1;
        if(now_index[min_index] >= proc_start_index[min_index] + proc_ele_num[min_index]){
            now_index[min_index] = sum + 1;
        }
        sort_num = sort_num + 1;
    }

    printf("rank[%d]: ", rank);
    for(int i = 0; i < 10; i++){
        printf("%f ", sort_proc_data[i]);
    }
    printf("\n");
    // 测试结果是否准确
    if(CHECK){
        uint64_t *local_sum = (uint64_t*)malloc(sizeof(uint64_t) * size);
        uint64_t *start_index = (uint64_t*)malloc(sizeof(uint64_t) * size);
        float *sorted_data = (float*)malloc(sizeof(float) * n);
        if(rank == 0){
            for(int i = 1; i < size; i++){
                MPI_Recv(&(local_sum[i]), 1, MPI_UINT64_T, i, 345, MPI_COMM_WORLD, &status);
            }
            local_sum[0] = sum;
            start_index[0] = 0;
            for(int i = 1; i < size; i++){
                start_index[i] = start_index[i - 1] + local_sum[i - 1]; 
            }
            printf("rank[%d]: ", rank);
            for(int i = 0; i < size; i++){
                printf("%lu ", local_sum[i]);
            }
            printf("\n");
            for(int i = 1; i < size; i++){
                MPI_Recv(&(sorted_data[start_index[i]]), local_sum[i], MPI_FLOAT, i, 254, MPI_COMM_WORLD, &status);
            }
            memcpy(sorted_data, sort_proc_data, sum * sizeof(float));
            qsort(arr, n, sizeof(float), compare);
            int correct = memcmp(sorted_data, arr, sizeof(float) * n);
            char* correct_str = (correct == 0 ? "true" : "false");
            printf("rank[%d]: sort correct: %s\n", rank, correct_str);
        } else{
            MPI_Send(sort_proc_data, sum, MPI_FLOAT, 0, 254, MPI_COMM_WORLD);
        }
        
        // MPI_Bcast(start_index, size, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        printf("rank[%d] start_index: ", rank);
        for(int i = 0; i < size; i++){
            printf("%lu ", start_index[i]);
        }
        printf("\n");
        MPI_Barrier(MPI_COMM_WORLD);
        

        
        free(local_sum);
        free(start_index);
        free(sorted_data);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){
        free(total_local_crucial_arr);
    }
    free(local_arr);
    free(local_crucial_arr);
    free(root_crucial_arr);
    free(last_ele_index);
    free(ele_number);
    free(proc_ele_num);
    free(proc_start_index);
    free(proc_data);
    free(sort_proc_data);
    free(now_index);

}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int rank, size;
    int is_test = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank==0){
        if(argc < 2){
            printf("Usage: %s <number_str>\n", argv[0]);
            exit -1;
        }
    }

    uint64_t N = get_Num_from_cmd(argc, argv, rank, is_test);
    MPI_Barrier(MPI_COMM_WORLD);

    float *arr = NULL;
    float seed = 0.344;
    clock_t start_time, end_time;
    if(rank == 0){
        arr = (float*)malloc(N * sizeof(float));
        matrix_gen(arr, N, seed);
        start_time = clock();
    }

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Process %d out of %d is on %s\n", rank, size, processor_name);

    // 进行排序
    PSRS(arr, N, rank, size);

    if(rank == 0){
        end_time = clock();
        // 计算排序时间
        double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        // 打印排序时间
        printf("PSRS sort took %f seconds to sort %s elements.\n", time_taken, argv[1]);
        if(is_test){
            for(uint64_t i = 0; i < N; i++){
                printf("%f ", arr[i]);
            }
            printf("\n");
        }
        free(arr);
    }
    
    MPI_Finalize();
    return 0;
}