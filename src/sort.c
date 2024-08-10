#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <time.h>

#include "utils.h"
#include "sort.h"
#include "setting.h"


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
            printf("rank[0]: step 1 time: %f s.\n", time_taken);
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
    
    // 保证所有进程都完成了局部排序
    MPI_Barrier(MPI_COMM_WORLD);
    #ifdef STEPTIME
        two_step_end = clock();
        time_taken = (double)(two_step_end - two_step_start) / CLOCKS_PER_SEC;
        printf("rank[%d]: step 2 time: %f s.\n", rank, time_taken);
        if(rank == 0){
            three_step_start = clock();
        }
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
    MPI_Barrier(MPI_COMM_WORLD);
    #ifdef STEPTIME
        if(rank == 0){
            three_step_end = clock();
            time_taken = (double)(three_step_end - three_step_start) / CLOCKS_PER_SEC;
            printf("rank[%d]: step 3 time: %f s.\n", rank, time_taken);
        }
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
        if(rank == 0){
            four_step_start = clock();
        }
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

    } else{
        MPI_Send_Large(local_arr, local_n, MPI_FLOAT, 0, rank * 1000, MPI_COMM_WORLD, sizeof(float));
        #ifdef CACULATE_TRANSFER_DATA
            send_count += local_n * sizeof(float);
        #endif
    }

    #ifdef STEPTIME
        if(rank == 0){
            four_step_end = clock();
            time_taken = (double)(four_step_end - four_step_start) / CLOCKS_PER_SEC;
            printf("rank[%d]: step 4 time: %f s.\n", rank, time_taken);
        } 
    #endif

    #ifdef CACULATE_TRANSFER_DATA
        printf("rank[%d]: ", rank);
        printf("send message count: %lu Byte.\n", send_count);
        printf("rank[%d]: ", rank);
        printf("recv message count: %lu Byte.\n", recv_count);
        MPI_Barrier(MPI_COMM_WORLD);
    #endif

    free(local_arr);
}

void PSRS(float *arr, uint64_t n, int rank, int size){
    // 定义两个变量来存储发送总量和接收总量
    uint64_t send_count = 0;
    uint64_t recv_count = 0;
    // 记录每个步骤排序时间
    clock_t step_start[7];
    clock_t step_end[7];

    uint64_t local_n = n / size;
    #ifdef DEBUG
        if(rank == 0){
            printf("local_n: %lu\n", local_n);

            if(local_n > (uint64_t)INT_MAX){
                printf("Warning: data loss may occur during conversion.\n");
            }
        }
    #endif

    float *local_arr = (float*)malloc(local_n * sizeof(float));
    MPI_Status status;

    #ifdef STEPTIME
        if(rank == 0){
            step_start[0] = clock();
        }
    #endif
    if(local_n <= (uint64_t)INT_MAX){
        #ifdef DEBUG
            if(rank == 0){
                printf("local_n{%lu} <= INT_MAX\n", local_n);
            }
        #endif
        MPI_Scatter(arr, (int)local_n, MPI_FLOAT, local_arr, (int)local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    } else{
        if(rank == 0){
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
        if(rank == 0){
            step_end[0] = clock();
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
    #ifdef STEPTIME
        step_start[1] = clock();
    #endif
    qsort(local_arr, local_n, sizeof(float), compare);
    #ifdef DEBUG
        printf("[rank %d]: local quick sort compelete!\n", rank);
    #endif
    // 保证所有进程都完成了局部排序
    MPI_Barrier(MPI_COMM_WORLD);
    #ifdef STEPTIME
        step_end[1] = clock();
    #endif
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
        #ifdef STEPTIME
            step_start[2] = clock();
        #endif
        local_crucial_arr = (float*)malloc(sizeof(float) * (size - 1));
        uint64_t w = local_n / size;
        for(int i = 1; i < size; i++){
            local_crucial_arr[i - 1] = local_arr[(uint64_t)i * w];
        }
    }
    #ifdef DEBUG
        if(rank == 0){
            for(int i = 0; i < (size - 1); i++){
                printf("%f ", local_crucial_arr[i]);
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    #endif
    // 将局部数组中的关键元素全发给0进程，然后选出所有进程的关键元素。
    float *total_local_crucial_arr;
    if(rank != 0){
        MPI_Send(local_crucial_arr, size - 1, MPI_FLOAT, 0, rank, MPI_COMM_WORLD);
        #ifdef CACULATE_TRANSFER_DATA
            send_count += (size - 1) * sizeof(float);
        #endif
    } else{
        total_local_crucial_arr = (float*)malloc(sizeof(float) * (size * (size - 1)));
        memcpy(total_local_crucial_arr, local_crucial_arr, (size - 1) * sizeof(float));
        #ifdef DEBUG
            for(int i = 0; i < (size - 1); i++){
                printf("%f ", total_local_crucial_arr[i]);
            }
            printf("\n");
        #endif
        for(int i = 1; i < size; i++){
            MPI_Recv(&(total_local_crucial_arr[i * (size - 1)]), size - 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, &status);
            #ifdef CACULATE_TRANSFER_DATA
                recv_count += sizeof(float) * (size - 1);
            #endif
        }
    }
    #ifdef STEPTIME
        step_end[2] = clock();
    #endif
    #ifdef DEBUG
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0){
            for(int i = 0; i < (size - 1) * size; i++){
                printf("%f ", total_local_crucial_arr[i]);
            }
            printf("\n");
        }
    #endif
    
    // 这里我们需要从（size - 1) * size中选出 size - 1 个元素
    #ifdef STEPTIME
        step_start[3] = clock();
    #endif
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
        #ifdef DEBUG
            for(int i = 0; i < size * (size - 1); i++){
                printf("%f ", total_local_crucial_arr_tmp[i]);
            }
            printf("\n");
        #endif
        free(total_local_crucial_arr_tmp);
        free(now_index);
    }
    MPI_Bcast(root_crucial_arr, size - 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    #ifdef STEPTIME
        step_end[3] = clock();
    #endif

    #ifdef CACULATE_TRANSFER_DATA
        if(rank != 0){
            recv_count += sizeof(float) * (size - 1);
        } else{
            send_count += sizeof(float) * (size - 1) * (size - 1);
        }
    #endif
    #ifdef DEBUG
        printf("rank[%d]: ", rank);
        for(int i = 0; i < size - 1; i++){
            printf("%f ", root_crucial_arr[i]);
        }
        printf("\n");
    #endif

    #ifdef STEPTIME
        step_start[4] = clock();
    #endif
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

    #ifdef DEBUG
        printf("rank[%d]: ", rank);
        for(int i = 0; i < size + 1; i++){
            printf("%lu ", last_ele_index[i]);
        }
        printf("\n");
    #endif

    #ifdef STEPTIME
        step_end[4] = clock();
        if(rank == 0){
            step_start[5] = clock();
        }
    #endif

    // 将每个进程i第j段的数量发给j进程，并接收其他进程第i段的数量。
    uint64_t *ele_number = (uint64_t*)malloc(sizeof(uint64_t) * (size));
    for(int i = 0; i < size; i++){
        ele_number[i] = last_ele_index[i+1] - last_ele_index[i];
    }

    #ifdef DEBUG
        printf("rank[%d]: ", rank);
        for(int i = 0; i < size; i++){
            printf("%lu ", ele_number[i]);
        }
        printf("\n");
    #endif
    // 设置一个变量接收分配给进程的数据大小
    uint64_t *proc_ele_num = (uint64_t*)malloc(sizeof(uint64_t) * (size));
    float *proc_data;
    for(int i = 0; i < size; i++){
        if(i != rank){
            MPI_Sendrecv(&(ele_number[i]), 1, MPI_UINT64_T, i, 0,
                &(proc_ele_num[i]), 1, MPI_UINT64_T, i, 0, MPI_COMM_WORLD, &status);
            #ifdef CACULATE_TRANSFER_DATA
                send_count += sizeof(uint64_t) * 1;
                recv_count += sizeof(uint64_t) * 1;
            #endif
        } else {
            proc_ele_num[i] = ele_number[i];
        }
    }

    #ifdef DEBUG
        if(rank == 0){
            for(int i = 0; i < size; i++){
                printf("%lu ", proc_ele_num[i]);
            }
            printf("\n");
        }
    #endif
    uint64_t sum = 0;
    for(int i = 0; i < size; i++){
        sum += proc_ele_num[i];
    }

    #ifdef DEBUG
        printf("rank[%d]: ", rank);
        printf("sum = %lu\n", sum);
    #endif

    #ifdef CHECK
        if(rank != 0){
            MPI_Send(&sum, 1, MPI_UINT64_T, 0, 345, MPI_COMM_WORLD);
        }
    #endif
    
    uint64_t *proc_start_index = malloc(sizeof(uint64_t) * (size));
    proc_start_index[0] = 0;
    for(int i = 1; i < size; i++){
        proc_start_index[i] = proc_start_index[i-1] + proc_ele_num[i-1];
    }

    #ifdef DEBUG
        printf("rank[%d]: ", rank);
        for(int i = 0; i < size; i++){
            printf("%lu ", proc_start_index[i]);
        }
        printf("\n");
    #endif

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
            #ifdef CACULATE_TRANSFER_DATA
                send_count += sizeof(float) * ele_number[i];
                recv_count += sizeof(float) * proc_ele_num[i];
            #endif
        } else{
            memcpy(&(proc_data[proc_start_index[i]]), &(local_arr[last_ele_index[i]]), ele_number[i] * sizeof(float));
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    #ifdef STEPTIME
        if(rank == 0){
            step_end[5] = clock();
        }
        step_start[6] = clock();
    #endif
    //每个进程进行局部归并排序
    float *sort_proc_data = (float*)malloc(sizeof(float) * sum);
    uint64_t *now_index = (uint64_t*)malloc(sizeof(uint64_t) * (size));
    for(int i = 0; i < size; i++){
        now_index[i] = proc_start_index[i];
    }
    int min_index = 0;
    float min_value = 0.;
    uint64_t sort_num = 0;
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
    #ifdef STEPTIME
        step_end[6] = clock();
        if(rank == 0){
            printf("rank[0]:  step 1 time: %f s.\n", (double)(step_end[0] - step_start[0]) / CLOCKS_PER_SEC);
        }
        for(int j = 1; j < 5; j++){
            printf("rank[%d]: step %d time: %f s.\n", rank, j + 1, (double)(step_end[j] - step_start[j]) / CLOCKS_PER_SEC);
        }
        if(rank == 0){
            printf("rank[0]:  step 6 time: %f s.\n", (double)(step_end[5] - step_start[5]) / CLOCKS_PER_SEC);
        }
        printf("rank[%d]: step 7 time: %f s.\n", rank, (double)(step_end[6] - step_start[6]) / CLOCKS_PER_SEC);
        MPI_Barrier(MPI_COMM_WORLD);
    #endif

    #ifdef DEBUG
        printf("rank[%d]: ", rank);
        for(int i = 0; i < 10; i++){
            printf("%f ", sort_proc_data[i]);
        }
        printf("\n");
    #endif
    // 测试结果是否准确
    #ifdef CHECK
        uint64_t *local_sum = (uint64_t*)malloc(sizeof(uint64_t) * size);
        uint64_t *start_index = (uint64_t*)malloc(sizeof(uint64_t) * size);
        // float *sorted_data = (float*)malloc(sizeof(float) * n);
        if(rank == 0){
            for(int i = 1; i < size; i++){
                MPI_Recv(&(local_sum[i]), 1, MPI_UINT64_T, i, 345, MPI_COMM_WORLD, &status);
            }
            local_sum[0] = sum;
            start_index[0] = 0;
            for(int i = 1; i < size; i++){
                start_index[i] = start_index[i - 1] + local_sum[i - 1]; 
            }
            #ifdef DEBUG
                printf("rank[%d]: ", rank);
                for(int i = 0; i < size; i++){
                    printf("%lu ", local_sum[i]);
                }
                printf("\n");
            #endif
            for(int i = 1; i < size; i++){
                MPI_Recv(&(arr[start_index[i]]), local_sum[i], MPI_FLOAT, i, 254, MPI_COMM_WORLD, &status);
            }
            memcpy(arr, sort_proc_data, sum * sizeof(float));

        } else{
            MPI_Send(sort_proc_data, sum, MPI_FLOAT, 0, 254, MPI_COMM_WORLD);
        }
        
        if(rank == 0){
            for(int i = 0; i < size; i++){
                printf("%lu ", start_index[i]);
            }
            printf("\n");
        }
        
        free(local_sum);
        free(start_index);
    #endif

    #ifdef CACULATE_TRANSFER_DATA
        printf("rank[%d]: ", rank);
        printf("send message count: %lu Byte.\n", send_count);
        printf("rank[%d]: ", rank);
        printf("recv message count: %lu Byte.\n", recv_count);
        MPI_Barrier(MPI_COMM_WORLD);
    #endif
    
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