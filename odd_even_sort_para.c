#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdint.h>
#include <time.h>

void matrix_gen(float *a, uint64_t N, float seed);
uint64_t get_Num_from_cmd(int argc, char** argv, int rank);

int compare (const void * a, const void * b)
{
    return ( *(float*)a - *(float*)b > 0 ? 1 :- 1);
}


void odd_even_sort_parallel(float *arr, uint64_t n, int rank, int size) {
    uint64_t local_n = n / size;
    // printf("%ld\n", n);
    // printf("%d\n", size);
    // printf("%lu\n", local_n);
    float *local_arr = (float *)malloc(local_n * sizeof(float));
    MPI_Scatter(arr, local_n, MPI_FLOAT, local_arr, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //局部排序(快排)
    qsort(local_arr, local_n, sizeof(float), compare);

    int sorted_local = 0;
    while (!sorted_local) {
        sorted_local = 1;

        // 奇数阶段和偶数阶段切换
        for (int phase = 0; phase < size; phase++) {
            int partner = (phase + rank) % 2 == 0 ? rank + 1 : rank - 1;
            if (partner >= 0 && partner < size) {
                float *recv_arr = (float *)malloc(local_n * sizeof(float));
                MPI_Sendrecv(local_arr, local_n, MPI_FLOAT, partner, 0,
                             recv_arr, local_n, MPI_FLOAT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                float *temp_arr = (float *)malloc(2 * local_n * sizeof(float));
                if (rank < partner) {
                    for (int i = 0; i < local_n; i++) temp_arr[i] = local_arr[i];
                    for (int i = 0; i < local_n; i++) temp_arr[i + local_n] = recv_arr[i];
                    //这里rank进程只需要前local_n个小的，partner只需要后面local_n个大的
                    float* first_point = temp_arr;
                    float* second_point = temp_arr + local_n;
                    uint64_t  local_index = 0;
                    while(local_index < local_n){
                        if (*first_point < *second_point){
                            local_arr[local_index] = *first_point;
                            first_point = first_point + 1;
                        } else{
                            local_arr[local_index] = *second_point;
                            second_point = second_point + 1;
                        }
                        local_index = local_index + 1;
                    }
                    //两个指针，偏移
                    //加一步local_arr赋值操作
                } else {
                    for (int i = 0; i < local_n; i++) temp_arr[i] = recv_arr[i];
                    for (int i = 0; i < local_n; i++) temp_arr[i + local_n] = local_arr[i];
                    //这里rank进程只需要后local_n个大的，partner只需要签名local_n个小的
                    float* first_point = temp_arr +  local_n - 1;
                    float* second_point = temp_arr + 2 * local_n - 1;
                    uint64_t  local_index = local_n - 1;
                    uint64_t local_sum = 0;
                    printf("%lu\n", local_index);
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
                    //两个指针，偏移
                    //合并
                    //加一步local_arr赋值操作
                }


                free(recv_arr);
                free(temp_arr);
            }

            // int global_sorted;
            // MPI_Allreduce(&sorted_local, &global_sorted, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
            // sorted_local = global_sorted;
        }
    }

    MPI_Gather(local_arr, local_n, MPI_FLOAT, arr, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    free(local_arr);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(rank==0){
        if (argc < 2){
            printf("Usage: %s <number_str>\n", argv[0]);
            return 1;
        }
    }
    uint64_t N = get_Num_from_cmd(argc, argv, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    // uint64_t N = 16;
    float *arr = NULL;
    float seed = 0.344;
    clock_t start_time, end_time;
    if (rank == 0) {
        
        arr = (float *)malloc(N* sizeof(float));
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
    MPI_Barrier(MPI_COMM_WORLD);
    odd_even_sort_parallel(arr, N, rank, size);

    if (rank == 0) {
        // printf("Sorted array:\n");
        // for (int i = 0; i < N; i++) printf("%f ", arr[i]);
        // printf("\n");
        end_time = clock();
        // 计算排序时间
        double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        // 打印排序时间
        printf("qsort took %f seconds to sort %s elements.\n", time_taken, argv[1]);
        
        free(arr);
    }

    MPI_Finalize();
    return 0;
}
