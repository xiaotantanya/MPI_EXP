#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>
#include <string.h>

#include "utils.h"
#include "setting.h"
#include "sort.h"


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
