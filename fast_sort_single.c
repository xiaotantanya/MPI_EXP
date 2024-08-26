#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

#include "setting.h"
#include "utils.h"

int main(int argc, char** argv){

    if (argc < 2){
        printf("Usage: %s <number_str>\n", argv[0]);
        return 1;
    }

    int is_test = 0;
    uint64_t N = get_Num_from_cmd(argc, argv, 0, is_test);
    float seed = 0.344;
    
    // 异常检测
    if (N == 1){
        return 1;
    }

    float* arr = (float*)malloc(N * sizeof(float));
    matrix_gen(arr, N, seed);
    if(initial_data == 1){
        qsort(arr, N, sizeof(float), re_compare);
        printf("A\n");
    } else if(initial_data == 2){
        qsort(arr, N, sizeof(float), compare);
        printf("B\n");
    } else if(initial_data == 3){
        qsort(arr, N, sizeof(float), compare);
        uint64_t exchange_num  = 0;
        uint64_t exchange_total = (uint64_t)N / (uint64_t)100;
        while(exchange_num < exchange_total){
            uint64_t rand_index = (uint64_t)(((float)rand() / (uint64_t)RAND_MAX) * (N - 1));
            if(arr[rand_index] != arr[rand_index + 1]){
                float temp = arr[rand_index];
                arr[rand_index] = arr[rand_index + 1];
                arr[rand_index + 1] = temp;
                exchange_num += 1;
                // printf("%lu\n", exchange_num);
            }
        }
        printf("C\n");
    }
    // 记录排序前的时间
    clock_t start_time = clock();

    qsort(arr, N, sizeof(float), compare);

    // 记录排序后的时间
    clock_t end_time = clock();

    // 计算排序时间
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // 打印排序时间
    printf("qsort took %f seconds to sort %s elements.\n", time_taken, argv[1]);

    // #ifdef CHECK
    //     char check_file[20] = "./sorted_";
    //     strcat(check_file, argv[1]);
    //     printf("check file: %s\n", check_file);
    //     gen_correct_file(check_file, N, 0, 1, 0, seed);
    //     check_correct(arr, check_file, N);
    // #endif

    // 释放内存
    free(arr);
}
