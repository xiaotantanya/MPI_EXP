#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

void matrix_gen(float *a, uint64_t N, float seed);
uint64_t get_Num_from_cmd(int argc, char** argv);
void check_correct(float* sort_arr, char* check_file, uint64_t N);
void gen_correct_file(char* check_file, uint64_t N, int rank, int force, int is_test, float seed);

int compare (const void * a, const void * b)
{
    return ( *(float*)a - *(float*)b > 0 ? 1 :- 1);
}

int main(int argc, char** argv){
    if (argc < 2){
        printf("Usage: %s <number_str>\n", argv[0]);
        return 1;
    }

    uint64_t N = get_Num_from_cmd(argc, argv);
    float seed = 0.344;
    // 异常检测
    if (N == 1){
        return 1;
    }

    float* a = (float*)malloc(N * sizeof(float));
    matrix_gen(a, N, seed);

    // 记录排序前的时间
    clock_t start_time = clock();

    qsort(a, N, sizeof(float), compare);

    // 记录排序后的时间
    clock_t end_time = clock();

    // 计算排序时间
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // 打印排序时间
    printf("qsort took %f seconds to sort %s elements.\n", time_taken, argv[1]);


    char check_file[20] = "./sorted_";
    strcat(check_file, argv[1]);
    printf("check file: %s\n", check_file);
    gen_correct_file(check_file, N, 0, 1, 0, seed);
    check_correct(a, check_file, N);
    // 释放内存
    free(a);
}