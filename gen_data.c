#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

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

int compare_local (const void * a, const void * b)
{
    return ( *(float*)a - *(float*)b > 0 ? 1 :- 1);
}

void gen_correct_file(char* check_file, uint64_t N, int rank, int force, int is_test, float seed){
    if(rank == 0){
        if(force || access(check_file, 0)){
            float* arr_data = (float*)malloc(sizeof(float) * N);
            matrix_gen(arr_data, N, seed);
            qsort(arr_data, N, sizeof(float), compare_local);

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
    for(int i = 0; i < 20; i++){
        printf("%f ", read_arr[1024 * 1024 * i]);
    }
    printf("\n");
    int correct = memcmp(sort_arr, read_arr, sizeof(float) * N);
    char* correct_str = (correct == 0 ? "true" : "false");
    printf("rank[%d]: sort correct: %s\n", 0, correct_str);
    free(read_arr);
}