#include <stdint.h>
#include <stdio.h>
#include <string.h>

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

uint64_t get_Num_from_cmd(int argc, char** argv, int rank){
    char * Num_str = argv[1];
    const uint64_t M_num = 1024 * 1024;
    uint64_t N = 0;
    if(strcmp(Num_str, "256M")==0){
        N = 256 * M_num;
    } else if(strcmp(Num_str, "1G")==0){
        N = 1024 * M_num;
    } else if(strcmp(Num_str, "4G")==0){
        N = 4 * 1024 * M_num;
    } else{
        printf("Number Error !\n");
        return 1;
    }
    if (rank ==0) {
        printf("N = %lu\n", N);
    }
    return N;
}